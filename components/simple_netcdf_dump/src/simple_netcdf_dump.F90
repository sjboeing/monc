module simple_netcdf_dump_mod
 !==============================================================================
 ! simple_netcdf_dump_mod
 ! - Coordinates: gx0/gx1/gy0/gy1 from state%local_grid%start/end(*)
 ! - Per-rank NetCDF files in per-node/per-set directories
 ! - Robust open-or-create logic; coordinates/refs written once per file
 ! The content of this file has been produced with the assistance of
 ! University of Leeds, m365 co-pilot
 ! CRITICAL: time must be written last so tlen reflects completed records
 ! IMPORTANT: file_index is derived from global_slot
 ! DESIGN:
 ! Output is organised into discrete time slots:
 !   global_slot = compute_global_slot(t)
 !   (file_index, r) = split_global_slot(global_slot)
 !
 ! We write at most one record per slot.
 !
 ! INVARIANT:
 !   global_slot uniquely determines (file_index, r)
 !==============================================================================
 use, intrinsic :: iso_fortran_env, only : int32, real64
 use netcdf
 use datadefn_mod, only : default_precision
 use monc_component_mod, only : component_descriptor_type
 use state_mod, only : model_state_type
 use prognostics_mod, only : prognostic_field_type
 use grids_mod, only : x_index, y_index, z_index
 use optionsdatabase_mod, only : options_get_real, options_get_integer, &
   options_get_string, options_get_logical
 use logging_mod, only : log_master_log, log_info, log_warn
 use q_indices_mod, only : q_metadata_type, get_indices_descriptor
 use monc_simple_io_utils_mod, only : &
   ensure_dir, make_setdir, make_rowdir, &
   micro_stagger_sleep, &
   inq_varid_if_present, verify_coords_3d, &
   get_q_cf_metadata, cf_define_3d_var, &
   compute_global_slot, split_global_slot, &
   compute_xy_coords, compute_dual_xy_coords, &
   sanitise_registry_name, &
   def_var_1d, def_var_1d_plain, def_time_1d, put_global_att, &
   lowercase_inplace, &
   c_usleep, join_path, pad3, &
   io_jitter, &
   precreate_rows_global, precreate_sets_for_this_row, global_barrier, &
   io_fatal
 implicit none
 private
 public :: simple_netcdf_dump_get_descriptor

 !==================== Options ====================
 real(default_precision) :: dump3d_interval = 0.0_default_precision
 integer :: dump3d_per_file = 0
 character(len=:), allocatable :: base_filename, dump3d_dir
 real(default_precision) :: dxx = 0.0_default_precision
 real(default_precision) :: dyy = 0.0_default_precision

 ! Micro-stagger
 logical :: opt_micro_stagger = .false.
 real(default_precision) :: micro_stagger_max_us = 200.0_default_precision

 ! Burst
 logical :: burst_enabled = .false.
 real(default_precision) :: burst_start = 0.0_default_precision
 real(default_precision) :: burst_end   = 0.0_default_precision
 real(default_precision) :: burst_interval = 0.0_default_precision
 integer :: dump3d_burst_per_file = 0

 !==================== I/O policy ==================
 real(real64), parameter :: fillv3d = -9.96921d36

 !==================== Sink ========================
 type :: sink_3d_type
   character(len=:), allocatable :: root_dir
   character(len=512) :: final_fname
   integer :: ncid = -1

   ! Dimensions
   integer :: dim_x=-1, dim_y=-1, dim_z=-1, dim_t=-1

   ! Coordinate variables (centers)
   integer :: var_x=-1, var_y=-1, var_z=-1, var_time=-1

   ! NEW: dual coordinate variables
   integer :: var_xn=-1, var_yn=-1, var_zn=-1

   ! Vertical reference fields
   integer :: var_thref=-1, var_rho=-1, var_rhon=-1, var_prefn=-1

   ! Core 3D fields
   integer :: vid_u=-1, vid_v=-1, vid_w=-1, vid_th=-1

   ! Added tke and dissipation 
   integer :: vid_sgs_tke=-1, vid_sgs_dissipation=-1

   ! Scalar tracers
   integer, allocatable :: q_idx(:), q_vid(:)
   character(len=64), allocatable :: q_vnm(:)

   logical :: opened = .false.
   logical :: wrote_coords = .false.

   real(default_precision) :: write_interval = 0.0_default_precision
   real(default_precision) :: window_anchor  = 0.0_default_precision
   integer :: file_index = 0
   integer :: last_slot = -1
   integer :: per_file_count = 0
 end type sink_3d_type

 type(sink_3d_type) :: sink_norm, sink_burst

 !==================== Work buffer =================
 logical, save :: ws3d_init = .false.
 integer, save :: ws_nx = 0, ws_ny = 0, ws_nz = 0
 real(default_precision), allocatable, save :: buf_zyx_ws(:,:,:)

 !==== for directory pre-creation ====
 real(default_precision) :: termination_time = 0.0_default_precision
 real(default_precision) :: tmax_burst      = 0.0_default_precision

contains

 function simple_netcdf_dump_get_descriptor() result(desc)
   type(component_descriptor_type) :: desc
   desc%name = 'simple_netcdf_dump'
   desc%version = 1.1_default_precision
   desc%initialisation => initialise
   desc%timestep       => timestep
   desc%finalisation   => finalise
 end function simple_netcdf_dump_get_descriptor

 subroutine initialise(state)
   type(model_state_type), target, intent(inout) :: state
   character(len=256) :: tmp
   integer :: rx, ry
   character(len=256) :: msg
   real(default_precision) :: normal_file_window, burst_file_window

   dump3d_interval = options_get_real(state%options_database, "dump_3d_interval")
   dump3d_per_file = options_get_integer(state%options_database, "dump_3d_per_file")
   dxx = options_get_real(state%options_database, "dxx")
   dyy = options_get_real(state%options_database, "dyy")
   termination_time = options_get_real(state%options_database, "termination_time")
   if (termination_time < 0.0_default_precision) &
     call io_fatal('3D:initialise: termination_time must be >= 0', -1)

   tmp = options_get_string(state%options_database, "dump_3d_base_filename")
   if (len_trim(tmp) == 0) call io_fatal('3D:initialise: dump_3d_base_filename not set', -1)
   base_filename = trim(tmp)

   tmp = options_get_string(state%options_database, "dump_3d_dir")
   if (len_trim(tmp) == 0) call io_fatal('3D:initialise: dump_3d_dir not set', -1)
   dump3d_dir = trim(tmp)
   call ensure_dir(dump3d_dir)

   opt_micro_stagger     = options_get_logical(state%options_database, "dump_3d_micro_stagger")
   micro_stagger_max_us  = options_get_real   (state%options_database, "dump_3d_micro_stagger_max_us")

   burst_enabled = options_get_logical(state%options_database, "dump_3d_burst_enabled")
   if (burst_enabled) then
     burst_start   = options_get_real   (state%options_database, "dump_3d_burst_start")
     burst_end     = options_get_real   (state%options_database, "dump_3d_burst_end")
     burst_interval= options_get_real   (state%options_database, "dump_3d_burst_interval")
     dump3d_burst_per_file = options_get_integer(state%options_database, "dump_3d_burst_per_file")
     call ensure_dir(trim(dump3d_dir)//"_burst")
   end if

   if (dump3d_interval <= 0.0_default_precision) call io_fatal('3D:initialise: dump_3d_interval must be > 0', -1)
   if (dump3d_per_file <= 0) call io_fatal('3D:initialise: dump_3d_per_file must be > 0', -1)
   if (burst_enabled) then
     if (burst_interval <= 0.0_default_precision) call io_fatal('3D:initialise: dump_3d_burst_interval must be > 0', -1)
     if (dump3d_burst_per_file <= 0) call io_fatal('3D:initialise: dump_3d_burst_per_file must be > 0', -1)
   end if

   normal_file_window = dump3d_per_file * dump3d_interval
   if (burst_enabled) burst_file_window = dump3d_burst_per_file * burst_interval

   ! Normal sink
   sink_norm%root_dir       = dump3d_dir
   sink_norm%write_interval = dump3d_interval
   sink_norm%window_anchor  = 0.0_default_precision
   sink_norm%opened         = .false.
   sink_norm%wrote_coords   = .false.
   sink_norm%file_index     = 0
   sink_norm%last_slot      = -1
   sink_norm%per_file_count = dump3d_per_file

   ! Burst sink
   if (burst_enabled) then
     sink_burst%root_dir       = trim(dump3d_dir)//"_burst"
     sink_burst%write_interval = burst_interval
     sink_burst%window_anchor  = burst_start
     sink_burst%opened         = .false.
     sink_burst%wrote_coords   = .false.
     sink_burst%file_index     = 0
     sink_burst%last_slot      = -1
     sink_burst%per_file_count = dump3d_burst_per_file
   end if

   write(msg,'(A,ES12.5,A,I0,A,ES12.5)') 'simple_netcdf_dump: normal interval= ', &
        dump3d_interval, ', per_file=', dump3d_per_file, ', window=', normal_file_window
   call log_master_log(log_info, trim(msg))
   if (burst_enabled) then
     write(msg,'(A,ES12.5,A,I0,A,ES12.5)') 'burst interval= ', burst_interval, ', per_file=', &
          dump3d_burst_per_file, ', window=', burst_file_window
     call log_master_log(log_info, trim(msg))
   end if

   !=== Row-first precreation with two global barriers ===
   rx = state%parallel%my_coords(x_index)
   ry = state%parallel%my_coords(y_index)

   ! 1) Create row dirs (row leaders only)
   call precreate_rows_global(dump3d_dir, ry, rx)
   if (burst_enabled) call precreate_rows_global(trim(dump3d_dir)//"_burst", ry, rx)
   call global_barrier()

   ! 2) Create setNNNN under row (row leaders only)
   call precreate_sets_for_this_row(dump3d_dir, ry, rx, 0.0_real64, &
        real(normal_file_window,real64), real(termination_time,real64))
   if (burst_enabled) then
     tmax_burst = max(burst_start, burst_end - 1.0e-9_default_precision * burst_interval)
     call precreate_sets_for_this_row(trim(dump3d_dir)//"_burst", ry, rx, &
          real(burst_start,real64), real(burst_file_window,real64), real(tmax_burst,real64))
   end if
   call global_barrier()
 end subroutine initialise

 subroutine finalise(state)
   type(model_state_type), target, intent(inout) :: state
   integer :: ierr
   if (sink_norm%opened) then
     call io_jitter()
     if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
     ierr = nf90_close(sink_norm%ncid)
     if (ierr /= NF90_NOERR) call io_fatal('3D:finalise:close(norm)', ierr)
   end if
   if (burst_enabled .and. sink_burst%opened) then
     call io_jitter()
     if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
     ierr = nf90_close(sink_burst%ncid)
     if (ierr /= NF90_NOERR) call io_fatal('3D:finalise:close(burst)', ierr)
   end if
   if (allocated(buf_zyx_ws)) deallocate(buf_zyx_ws)
   ws3d_init = .false.; ws_nx = 0; ws_ny = 0; ws_nz = 0
 end subroutine finalise

 subroutine timestep(state)
   type(model_state_type), target, intent(inout) :: state
   real(default_precision) :: t
   logical :: active_burst

   t = state%time

   ! NORMAL
   call timestep_one_sink(sink_norm, state, t)

   ! BURST
   active_burst = burst_enabled .and. &
                 (t >= burst_start) .and. (t < burst_end)

   if (active_burst) then
     call timestep_one_sink(sink_burst, state, t)
   endif

 end subroutine timestep

 subroutine timestep_one_sink(sink, state, t)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), intent(inout) :: state

   real(default_precision), intent(in) :: t

   integer :: global_slot, file_index, r
   real(real64) :: anchor64, interval64

   ! ---------------------------------------
   ! Map time -> global slot
   ! ---------------------------------------
   anchor64   = real(sink%window_anchor, real64)
   interval64 = real(sink%write_interval, real64)

   global_slot = compute_global_slot(real(t,real64), anchor64, interval64)

   ! ---------------------------------------
   ! Skip duplicate slot
   ! ---------------------------------------
   if (global_slot == sink%last_slot) return
   sink%last_slot = global_slot

   ! ---------------------------------------
   ! Split slot -> (file_index, r)
   ! ---------------------------------------
   call split_global_slot(global_slot, sink%per_file_count, file_index, r)

   ! ---------------------------------------
   ! Ensure correct file
   ! ---------------------------------------
   if (file_index /= sink%file_index) then
     if (sink%opened) call sink_close(sink)
     sink%file_index = file_index
   end if

   ! ---------------------------------------
   ! Delegate write decision
   ! ---------------------------------------
   call process_sink_write(sink, state, t, r)

 end subroutine timestep_one_sink

 subroutine process_sink_write(sink, state, t, r)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), intent(inout) :: state

   real(default_precision), intent(in) :: t
   integer, intent(in) :: r

   integer :: tlen, ierr

   ! ---------------------------------------
   ! Ensure file open
   ! ---------------------------------------
   call sink_ensure_open(sink, state, base_filename, dxx, dyy, fillv3d)

   ! ---------------------------------------
   ! Query current length
   ! ---------------------------------------
   ierr = nf90_inquire_dimension(sink%ncid, sink%dim_t, len=tlen)
   if (ierr /= NF90_NOERR) call io_fatal('process inquire_dimension', ierr)

   ! ---------------------------------------
   ! STRICT append rule
   ! ---------------------------------------
   if (r /= tlen + 1) return

   ! ---------------------------------------
   ! Perform write
   ! ---------------------------------------
   call sink_write_coords_and_refs(sink, state, dxx, dyy)
   call ensure_ws3d(state)

   call sink_put_xyz_t(sink, sink%vid_u, state%u%data, state, r, plus_scalar=state%ugal)
   call sink_put_xyz_t(sink, sink%vid_v, state%v%data, state, r, plus_scalar=state%vgal)
   call sink_put_xyz_t(sink, sink%vid_w, state%w%data, state, r)
   call sink_put_xyz_t(sink, sink%vid_th, state%th%data, state, r, plus_thref=.true.)
   call sink_put_xyz_t(sink, sink%vid_sgs_tke, state%sgs_tke%data, state, r)
   call sink_put_xyz_t(sink, sink%vid_sgs_dissipation, state%sgs_dissipation%data, state, r)

   if (allocated(sink%q_vid)) call sink_put_all_q(sink, state, r)

   if (opt_micro_stagger) &
     call micro_stagger_sleep(real(micro_stagger_max_us,kind=real64))
   call io_jitter()

   ! ---------------------------------------
   ! Commit time LAST
   ! ---------------------------------------
   ierr = nf90_put_var(sink%ncid, sink%var_time, real(t,real64), start=(/r/))
   if (ierr /= NF90_NOERR) call io_fatal('process: put_time', ierr)

 end subroutine process_sink_write

 subroutine sink_ensure_open(sink, state, base_filename, dxx, dyy, fillv)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), target, intent(inout) :: state
   character(len=*), intent(in) :: base_filename
   real(default_precision), intent(in) :: dxx, dyy
   real(real64), intent(in) :: fillv
   character(len=512) :: setdir, rowdir, fname
   integer :: rx, ry, ierr

   if (sink%opened) return

   rx = state%parallel%my_coords(x_index)
   ry = state%parallel%my_coords(y_index)

   rowdir = make_rowdir(sink%root_dir, ry)
   setdir = make_setdir(rowdir, sink%file_index)
   fname  = join_path(setdir, trim(base_filename)//'.x'//pad3(rx)//'.y'//pad3(ry)//'.nc')
   sink%final_fname = fname

   if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
   call io_jitter()
   ierr = nf90_create(trim(fname), ior(NF90_NOCLOBBER, NF90_64BIT_OFFSET), sink%ncid)
   if (ierr == NF90_EEXIST) ierr = nf90_open(trim(fname), NF90_WRITE, sink%ncid)
   if (ierr /= NF90_NOERR) call io_fatal('3D:sink_ensure_open:create/open('//trim(fname)//')', ierr)

   ierr = nf90_inq_dimid(sink%ncid,'time', sink%dim_t)
   if (ierr /= NF90_NOERR) then
     call sink_define_file(sink, state, fillv)
   else
     call inq_varid_if_present(sink%ncid,"x",sink%var_x)
     call inq_varid_if_present(sink%ncid,"y",sink%var_y)
     call inq_varid_if_present(sink%ncid,"z",sink%var_z)
     call inq_varid_if_present(sink%ncid,"xn",sink%var_xn)
     call inq_varid_if_present(sink%ncid,"yn",sink%var_yn)
     call inq_varid_if_present(sink%ncid,"zn",sink%var_zn)
     call inq_varid_if_present(sink%ncid,"time",sink%var_time)
     call inq_varid_if_present(sink%ncid,"u",sink%vid_u)
     call inq_varid_if_present(sink%ncid,"v",sink%vid_v)
     call inq_varid_if_present(sink%ncid,"w",sink%vid_w)
     call inq_varid_if_present(sink%ncid,"th",sink%vid_th)
     call inq_varid_if_present(sink%ncid,"sgs_tke",sink%vid_sgs_tke)
     call inq_varid_if_present(sink%ncid,"sgs_dissipation",sink%vid_sgs_dissipation)
     call define_all_q_vars(sink, state, .true.)
   end if

   sink%opened = .true.
 end subroutine sink_ensure_open

 subroutine sink_define_file(sink, state, fillv)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), intent(in) :: state
   real(real64), intent(in) :: fillv
   integer :: nx, ny, nz, ierr, dims4(4)
   character(len=64) :: coords_u, coords_v, coords_w, coords_c
   integer :: ix0, ix1, iy0, iy1, iz0, iz1

   ix0 = state%local_grid%local_domain_start_index(x_index)
   ix1 = state%local_grid%local_domain_end_index (x_index)
   iy0 = state%local_grid%local_domain_start_index(y_index)
   iy1 = state%local_grid%local_domain_end_index (y_index)
   iz0 = state%local_grid%local_domain_start_index(z_index)
   iz1 = state%local_grid%local_domain_end_index (z_index)

   nx = ix1 - ix0 + 1
   ny = iy1 - iy0 + 1
   nz = iz1 - iz0 + 1

   ierr = nf90_def_dim(sink%ncid,"x", nx, sink%dim_x)
   if (ierr /= NF90_NOERR) call io_fatal('3D:def_dim x', ierr)
   ierr = nf90_def_dim(sink%ncid,"y", ny, sink%dim_y)
   if (ierr /= NF90_NOERR) call io_fatal('3D:def_dim y', ierr)
   ierr = nf90_def_dim(sink%ncid,"z", nz, sink%dim_z)
   if (ierr /= NF90_NOERR) call io_fatal('3D:def_dim z', ierr)
   ierr = nf90_def_dim(sink%ncid,"time", NF90_UNLIMITED, sink%dim_t)
   if (ierr /= NF90_NOERR) call io_fatal('3D:def_dim time', ierr)

   call def_var_1d(sink%ncid,"x", sink%var_x, sink%dim_x, "projection_x_coordinate","m","X")
   call def_var_1d(sink%ncid,"y", sink%var_y, sink%dim_y, "projection_y_coordinate","m","Y")
   call def_var_1d(sink%ncid,"z", sink%var_z, sink%dim_z, "height","m","Z",positive="up")

   call def_time_1d(sink%ncid, sink%var_time, sink%dim_t)

   call def_var_1d(sink%ncid, "thref", sink%var_thref, sink%dim_z, "air_potential_temperature","K")
   call def_var_1d(sink%ncid, "rho",   sink%var_rho,   sink%dim_z, "air_density","kg m-3")
   call def_var_1d(sink%ncid, "rhon",  sink%var_rhon,  sink%dim_z, "air_density","kg m-3")
   call def_var_1d(sink%ncid, "prefn", sink%var_prefn, sink%dim_z, "air_pressure","Pa")

   call def_var_1d_plain(sink%ncid,"xn", sink%var_xn, sink%dim_x)
   call def_var_1d_plain(sink%ncid,"yn", sink%var_yn, sink%dim_y)
   call def_var_1d_plain(sink%ncid,"zn", sink%var_zn, sink%dim_z)

   dims4 = (/sink%dim_x, sink%dim_y, sink%dim_z, sink%dim_t/)

   coords_u = "time z y xn"
   coords_v = "time z yn x"
   coords_w = "time zn y x"
   coords_c = "time z y x"

   if (allocated(state%u%data)) &
     call cf_define_3d_var(sink%ncid,"u", dims4, sink%vid_u, "eastward_wind","m s-1", coords_u, fillv)
   if (allocated(state%v%data)) &
     call cf_define_3d_var(sink%ncid,"v", dims4, sink%vid_v, "northward_wind","m s-1", coords_v, fillv)
   if (allocated(state%w%data)) &
     call cf_define_3d_var(sink%ncid,"w", dims4, sink%vid_w, "upward_air_velocity","m s-1", coords_w, fillv)
   if (allocated(state%th%data)) &
     call cf_define_3d_var(sink%ncid,"th", dims4, sink%vid_th, "air_potential_temperature","K", coords_c, fillv)
   if (allocated(state%sgs_tke%data)) &
     call cf_define_3d_var(sink%ncid,"sgs_tke", dims4, sink%vid_sgs_tke, "sgs_turbulent_kinetic_energy","m2 s-2", coords_c, fillv)
   if (allocated(state%sgs_dissipation%data)) &
     call cf_define_3d_var(sink%ncid,"sgs_dissipation", dims4, sink%vid_sgs_dissipation, "sgs_dissipation_rate","m2 s-3", coords_c, fillv)

   call define_all_q_vars(sink, state, .false.)

   call put_global_att(sink%ncid,"Conventions","CF-1.8")
   call put_global_att(sink%ncid,"source","MONC LES")
   call put_global_att(sink%ncid,"title","MONC 3D output (CDF-2)")

   if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
   call io_jitter()
   ierr = nf90_enddef(sink%ncid)
   if (ierr /= NF90_NOERR) call io_fatal('3D:enddef', ierr)
   ierr = nf90_set_fill(sink%ncid, NF90_NOFILL, ierr)
   if (ierr /= NF90_NOERR) call io_fatal('3D:set_fill(NOFILL)', ierr)
 end subroutine sink_define_file

 subroutine sink_write_coords_and_refs(sink, state, dxx, dyy)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), intent(in) :: state
   real(default_precision), intent(in) :: dxx, dyy

   integer :: iz0, iz1, ierr
   integer :: ix0, ix1, iy0, iy1
   integer :: gx0, gx1, gy0, gy1
   real(real64), allocatable :: x(:), y(:), xn(:), yn(:)

   if (sink%wrote_coords) return

   gx0 = state%local_grid%start(x_index)
   gx1 = state%local_grid%end(x_index)
   gy0 = state%local_grid%start(y_index)
   gy1 = state%local_grid%end(y_index)
   ix0 = state%local_grid%local_domain_start_index(x_index)
   ix1 = state%local_grid%local_domain_end_index(x_index)
   iy0 = state%local_grid%local_domain_start_index(y_index)
   iy1 = state%local_grid%local_domain_end_index(y_index)
   iz0 = state%local_grid%local_domain_start_index(z_index)
   iz1 = state%local_grid%local_domain_end_index(z_index)
   
   call compute_xy_coords      (gx0,gx1,gy0,gy1, real(dxx,real64),real(dyy,real64), x,  y)
   call compute_dual_xy_coords (gx0,gx1,gy0,gy1, real(dxx,real64),real(dyy,real64), xn, yn)

   ! Write center coords
   ierr = nf90_put_var(sink%ncid, sink%var_x, x)
   if (ierr /= NF90_NOERR) call io_fatal('3D:put x', ierr)

   ierr = nf90_put_var(sink%ncid, sink%var_y, y)
   if (ierr /= NF90_NOERR) call io_fatal('3D:put y', ierr)

   ! Write dual coords
   ierr = nf90_put_var(sink%ncid, sink%var_xn, xn)
   if (ierr /= NF90_NOERR) call io_fatal('3D:put xn', ierr)

   ierr = nf90_put_var(sink%ncid, sink%var_yn, yn)
   if (ierr /= NF90_NOERR) call io_fatal('3D:put yn', ierr)

   ierr = nf90_put_var(sink%ncid, sink%var_zn, &
        state%global_grid%configuration%vertical%zn(iz0:iz1))
   if (ierr /= NF90_NOERR) call io_fatal('3D:put zn', ierr)

   ! Write center z-coords
   ierr = nf90_put_var(sink%ncid, sink%var_z, &
        state%global_grid%configuration%vertical%z(iz0:iz1))
   if (ierr /= NF90_NOERR) call io_fatal('3D:put z', ierr)

   if (sink%var_thref > 0 .and. allocated(state%global_grid%configuration%vertical%thref)) then
       ierr = nf90_put_var(sink%ncid, sink%var_thref, state%global_grid%configuration%vertical%thref(iz0:iz1))
       if (ierr /= NF90_NOERR) call io_fatal('3D:put thref', ierr)
   end if

   if (sink%var_rho > 0 .and. allocated(state%global_grid%configuration%vertical%rho)) then
       ierr = nf90_put_var(sink%ncid, sink%var_rho, state%global_grid%configuration%vertical%rho(iz0:iz1))
       if (ierr /= NF90_NOERR) call io_fatal('3D:put rho', ierr)
   end if

   if (sink%var_rhon > 0 .and. allocated(state%global_grid%configuration%vertical%rhon)) then
       ierr = nf90_put_var(sink%ncid, sink%var_rhon, state%global_grid%configuration%vertical%rhon(iz0:iz1))
       if (ierr /= NF90_NOERR) call io_fatal('3D:put rhon', ierr)
   end if

   if (sink%var_prefn > 0 .and. allocated(state%global_grid%configuration%vertical%prefn)) then
       ierr = nf90_put_var(sink%ncid, sink%var_prefn, state%global_grid%configuration%vertical%prefn(iz0:iz1))
       if (ierr /= NF90_NOERR) call io_fatal('3D:put prefn', ierr)
   end if

   sink%wrote_coords = .true.
 end subroutine sink_write_coords_and_refs

 subroutine sink_put_xyz_t(sink, vid, slab_kji, S, rec, plus_scalar, plus_thref)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), intent(in) :: S
   integer, intent(in) :: vid, rec
   real(default_precision), intent(in) :: slab_kji(:,:,:)
   real(default_precision), intent(in), optional :: plus_scalar
   logical, intent(in), optional :: plus_thref

   integer :: nx, ny, nz
   integer :: ix0, iy0, iz0
   integer :: start4(4), count4(4)
   integer :: ierr

   if (vid <= 0) return

   ix0 = S%local_grid%local_domain_start_index(x_index)
   iy0 = S%local_grid%local_domain_start_index(y_index)
   iz0 = S%local_grid%local_domain_start_index(z_index)
   nx = S%local_grid%local_domain_end_index(x_index) - ix0 + 1
   ny = S%local_grid%local_domain_end_index(y_index) - iy0 + 1
   nz = S%local_grid%local_domain_end_index(z_index) - iz0 + 1

   call pack_kji_to_zyx_add(buf_zyx_ws, slab_kji, S, nx, ny, nz, plus_scalar, plus_thref)

   start4 = (/1,1,1, rec/)
   count4 = (/nx, ny, nz, 1/)
   ierr = nf90_put_var(sink%ncid, vid, buf_zyx_ws, start=start4, count=count4)
   if (ierr /= NF90_NOERR) call io_fatal('3D:put data var', ierr)

 end subroutine sink_put_xyz_t

 subroutine sink_put_all_q(sink, S, rec)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), intent(in) :: S
   integer, intent(in) :: rec
   integer :: i
   if (.not. allocated(sink%q_vid)) return
   do i = 1, size(sink%q_vid)
     if (sink%q_vid(i) > 0 .and. allocated(S%q(sink%q_idx(i))%data)) then
       call sink_put_xyz_t(sink, sink%q_vid(i), S%q(sink%q_idx(i))%data, S, rec, plus_thref=.false.)
     end if
   end do
 end subroutine sink_put_all_q

 subroutine sink_put_time(sink, time_value, rec)
   type(sink_3d_type), intent(inout) :: sink
   real(default_precision), intent(in) :: time_value
   integer, intent(in) :: rec
   integer :: ierr
   real(default_precision) :: tbuf(1)
   if (sink%var_time <= 0) call io_fatal('3D:time varid invalid', -1)
   tbuf(1) = time_value
   ierr = nf90_put_var(sink%ncid, sink%var_time, tbuf, start=(/rec/), count=(/1/))
   if (ierr /= NF90_NOERR) call io_fatal('3D:put time', ierr)
 end subroutine sink_put_time
 
 subroutine sink_close(sink)
   type(sink_3d_type), intent(inout) :: sink
   integer :: ierr
   if (.not. sink%opened) return
   if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
   call io_jitter()
   ierr = nf90_close(sink%ncid)
   if (ierr /= NF90_NOERR) call io_fatal('3D:close', ierr)
   sink%opened = .false.
   sink%wrote_coords = .false.
   if (allocated(sink%q_idx)) deallocate(sink%q_idx, sink%q_vid, sink%q_vnm)
 end subroutine sink_close

 subroutine ensure_ws3d(S)
   type(model_state_type), intent(in) :: S
   integer :: ix0,ix1,iy0,iy1,iz0,iz1,nx,ny,nz

   ix0 = S%local_grid%local_domain_start_index(x_index)
   ix1 = S%local_grid%local_domain_end_index (x_index)
   iy0 = S%local_grid%local_domain_start_index(y_index)
   iy1 = S%local_grid%local_domain_end_index (y_index)
   iz0 = S%local_grid%local_domain_start_index(z_index)
   iz1 = S%local_grid%local_domain_end_index (z_index)

   nx = ix1-ix0+1
   ny = iy1-iy0+1
   nz = iz1-iz0+1

   if (.not. ws3d_init .or. nx /= ws_nx .or. ny /= ws_ny .or. nz /= ws_nz) then
     if (allocated(buf_zyx_ws)) deallocate(buf_zyx_ws)
     allocate(buf_zyx_ws(nx,ny,nz))
     ws_nx = nx
     ws_ny = ny
     ws_nz = nz
     ws3d_init = .true.
   end if
 end subroutine ensure_ws3d

 subroutine pack_kji_to_zyx_add(buf_zyx, slab_kji, S, nx, ny, nz, plus_scalar, plus_thref)
   type(model_state_type), intent(in) :: S
   real(default_precision), intent(out) :: buf_zyx(:,:,:)
   real(default_precision), intent(in) :: slab_kji(:,:,:)
   integer, intent(in) :: nx, ny, nz
   real(default_precision), intent(in), optional :: plus_scalar
   logical, intent(in), optional :: plus_thref

   real(default_precision) :: ps
   logical :: addtr
   integer :: ix0, iy0, iz0, i, j, k

   ps = 0.0_default_precision
   if (present(plus_scalar)) ps = plus_scalar
  
   addtr = .false.
   if (present(plus_thref)) addtr = plus_thref

   ix0 = S%local_grid%local_domain_start_index(x_index)
   iy0 = S%local_grid%local_domain_start_index(y_index)
   iz0 = S%local_grid%local_domain_start_index(z_index)

   if (addtr .and. allocated(S%global_grid%configuration%vertical%thref)) then
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
       buf_zyx(i,j,k) = slab_kji(iz0+k-1, iy0+j-1, ix0+i-1) + ps + &
            S%global_grid%configuration%vertical%thref(iz0+k-1)
     end do
     end do
     end do
   else
     do k = 1, nz
     do j = 1, ny
     do i = 1, nx
       buf_zyx(i,j,k) = slab_kji(iz0+k-1, iy0+j-1, ix0+i-1) + ps
     end do
     end do
     end do
   end if

 end subroutine pack_kji_to_zyx_add

 subroutine define_all_q_vars(sink, state, append_mode)
   type(sink_3d_type), intent(inout) :: sink
   type(model_state_type), intent(in) :: state
   logical, intent(in) :: append_mode
   integer :: i, nq, nalloc, vid, ierr, dims4(4)
   type(q_metadata_type) :: meta
   character(len=64) :: nm, nm_in
   character(len=128) :: stdn, units
   character(len=64) :: coords_q

   if (allocated(sink%q_idx)) deallocate(sink%q_idx, sink%q_vid, sink%q_vnm)
   if (.not. allocated(state%q)) return

   nq = min(state%number_q_fields, size(state%q))
   nalloc = 0
   do vid = 1, nq
     if (allocated(state%q(vid)%data)) nalloc = nalloc + 1
   end do
   if (nalloc == 0) return

   allocate(sink%q_idx(nalloc), sink%q_vid(nalloc), sink%q_vnm(nalloc))

   dims4 = (/sink%dim_x, sink%dim_y, sink%dim_z, sink%dim_t/)
   coords_q = "time z y x"

   i = 0
   do vid = 1, nq
     if (.not. allocated(state%q(vid)%data)) cycle
     meta = get_indices_descriptor(vid)

     nm = ' '
     if (len_trim(meta%name) > 0) then
       nm_in = trim(meta%name)
       call sanitise_registry_name(nm_in, nm)
       call lowercase_inplace(nm)
     else
       call sanitise_registry_name(nm, nm)
       call lowercase_inplace(nm)
     end if

     if (len_trim(nm) == 0) write(nm,'("q",I0)') vid
     nm = nm(1:min(len(nm),64))

     i = i + 1
     sink%q_idx(i) = vid
     sink%q_vnm(i) = nm

     if (append_mode) then
       ierr = nf90_inq_varid(sink%ncid, nm, sink%q_vid(i))
       if (ierr /= NF90_NOERR) sink%q_vid(i) = -1
     else
       call get_q_cf_metadata_safe(meta%name, stdn, units)
       call cf_define_3d_var(sink%ncid, nm, dims4, sink%q_vid(i), &
            trim(stdn), trim(units), coords_q, fillv3d)
     end if
   end do

 contains

   subroutine get_q_cf_metadata_safe(monc_name_in, stdn, units)
     character(len=*), intent(in) :: monc_name_in
     character(len=*), intent(out):: stdn, units
     character(len=128) :: tmp
     tmp = monc_name_in
     call get_q_cf_metadata(trim(tmp), stdn, units)
   end subroutine get_q_cf_metadata_safe
 end subroutine define_all_q_vars

end module simple_netcdf_dump_mod
