module simple_2d_diagnostics_mod
 !==============================================================================
 ! simple_2d_diagnostics_mod (uses start/end grid indices for coords & verify)
 ! - Coordinates: ix0/ix1/iy0/iy1 taken from state%local_grid%start/end(*)
 ! - Append-mode verify uses the same gx*/gy* convention.
 ! - NetCDF calls check ierr to avoid silent failures.
 ! - Physics/diagnostics math unchanged.
 ! The content of this file has been produced with the assistance of
 ! University of Leeds, m365 co-pilot
 !==============================================================================
 use, intrinsic :: iso_fortran_env, only: int32, real64
 use netcdf
 use datadefn_mod, only: default_precision
 use monc_component_mod, only: component_descriptor_type
 use state_mod, only: model_state_type
 use grids_mod, only: x_index, y_index, z_index
 use optionsdatabase_mod, only: options_get_real, options_get_integer, &
     options_get_string, options_get_logical
 use q_indices_mod, only: q_metadata_type, get_indices_descriptor
 use science_constants_mod, only: cp, rlvap
 use monc_simple_io_utils_mod, only: &
   lowercase_inplace, &
   ensure_dir, make_setdir, make_rowdir, &
   micro_stagger_sleep, &
   verify_coords_2d, cf_define_2d_var, &
   compute_file_index, next_file_open_time, compute_xy_coords, &
   def_var_1d, def_time_1d, put_global_att, inq_varid_if_present, &
   c_usleep, join_path, pad3, &
   io_jitter, &
   precreate_rows_global, precreate_sets_for_this_row, global_barrier, &
   compute_rec_index, default_time_tolerance, io_fatal
 implicit none
 private
 public :: simple_2d_diagnostics_get_descriptor

 !==================== Options (global) ====================
 real(default_precision) :: dump2d_interval = 0.0_default_precision
 integer :: dump2d_per_file = 0
 logical :: dump2d_write_t0 = .false.
 character(len=:), allocatable :: base_filename, dump2d_dir
 real(default_precision) :: dxx = 0.0_default_precision
 real(default_precision) :: dyy = 0.0_default_precision
 real(default_precision), parameter :: qlcrit = 1.0e-6_default_precision

 !==================== Fill value ====================
 real(default_precision), parameter :: fillv2d = -9.96921d36

 ! Micro-stagger
 logical :: opt_micro_stagger = .false.
 real(default_precision) :: micro_stagger_max_us = 200.0_default_precision

 ! Burst mode
 logical :: burst_enabled = .false.
 real(default_precision) :: burst_start = 0.0_default_precision
 real(default_precision) :: burst_end   = 0.0_default_precision
 real(default_precision) :: burst_interval = 0.0_default_precision
 integer :: dump2d_burst_per_file = 0

 !==================== Sink definition ====================
 type :: sink_2d_type
    character(len=:), allocatable :: root_dir
    character(len=512) :: final_fname
    integer :: ncid = -1
    integer :: dim_x=-1, dim_y=-1, dim_t=-1
    integer :: var_x=-1, var_y=-1, var_time=-1
    integer :: vid_vwp=-1, vid_lwp=-1
    integer :: vid_wmax=-1, vid_wmin=-1
    integer :: vid_qlmax=-1, vid_hqlmax=-1
    integer :: vid_clbas=-1, vid_cltop=-1
    integer :: vid_senhf=-1, vid_lathf=-1
    integer :: vid_rwp=-1, vid_swp=-1, vid_iwp=-1, vid_gwp=-1
    logical :: opened = .false.
    logical :: wrote_xy = .false.
    logical :: append_mode = .false.
    integer :: rec_index = -1
    logical :: skip_this_write = .false.
    logical :: continuation_checked = .false.
    real(default_precision) :: write_interval = 0.0_default_precision
    real(default_precision) :: file_window = 0.0_default_precision
    real(default_precision) :: window_anchor = 0.0_default_precision
    real(default_precision) :: next_write_time = huge(1.0_default_precision)
    real(default_precision) :: next_file_open_time = huge(1.0_default_precision)
    integer :: file_index = 0
    integer :: per_file_count = 0
    logical :: wrote_t0 = .false.
    logical :: coords_verified_once = .false.
 end type sink_2d_type

 type(sink_2d_type) :: sink_norm, sink_burst

 !==================== q indices ====================
 integer :: iq_qv=-1, iq_ql=-1
 integer :: iq_r=-1, iq_s=-1, iq_i=-1, iq_g=-1

 !==================== workspaces ====================
 logical :: ws_init = .false.
 integer :: ws_nx=0, ws_ny=0
 real(default_precision), allocatable :: vwp_ws(:,:), lwp_ws(:,:), qlmax_ws(:,:), hqlmax_ws(:,:)
 real(default_precision), allocatable :: clbas_ws(:,:), cltop_ws(:,:), wmax_ws(:,:), wmin_ws(:,:)
 real(default_precision), allocatable :: sen_ws(:,:), lat_ws(:,:)
 real(default_precision), allocatable :: rwp_ws(:,:), swp_ws(:,:), iwp_ws(:,:), gwp_ws(:,:)

 !==== vairables needed for directory pre-creation ====
 real(default_precision) :: termination_time = 0.0_default_precision
 real(default_precision) :: tmax_burst      = 0.0_default_precision

contains

 function simple_2d_diagnostics_get_descriptor() result(desc)
   type(component_descriptor_type) :: desc
   desc%name = 'simple_2d_diagnostics'
   desc%version = 2.1_default_precision
   desc%initialisation => initialise
   desc%timestep       => timestep
   desc%finalisation   => finalise
 end function simple_2d_diagnostics_get_descriptor

 subroutine initialise(state)
   type(model_state_type), target, intent(inout) :: state
   character(len=256) :: tmp
   integer :: rx, ry
   real(default_precision) :: normal_file_window, burst_file_window

   dump2d_interval = options_get_real(state%options_database,'dump_2d_interval')
   dump2d_per_file = options_get_integer(state%options_database,'dump_2d_per_file')
   dump2d_write_t0 = options_get_logical(state%options_database,'dump_2d_write_t0')
   dxx = options_get_real(state%options_database,'dxx')
   dyy = options_get_real(state%options_database,'dyy')
   termination_time = options_get_real(state%options_database,'termination_time')
   if (termination_time < 0.0_default_precision) &
     call io_fatal('2D:initialise: termination_time must be >= 0', -1)

   tmp = options_get_string(state%options_database,'dump_2d_base_filename')
   if (len_trim(tmp)==0) call io_fatal('2D:initialise: dump_2d_base_filename not set', -1)
   base_filename = trim(tmp)

   tmp = options_get_string(state%options_database,'dump_2d_dir')
   if (len_trim(tmp)==0) call io_fatal('2D:initialise: dump_2d_dir not set', -1)
   dump2d_dir = trim(tmp)
   call ensure_dir(dump2d_dir)

   opt_micro_stagger     = options_get_logical(state%options_database,'dump_2d_micro_stagger')
   micro_stagger_max_us  = options_get_real   (state%options_database,'dump_2d_micro_stagger_max_us')

   burst_enabled = options_get_logical(state%options_database,'dump_2d_burst_enabled')
   if (burst_enabled) then
     burst_start = options_get_real   (state%options_database,'dump_2d_burst_start')
     burst_end   = options_get_real   (state%options_database,'dump_2d_burst_end')
     burst_interval = options_get_real(state%options_database,'dump_2d_burst_interval')
     dump2d_burst_per_file = options_get_integer(state%options_database,'dump_2d_burst_per_file')
     call ensure_dir(trim(dump2d_dir)//'_burst')
   end if

   if (dump2d_interval <= 0.0_default_precision) call io_fatal('2D:initialise: dump_2d_interval must be > 0', -1)
   if (dump2d_per_file <= 0) call io_fatal('2D:initialise: dump_2d_per_file must be > 0', -1)
   if (burst_enabled) then
     if (burst_interval <= 0.0_default_precision) call io_fatal('2D:initialise: dump_2d_burst_interval must be > 0', -1)
     if (dump2d_burst_per_file <= 0) call io_fatal('2D:initialise: dump_2d_burst_per_file must be > 0', -1)
   end if

   call map_q_indices(state)

   normal_file_window = dump2d_per_file * dump2d_interval
   if (burst_enabled) burst_file_window = dump2d_burst_per_file * burst_interval

   ! Normal sink
   sink_norm%root_dir       = dump2d_dir
   sink_norm%write_interval = dump2d_interval
   sink_norm%file_window    = normal_file_window
   sink_norm%window_anchor  = 0.0_default_precision
   sink_norm%next_file_open_time = next_file_open_time(0, sink_norm%window_anchor, normal_file_window)
   sink_norm%opened         = .false.
   sink_norm%wrote_xy       = .false.
   sink_norm%append_mode    = .false.
   sink_norm%file_index     = 0
   sink_norm%per_file_count = dump2d_per_file
   sink_norm%wrote_t0       = .false.
   sink_norm%continuation_checked = .false.
   sink_norm%next_write_time = merge(0.0_default_precision, sink_norm%write_interval, dump2d_write_t0)

   ! Burst sink
   if (burst_enabled) then
     sink_burst%root_dir       = trim(dump2d_dir)//'_burst'
     sink_burst%write_interval = burst_interval
     sink_burst%file_window    = burst_file_window
     sink_burst%window_anchor  = burst_start
     sink_burst%next_file_open_time = next_file_open_time(0, burst_start, burst_file_window)
     sink_burst%next_write_time = burst_start
     sink_burst%opened         = .false.
     sink_burst%wrote_xy       = .false.
     sink_burst%append_mode    = .false.
     sink_burst%file_index     = 0
     sink_burst%per_file_count = dump2d_burst_per_file
     sink_burst%wrote_t0       = .false.
     sink_burst%continuation_checked = .false.
   end if

   !=== Row-first precreation with two global barriers ===
   rx = state%parallel%my_coords(x_index)
   ry = state%parallel%my_coords(y_index)
 
   call precreate_rows_global(dump2d_dir, ry, rx)
   if (burst_enabled) call precreate_rows_global(trim(dump2d_dir)//'_burst', ry, rx)
   call global_barrier()

   call precreate_sets_for_this_row(dump2d_dir, ry, rx, 0.0_real64, &
        real(normal_file_window,real64), real(termination_time,real64))
   if (burst_enabled) then
     tmax_burst = max(burst_start, burst_end - 1.0e-9_default_precision * burst_interval)
     call precreate_sets_for_this_row(trim(dump2d_dir)//'_burst', ry, rx, &
          real(burst_start,real64), real(burst_file_window,real64), real(tmax_burst,real64))
   end if
   call global_barrier()
 end subroutine initialise

 subroutine finalise(state)
   type(model_state_type), target, intent(inout) :: state
   if (sink_norm%opened) call sink_close(sink_norm)
   if (burst_enabled .and. sink_burst%opened) call sink_close(sink_burst)
   if (allocated(vwp_ws))  deallocate(vwp_ws)
   if (allocated(lwp_ws))  deallocate(lwp_ws)
   if (allocated(qlmax_ws))deallocate(qlmax_ws)
   if (allocated(hqlmax_ws))deallocate(hqlmax_ws)
   if (allocated(clbas_ws))deallocate(clbas_ws)
   if (allocated(cltop_ws))deallocate(cltop_ws)
   if (allocated(wmax_ws)) deallocate(wmax_ws)
   if (allocated(wmin_ws)) deallocate(wmin_ws)
   if (allocated(sen_ws))  deallocate(sen_ws)
   if (allocated(lat_ws))  deallocate(lat_ws)
   if (allocated(rwp_ws))  deallocate(rwp_ws)
   if (allocated(swp_ws))  deallocate(swp_ws)
   if (allocated(iwp_ws))  deallocate(iwp_ws)
   if (allocated(gwp_ws))  deallocate(gwp_ws)
   ws_init = .false.; ws_nx = 0; ws_ny = 0 
 end subroutine finalise
 
subroutine timestep(state)
  use mpi
  type(model_state_type), target, intent(inout) :: state
  real(default_precision) :: t
  integer(int32) :: idx_norm, idx_burst
  logical :: want_write_norm, want_write_burst
  logical :: rollover_local, rollover_all
  integer :: ierr

  ! >>> ADDED FOR SLOT COMPLETION PATCH <<<
  real(real64) :: win_start64, interval64
  integer :: r_now, tlen_local, ierr

  t = state%time

  ! ================= NORMAL =================
  rollover_local = (t >= sink_norm%next_file_open_time)

  call MPI_Allreduce(rollover_local, rollover_all, 1, MPI_LOGICAL, MPI_LOR, &
                     MPI_COMM_WORLD, ierr)

  if (rollover_all) then
    call global_barrier()

    call sink_close_if_open(sink_norm)

    idx_norm = compute_file_index(t, &
         sink_norm%window_anchor, sink_norm%file_window)

    sink_norm%file_index = idx_norm
    sink_norm%next_file_open_time = next_file_open_time( &
         idx_norm, sink_norm%window_anchor, sink_norm%file_window)
  end if

  want_write_norm = (t >= sink_norm%next_write_time) .or. &
                    (dump2d_write_t0 .and. .not. state%continuation_run .and. &
                     .not. sink_norm%wrote_t0 .and. t == 0.0_default_precision)

  ! Ensure one write attempt after restart (existing logic)
  if (state%continuation_run .and. .not. sink_norm%continuation_checked) then
    want_write_norm = .true.
  end if

  ! >>> CRITICAL PATCH: enforce slot completion <<<
  if (state%continuation_run) then

    win_start64 = real(sink_norm%window_anchor + &
                       sink_norm%file_index * sink_norm%file_window, real64)
    interval64  = real(sink_norm%write_interval, real64)

    r_now = compute_rec_index(real(t,real64), win_start64, interval64)

    tlen_local = 0
    if (sink_norm%opened) then
      ierr=nf90_inquire_dimension(sink_norm%ncid, &
                                 sink_norm%dim_t, len=tlen_local)
      if (ierr /= NF90_NOERR) call io_fatal('2D:timestep slot completion', ierr)
    end if

    if (r_now > tlen_local) then
      want_write_norm = .true.
    end if

  end if

  call process_sink_write(sink_norm, state, t, want_write_norm)

  ! ================= BURST ==================
  rollover_local = (t >= sink_burst%next_file_open_time)

  call MPI_Allreduce(rollover_local, rollover_all, 1, MPI_LOGICAL, MPI_LOR, &
                     MPI_COMM_WORLD, ierr)

  if (rollover_all) then
    call global_barrier()

    call sink_close_if_open(sink_burst)

    idx_burst = compute_file_index(t, &
         sink_burst%window_anchor, sink_burst%file_window)

    sink_burst%file_index = idx_burst
    sink_burst%next_file_open_time = next_file_open_time( &
         idx_burst, sink_burst%window_anchor, sink_burst%file_window)
  end if

  want_write_burst = burst_enabled .and. &
                     t >= burst_start .and. t < burst_end .and. &
                     t >= sink_burst%next_write_time

  ! Ensure one write attempt after restart (existing logic)
  if (state%continuation_run .and. .not. sink_burst%continuation_checked) then
    if (burst_enabled .and. t >= burst_start .and. t < burst_end) then
      want_write_burst = .true.
    end if
  end if

  ! >>> CRITICAL PATCH (burst) <<<
  if (state%continuation_run) then

    if (burst_enabled .and. t >= burst_start .and. t < burst_end) then

      win_start64 = real(sink_burst%window_anchor + &
                         sink_burst%file_index * sink_burst%file_window, real64)
      interval64  = real(sink_burst%write_interval, real64)

      r_now = compute_rec_index(real(t,real64), win_start64, interval64)

      tlen_local = 0
      if (sink_burst%opened) then
        ierr=nf90_inquire_dimension(sink_burst%ncid, &
                                   sink_burst%dim_t, len=tlen_local)
        if (ierr /= NF90_NOERR) call io_fatal('2D:timestep slot completion', ierr)
      end if

      if (r_now > tlen_local) then
        want_write_burst = .true.
      end if

    end if

  end if

  call process_sink_write(sink_burst, state, t, want_write_burst)

end subroutine timestep

 subroutine sink_ensure_open(sink, state, base_filename, t, dxx, dyy, fillv)
   type(sink_2d_type), intent(inout) :: sink
   type(model_state_type), target, intent(inout) :: state
   character(len=*), intent(in) :: base_filename
   real(default_precision), intent(in) :: t, dxx, dyy
   real(default_precision), intent(in) :: fillv
   character(len=512) :: setdir, rowdir, fname
   integer :: rx, ry, ierr

   if (sink%opened) return

   rx = state%parallel%my_coords(x_index)
   ry = state%parallel%my_coords(y_index)

   rowdir = make_rowdir(sink%root_dir, ry)
   setdir = make_setdir(rowdir, sink%file_index)
   fname  = join_path(setdir, trim(base_filename)//'.x'//pad3(rx)//'.y'//pad3(ry)//'.nc')
   sink%final_fname = fname

   call io_jitter()
   if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
   ierr = nf90_create(trim(fname), ior(NF90_NOCLOBBER, NF90_64BIT_OFFSET), sink%ncid)
   if (ierr == NF90_EEXIST) ierr = nf90_open(trim(fname), NF90_WRITE, sink%ncid)
   if (ierr /= NF90_NOERR) call io_fatal('2D:sink_ensure_open:create/open('//trim(fname)//')', ierr)

   ierr = nf90_inq_dimid(sink%ncid, 'time', sink%dim_t)
   if (ierr /= NF90_NOERR) then
     call sink_define_file(sink, state, fillv)
   else
     call inq_varid_if_present(sink%ncid,"x",sink%var_x)
     call inq_varid_if_present(sink%ncid,"y",sink%var_y)
     call inq_varid_if_present(sink%ncid,"time",sink%var_time)
     call inq_varid_if_present(sink%ncid,"vwp",sink%vid_vwp)
     call inq_varid_if_present(sink%ncid,"lwp",sink%vid_lwp)
     call inq_varid_if_present(sink%ncid,"wmax",sink%vid_wmax)
     call inq_varid_if_present(sink%ncid,"wmin",sink%vid_wmin)
     call inq_varid_if_present(sink%ncid,"qlmax",sink%vid_qlmax)
     call inq_varid_if_present(sink%ncid,"hqlmax",sink%vid_hqlmax)
     call inq_varid_if_present(sink%ncid,"clbas",sink%vid_clbas)
     call inq_varid_if_present(sink%ncid,"cltop",sink%vid_cltop)
     call inq_varid_if_present(sink%ncid,"senhf",sink%vid_senhf)
     call inq_varid_if_present(sink%ncid,"lathf",sink%vid_lathf)
     call inq_varid_if_present(sink%ncid,"rwp",sink%vid_rwp)
     call inq_varid_if_present(sink%ncid,"swp",sink%vid_swp)
     call inq_varid_if_present(sink%ncid,"iwp",sink%vid_iwp)
     call inq_varid_if_present(sink%ncid,"gwp",sink%vid_gwp)
   end if

   sink%opened = .true.
 end subroutine sink_ensure_open
 
 subroutine sink_define_file(sink, state, fillv)
   type(sink_2d_type), intent(inout) :: sink
   type(model_state_type), intent(in) :: state
   real(default_precision), intent(in) :: fillv
   integer :: nx, ny, ierr, dims3(3)
   integer :: ix0, ix1, iy0, iy1
   logical :: have_r, have_s, have_i, have_g

   ix0 = state%local_grid%local_domain_start_index(x_index)
   ix1 = state%local_grid%local_domain_end_index (x_index)
   iy0 = state%local_grid%local_domain_start_index(y_index)
   iy1 = state%local_grid%local_domain_end_index (y_index)

   nx = ix1 - ix0 + 1
   ny = iy1 - iy0 + 1

   ierr = nf90_def_dim(sink%ncid,'x',nx, sink%dim_x);    if (ierr /= NF90_NOERR) call io_fatal('2D:def_dim x', ierr)
   ierr = nf90_def_dim(sink%ncid,'y',ny, sink%dim_y);    if (ierr /= NF90_NOERR) call io_fatal('2D:def_dim y', ierr)
   ierr = nf90_def_dim(sink%ncid,'time',NF90_UNLIMITED, sink%dim_t)
   if (ierr /= NF90_NOERR) call io_fatal('2D:def_dim time', ierr)

   call def_var_1d(sink%ncid,'x', sink%var_x, sink%dim_x, 'projection_x_coordinate','m','X')
   call def_var_1d(sink%ncid,'y', sink%var_y, sink%dim_y, 'projection_y_coordinate','m','Y')
   call def_time_1d(sink%ncid, sink%var_time, sink%dim_t)

   dims3 = (/sink%dim_x, sink%dim_y, sink%dim_t/)

   call cf_define_2d_var(sink%ncid,'vwp', dims3, sink%vid_vwp, &
        'atmosphere_mass_content_of_water_vapor','kg m-2', fillv)
   call cf_define_2d_var(sink%ncid,'lwp', dims3, sink%vid_lwp, &
        'atmosphere_mass_content_of_cloud_liquid_water','kg m-2', fillv)
   call cf_define_2d_var(sink%ncid,'wmax', dims3, sink%vid_wmax, &
        'upward_air_velocity','m s-1', fillv)
   call cf_define_2d_var(sink%ncid,'wmin', dims3, sink%vid_wmin, &
        'upward_air_velocity','m s-1', fillv)
   call cf_define_2d_var(sink%ncid,'qlmax', dims3, sink%vid_qlmax,'','kg kg-1', fillv)
   call cf_define_2d_var(sink%ncid,'hqlmax',dims3, sink%vid_hqlmax,'height','m', fillv)
   call cf_define_2d_var(sink%ncid,'clbas', dims3, sink%vid_clbas,'height','m', fillv)
   call cf_define_2d_var(sink%ncid,'cltop', dims3, sink%vid_cltop,'height','m', fillv)
   call cf_define_2d_var(sink%ncid,'senhf', dims3, sink%vid_senhf,'surface_upward_sensible_heat_flux','W m-2', fillv)
   call cf_define_2d_var(sink%ncid,'lathf', dims3, sink%vid_lathf,'surface_upward_latent_heat_flux','W m-2', fillv)

   have_r = (iq_r > 0) .and. allocated(state%q(iq_r)%data)
   have_s = (iq_s > 0) .and. allocated(state%q(iq_s)%data)
   have_i = (iq_i > 0) .and. allocated(state%q(iq_i)%data)
   have_g = (iq_g > 0) .and. allocated(state%q(iq_g)%data)

   if (have_r) call cf_define_2d_var(sink%ncid,'rwp',dims3,sink%vid_rwp, &
        'atmosphere_mass_content_of_rain_water','kg m-2', fillv)
   if (have_s) call cf_define_2d_var(sink%ncid,'swp',dims3,sink%vid_swp, &
        'atmosphere_mass_content_of_snow','kg m-2', fillv)
   if (have_i) call cf_define_2d_var(sink%ncid,'iwp',dims3,sink%vid_iwp, &
        'atmosphere_mass_content_of_cloud_ice','kg m-2', fillv)
   if (have_g) call cf_define_2d_var(sink%ncid,'gwp',dims3,sink%vid_gwp, &
        'atmosphere_mass_content_of_graupel','kg m-2', fillv)

   call put_global_att(sink%ncid,'Conventions','CF-1.8')
   call put_global_att(sink%ncid,'source','MONC LES')
   call put_global_att(sink%ncid,'title','MONC 2D diagnostics (CDF-2)')
   call io_jitter()
   if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
   ierr = nf90_enddef(sink%ncid)
   if (ierr /= NF90_NOERR) call io_fatal('2D:enddef', ierr)
   ierr = nf90_set_fill(sink%ncid, NF90_NOFILL, ierr)
   if (ierr /= NF90_NOERR) call io_fatal('2D:set_fill(NOFILL)', ierr)
 end subroutine sink_define_file

 subroutine process_sink_write(sink, state, t, want_write)
   type(sink_2d_type), intent(inout) :: sink
   type(model_state_type), intent(inout) :: state
   real(default_precision), intent(in) :: t
   logical, intent(in) :: want_write
   integer :: rec
   ! NEW locals (minimal)
   integer :: r, tlen, ierr, r_here
   real(real64) :: win_start64, interval64, tol, t_here(1)

   if (.not. want_write) return
   call ensure_ws(state)
   call compute_all_diagnostics(state)

   call sink_ensure_open(sink, state, base_filename, t, dxx, dyy, fillv2d)

   ! Recompute record index for THIS write
   win_start64 = real(sink%window_anchor + sink%file_index * sink%file_window, real64)
   interval64  = real(sink%write_interval, real64)
   r = compute_rec_index(real(t,real64), win_start64, interval64)
   r = min(max(1, r), sink%per_file_count)
   sink%rec_index = r

   ! One-time continuation guard (per simulation per sink)
   sink%skip_this_write = .false.
   if (state%continuation_run .and. .not. sink%continuation_checked) then
     ierr = nf90_inquire_dimension(sink%ncid, sink%dim_t, len=tlen)
     if (ierr == NF90_NOERR .and. r <= tlen) then
       ierr = nf90_get_var(sink%ncid, sink%var_time, t_here, &
                         start=(/r/), count=(/1/))
       if (ierr == NF90_NOERR .and. t_here(1) == t_here(1)) then
         r_here = compute_rec_index(t_here(1), win_start64, interval64)
         if (r_here == r) then
           sink%skip_this_write = .true.
         end if
      end if
     end if
     sink%continuation_checked = .true.
   end if

   if (sink%skip_this_write) then
     if (dump2d_write_t0 .and. .not. state%continuation_run .and. t == 0.0_default_precision) then
       sink%wrote_t0 = .true.
       sink%next_write_time = sink%write_interval
     else
       sink%next_write_time = sink%next_write_time + sink%write_interval
     end if
     return
   end if

   call sink_write_xy(sink, state, dxx, dyy)
   rec = sink%rec_index
   call sink_put_all_fields(sink, rec)
   call io_jitter()
   if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
   call sink_put_time(sink, t, rec)
 
   if (dump2d_write_t0 .and. t == 0.0_default_precision) then
     sink%wrote_t0 = .true.
     sink%next_write_time = sink%write_interval
   else
     sink%next_write_time = sink%next_write_time + sink%write_interval
   end if
 end subroutine process_sink_write

 subroutine sink_write_xy(sink, state, dxx, dyy)
   type(sink_2d_type), intent(inout) :: sink
   type(model_state_type), intent(in) :: state
   real(default_precision), intent(in) :: dxx, dyy
   integer :: gx0, gx1, gy0, gy1, ierr
   real(real64), allocatable :: x(:), y(:)
   logical :: ok

   if (sink%wrote_xy) return

   gx0 = state%local_grid%start(x_index)
   gx1 = state%local_grid%end(x_index)
   gy0 = state%local_grid%start(y_index)
   gy1 = state%local_grid%end(y_index)

   call compute_xy_coords(gx0, gx1, gy0, gy1, real(dxx,real64), real(dyy,real64), x, y)
   ierr = nf90_put_var(sink%ncid, sink%var_x, x); if (ierr /= NF90_NOERR) call io_fatal('2D:put x', ierr)
   ierr = nf90_put_var(sink%ncid, sink%var_y, y); if (ierr /= NF90_NOERR) call io_fatal('2D:put y', ierr)
   sink%wrote_xy = .true.
 end subroutine sink_write_xy

 subroutine sink_put_all_fields(sink, rec)
   type(sink_2d_type), intent(inout) :: sink
   integer, intent(in) :: rec
   call sink_put_xy_t(sink, sink%vid_vwp, vwp_ws, rec)
   call sink_put_xy_t(sink, sink%vid_lwp, lwp_ws, rec)
   call sink_put_xy_t(sink, sink%vid_qlmax, qlmax_ws, rec)
   call sink_put_xy_t(sink, sink%vid_hqlmax, hqlmax_ws, rec)
   call sink_put_xy_t(sink, sink%vid_clbas, clbas_ws, rec)
   call sink_put_xy_t(sink, sink%vid_cltop, cltop_ws, rec)
   call sink_put_xy_t(sink, sink%vid_wmax, wmax_ws, rec)
   call sink_put_xy_t(sink, sink%vid_wmin, wmin_ws, rec)
   call sink_put_xy_t(sink, sink%vid_senhf, sen_ws,  rec)
   call sink_put_xy_t(sink, sink%vid_lathf, lat_ws,  rec)
   if (sink%vid_rwp>0) call sink_put_xy_t(sink, sink%vid_rwp, rwp_ws, rec)
   if (sink%vid_swp>0) call sink_put_xy_t(sink, sink%vid_swp, swp_ws, rec)
   if (sink%vid_iwp>0) call sink_put_xy_t(sink, sink%vid_iwp, iwp_ws, rec)
   if (sink%vid_gwp>0) call sink_put_xy_t(sink, sink%vid_gwp, gwp_ws, rec)
 end subroutine sink_put_all_fields

 subroutine sink_put_xy_t(sink, vid, a2, rec)
   type(sink_2d_type), intent(inout) :: sink
   integer, intent(in) :: vid, rec
   real(default_precision), intent(in) :: a2(:,:)
   integer :: ierr

   if (vid <= 0) return

   ierr = nf90_put_var(sink%ncid, vid, a2, &
        start=(/1,1,rec/), count=(/size(a2,1), size(a2,2), 1/))
   if (ierr /= NF90_NOERR) call io_fatal('2D:put data var', ierr)
 end subroutine sink_put_xy_t

 subroutine sink_put_time(sink, time_value, rec)
   type(sink_2d_type), intent(inout) :: sink
   real(default_precision), intent(in) :: time_value
   integer, intent(in) :: rec
   integer :: ierr
   real(default_precision) :: tbuf(1)
   if (sink%var_time <= 0) call io_fatal('2D:time varid invalid', -1)
   tbuf(1) = time_value
   ierr = nf90_put_var(sink%ncid, sink%var_time, tbuf, start=(/rec/), count=(/1/))
   if (ierr /= NF90_NOERR) call io_fatal('2D:put time', ierr)
 end subroutine sink_put_time

 subroutine sink_close_if_open(sink)
   type(sink_2d_type), intent(inout) :: sink
   if (sink%opened) call sink_close(sink)
 end subroutine sink_close_if_open

 subroutine sink_close(sink)
   type(sink_2d_type), intent(inout) :: sink
   integer :: ierr
   if (.not. sink%opened) return
   call io_jitter()
   if (opt_micro_stagger) call micro_stagger_sleep(real(micro_stagger_max_us, kind=real64))
   ierr = nf90_close(sink%ncid)
   if (ierr /= NF90_NOERR) call io_fatal('2D:close', ierr)
   sink%opened = .false.
   sink%wrote_xy = .false.
   sink%append_mode= .false.
 end subroutine sink_close

 subroutine compute_all_diagnostics(S)
   type(model_state_type), intent(in) :: S
   call compute_vwp(S)
   call compute_ql_family(S)
   call compute_w_extrema(S)
   call compute_latent(S)
   call compute_sensible(S)
   if (iq_r>0 .and. allocated(S%q(iq_r)%data)) call compute_qx_to_wpath(S, iq_r, rwp_ws)
   if (iq_s>0 .and. allocated(S%q(iq_s)%data)) call compute_qx_to_wpath(S, iq_s, swp_ws)
   if (iq_i>0 .and. allocated(S%q(iq_i)%data)) call compute_qx_to_wpath(S, iq_i, iwp_ws)
   if (iq_g>0 .and. allocated(S%q(iq_g)%data)) call compute_qx_to_wpath(S, iq_g, gwp_ws)
 end subroutine compute_all_diagnostics

subroutine compute_vwp(S)
   use grids_mod, only : x_index, y_index, z_index
   type(model_state_type), intent(in) :: S
   integer :: ix0,ix1,iy0,iy1,iz0,iz1,nx,ny,nz
   integer :: i,j,k, ii,jj,kk, k0,k1
   logical :: have_qv
   real(default_precision), allocatable :: dz(:), rho(:)
   real(default_precision) :: sumv

   if (iq_qv < 1 .or. iq_qv > size(S%q)) return
   have_qv = allocated(S%q(iq_qv)%data)
   if (.not. have_qv) return

   ix0=S%local_grid%local_domain_start_index(x_index)
   ix1=S%local_grid%local_domain_end_index(x_index)
   iy0=S%local_grid%local_domain_start_index(y_index)
   iy1=S%local_grid%local_domain_end_index(y_index)
   iz0=S%local_grid%local_domain_start_index(z_index)
   iz1=S%local_grid%local_domain_end_index(z_index)

   nx=ix1-ix0+1
   ny=iy1-iy0+1
   nz=iz1-iz0+1

   k0=iz0
   k1=iz1

   allocate(dz(nz), rho(nz))
   do k = 1, nz
    if (k < nz) then
      dz(k) = S%global_grid%configuration%vertical%zn(k0+k) - S%global_grid%configuration%vertical%zn(k0+k-1)
    else
      dz(k) = S%global_grid%configuration%vertical%zn(k0+k-1) - S%global_grid%configuration%vertical%zn(k0+k-2)
    end if
   end do

   rho = S%global_grid%configuration%vertical%rho(k0:k1)
   vwp_ws = 0.0_default_precision

   do j = 1, ny
    jj=iy0+j-1
    do i = 1, nx
     ii=ix0+i-1
     sumv = 0.0_default_precision
     do k = 1, nz
      kk=iz0+k-1
      sumv = sumv + rho(k) * S%q(iq_qv)%data(kk,jj,ii) * dz(k)
     end do
     vwp_ws(i,j) = sumv
    end do
   end do

   deallocate(dz, rho)
 end subroutine compute_vwp

 subroutine compute_ql_family(S)
   use grids_mod, only : x_index, y_index, z_index
   type(model_state_type), intent(in) :: S
   integer :: ix0,ix1,iy0,iy1,iz0,iz1,nx,ny,nz
   integer :: i,j,k, ii,jj,kk,kbase,ktop,k0,k1
   logical :: have_ql
   real(default_precision), allocatable :: dz(:), rho(:), z(:)
   real(default_precision) :: suml, qmax, qlk

   if (iq_ql < 1 .or. iq_ql > size(S%q)) return
   have_ql = allocated(S%q(iq_ql)%data)
   if (.not. have_ql) return

   ix0=S%local_grid%local_domain_start_index(x_index)
   ix1=S%local_grid%local_domain_end_index(x_index)
   iy0=S%local_grid%local_domain_start_index(y_index)
   iy1=S%local_grid%local_domain_end_index(y_index)
   iz0=S%local_grid%local_domain_start_index(z_index)
   iz1=S%local_grid%local_domain_end_index(z_index)

   nx=ix1-ix0+1
   ny=iy1-iy0+1
   nz=iz1-iz0+1
   k0=iz0
   k1=iz1

   allocate(dz(nz), rho(nz), z(nz))

   do k = 1, nz
    if (k < nz) then
      dz(k) = S%global_grid%configuration%vertical%zn(k0+k) - S%global_grid%configuration%vertical%zn(k0+k-1)
    else
      dz(k) = S%global_grid%configuration%vertical%zn(k0+k-1) - S%global_grid%configuration%vertical%zn(k0+k-2)
    end if
   end do

   rho = S%global_grid%configuration%vertical%rho(k0:k1)
   z = S%global_grid%configuration%vertical%z(k0:k1)

   lwp_ws = 0.0_default_precision
   qlmax_ws = 0.0_default_precision
   hqlmax_ws = 0.0_default_precision
   clbas_ws = 0.0_default_precision
   cltop_ws = 0.0_default_precision

   do j = 1, ny
    jj=iy0+j-1
    do i = 1, nx
     ii=ix0+i-1

     suml = 0.0_default_precision
     qmax = 0.0_default_precision
     kbase = 0
     ktop = 0

     do k = 1, nz
      kk=iz0+k-1
      qlk = S%q(iq_ql)%data(kk,jj,ii)
      suml = suml + rho(k) * qlk * dz(k)
      if (qlk > qmax) then
        qmax = qlk
        hqlmax_ws(i,j) = z(k)
      end if
      if (kbase == 0 .and. qlk > qlcrit) kbase = k
      if (qlk > qlcrit) ktop = k
     end do

     lwp_ws(i,j) = suml
     qlmax_ws(i,j) = qmax

     if (kbase > 0) then
       clbas_ws(i,j) = z(kbase)
       cltop_ws(i,j) = z(ktop)
     end if

    end do
   end do

   deallocate(dz, rho, z)
 end subroutine compute_ql_family

 subroutine compute_w_extrema(S)
   use grids_mod, only : x_index, y_index, z_index
   type(model_state_type), intent(in) :: S
   integer :: ix0,ix1,iy0,iy1,iz0,iz1,nx,ny,nz
   integer :: i,j,k, ii,jj,kk
   real(default_precision) :: val, wmin_ij, wmax_ij

   if (.not. allocated(S%w%data)) return

   ix0=S%local_grid%local_domain_start_index(x_index)
   ix1=S%local_grid%local_domain_end_index(x_index)
   iy0=S%local_grid%local_domain_start_index(y_index)
   iy1=S%local_grid%local_domain_end_index(y_index)
   iz0=S%local_grid%local_domain_start_index(z_index)
   iz1=S%local_grid%local_domain_end_index(z_index)

   nx=ix1-ix0+1
   ny=iy1-iy0+1
   nz=iz1-iz0+1

   do j = 1, ny
    jj=iy0+j-1
    do i = 1, nx
     ii=ix0+i-1
     wmin_ij = huge(1.0_default_precision)
     wmax_ij = -huge(1.0_default_precision)
     do k = 1, nz
      kk=iz0+k-1
      val = S%w%data(kk,jj,ii)
      if (val < wmin_ij) wmin_ij = val
      if (val > wmax_ij) wmax_ij = val
     end do
     wmin_ws(i,j) = wmin_ij
     wmax_ws(i,j) = wmax_ij
    end do
   end do
 end subroutine compute_w_extrema

 subroutine compute_latent(S)
   use grids_mod, only : x_index, y_index
   use science_constants_mod, only : rlvap
   type(model_state_type), intent(in) :: S
   integer :: ix0,ix1,iy0,iy1,nx,ny,i,j,ii,jj,iqv
   real(default_precision) :: rdzn2, rhon1, K1, q1, q2

   if (.not. S%use_surface_boundary_conditions) return
   iqv = S%water_vapour_mixing_ratio_index
   if (iqv < 1 .or. .not. allocated(S%q(iqv)%data)) return
   if (.not. allocated(S%diff_coefficient%data)) return

   ix0=S%local_grid%local_domain_start_index(x_index)
   ix1=S%local_grid%local_domain_end_index(x_index)
   iy0=S%local_grid%local_domain_start_index(y_index)
   iy1=S%local_grid%local_domain_end_index(y_index)

   nx=ix1-ix0+1
   ny=iy1-iy0+1

   rdzn2 = S%global_grid%configuration%vertical%rdzn(2)
   rhon1 = S%global_grid%configuration%vertical%rhon(1)

   lat_ws = 0.0_default_precision

   do j = 1, ny
    jj=iy0+j-1
    do i = 1, nx
     ii=ix0+i-1
     K1 = S%diff_coefficient%data(1,jj,ii)
     q1 = S%q(iqv)%data(1,jj,ii)
     q2 = S%q(iqv)%data(2,jj,ii)
     lat_ws(i,j) = K1 * rdzn2 * (q1 - q2) * rlvap * rhon1
    end do
   end do
 end subroutine compute_latent

 subroutine compute_sensible(S)
   use grids_mod, only : x_index, y_index
   use science_constants_mod, only : cp
   type(model_state_type), intent(in) :: S
   integer :: ix0,ix1,iy0,iy1,nx,ny,i,j,ii,jj
   real(default_precision) :: rdzn2, rhon1, dth1, K1, th1, th2

   if (.not. S%use_surface_boundary_conditions) return
   if (.not. allocated(S%diff_coefficient%data)) return
   if (.not. S%th%active .or. .not. allocated(S%th%data)) return

   ix0=S%local_grid%local_domain_start_index(x_index)
   ix1=S%local_grid%local_domain_end_index(x_index)
   iy0=S%local_grid%local_domain_start_index(y_index)
   iy1=S%local_grid%local_domain_end_index(y_index)

   nx=ix1-ix0+1
   ny=iy1-iy0+1

   rdzn2 = S%global_grid%configuration%vertical%rdzn(2)
   rhon1 = S%global_grid%configuration%vertical%rhon(1)
   dth1 = S%global_grid%configuration%vertical%dthref(1)

   sen_ws = 0.0_default_precision

   do j = 1, ny
    jj=iy0+j-1
    do i = 1, nx
     ii=ix0+i-1
     K1 = S%diff_coefficient%data(1,jj,ii)
     th1 = S%th%data(1,jj,ii)
     th2 = S%th%data(2,jj,ii)
     sen_ws(i,j) = K1 * rdzn2 * (th1 - th2 - dth1) * rhon1 * cp
    end do
   end do
 end subroutine compute_sensible

 subroutine compute_qx_to_wpath(S, iq, ws)
   use grids_mod, only : x_index, y_index, z_index
   type(model_state_type), intent(in) :: S
   integer, intent(in) :: iq
   real(default_precision), intent(inout) :: ws(:,:)
   integer :: ix0,ix1,iy0,iy1,iz0,iz1,nx,ny,nz,i,j,k,ii,jj,kk,k0,k1
   logical :: have_q
   real(default_precision), allocatable :: dz(:), rho(:)
   real(default_precision) :: sumx

   if (iq < 1 .or. iq > size(S%q)) return
   have_q = allocated(S%q(iq)%data)
   if (.not. have_q) return

   ix0=S%local_grid%local_domain_start_index(x_index)
   ix1=S%local_grid%local_domain_end_index(x_index)
   iy0=S%local_grid%local_domain_start_index(y_index)
   iy1=S%local_grid%local_domain_end_index(y_index)
   iz0=S%local_grid%local_domain_start_index(z_index)
   iz1=S%local_grid%local_domain_end_index(z_index)

   nx=ix1-ix0+1
   ny=iy1-iy0+1
   nz=iz1-iz0+1

   k0=iz0
   k1=iz1

   allocate(dz(nz), rho(nz))

   do k = 1, nz
    if (k < nz) then
      dz(k) = S%global_grid%configuration%vertical%zn(k0+k) - S%global_grid%configuration%vertical%zn(k0+k-1)
    else
      dz(k) = S%global_grid%configuration%vertical%zn(k0+k-1) - S%global_grid%configuration%vertical%zn(k0+k-2)
    end if
   end do

   rho = S%global_grid%configuration%vertical%rho(k0:k1)

   ws = 0.0_default_precision
   do j = 1, ny
    jj=iy0+j-1
    do i = 1, nx
     ii=ix0+i-1
     sumx = 0.0_default_precision
     do k = 1, nz
      kk=iz0+k-1
      sumx = sumx + rho(k) * S%q(iq)%data(kk,jj,ii) * dz(k)
     end do
     ws(i,j) = sumx
    end do
   end do

   deallocate(dz, rho)
 end subroutine compute_qx_to_wpath

 subroutine ensure_ws(S)
   type(model_state_type), intent(in) :: S
   integer :: ix0,ix1,iy0,iy1,iz0,iz1,nx,ny,nz

   ix0=S%local_grid%local_domain_start_index(x_index)
   ix1=S%local_grid%local_domain_end_index (x_index)
   iy0=S%local_grid%local_domain_start_index(y_index)
   iy1=S%local_grid%local_domain_end_index (y_index)
   iz0=S%local_grid%local_domain_start_index(z_index)
   iz1=S%local_grid%local_domain_end_index (z_index)

   nx=ix1-ix0+1
   ny=iy1-iy0+1
   nz=iz1-iz0+1

   if (.not. ws_init .or. nx/=ws_nx .or. ny/=ws_ny) then
     if (allocated(vwp_ws)) deallocate(vwp_ws)
     if (allocated(lwp_ws)) deallocate(lwp_ws)
     if (allocated(qlmax_ws)) deallocate(qlmax_ws)
     if (allocated(hqlmax_ws)) deallocate(hqlmax_ws)
     if (allocated(clbas_ws)) deallocate(clbas_ws)
     if (allocated(cltop_ws)) deallocate(cltop_ws)
     if (allocated(wmax_ws)) deallocate(wmax_ws)
     if (allocated(wmin_ws)) deallocate(wmin_ws)
     if (allocated(sen_ws)) deallocate(sen_ws)
     if (allocated(lat_ws)) deallocate(lat_ws)
     if (allocated(rwp_ws)) deallocate(rwp_ws)
     if (allocated(swp_ws)) deallocate(swp_ws)
     if (allocated(iwp_ws)) deallocate(iwp_ws)
     if (allocated(gwp_ws)) deallocate(gwp_ws)

     allocate(vwp_ws(nx,ny), lwp_ws(nx,ny), qlmax_ws(nx,ny), hqlmax_ws(nx,ny))
     allocate(clbas_ws(nx,ny), cltop_ws(nx,ny))
     allocate(wmax_ws(nx,ny), wmin_ws(nx,ny))
     allocate(sen_ws(nx,ny),  lat_ws(nx,ny))
     allocate(rwp_ws(nx,ny),  swp_ws(nx,ny),  iwp_ws(nx,ny),  gwp_ws(nx,ny))

     ws_nx=nx
     ws_ny=ny
     ws_init=.true.
   end if
 end subroutine ensure_ws

 subroutine map_q_indices(state)
   type(model_state_type), intent(in) :: state
   integer :: i, n
   type(q_metadata_type) :: meta
   character(len=256) :: nm

   iq_qv=-1; iq_ql=-1; iq_r=-1; iq_s=-1; iq_i=-1; iq_g=-1
   n = state%number_q_fields
   do i=1,n
     meta = get_indices_descriptor(i)
     nm = meta%name
     call lowercase_inplace(nm)
     select case(trim(nm))
     case('vapour');            iq_qv=i
     case('cloud_liquid_mass'); iq_ql=i
     case('rain_mass');         iq_r=i
     case('snow_mass');         iq_s=i
     case('ice_mass');          iq_i=i
     case('graupel_mass');      iq_g=i
     end select
   end do
 end subroutine map_q_indices

end module simple_2d_diagnostics_mod
