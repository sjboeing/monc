module monc_simple_io_utils_mod
 !
 ! Utilities shared by the simple 2‑D and 3‑D writers
 ! The content of this file has been produced with the assistance of
 ! University of Leeds, m365 co-pilot
 !
 use, intrinsic :: iso_fortran_env, only: int64, real64
 use, intrinsic :: iso_c_binding, only: c_int, c_char
 use mpi
 use netcdf
 implicit none
 private

 !================ PUBLIC API ================
 public :: lowercase_inplace, sanitise_registry_name
 public :: join_path, ensure_dir, mkdir_p
 public :: micro_stagger_sleep, io_jitter
 public :: inq_varid_if_present, get_dim_length, get_dim_length_by_name
 public :: verify_coords_2d, verify_coords_3d
 public :: get_q_cf_metadata
 public :: cf_define_2d_var, cf_define_3d_var
 public :: compute_file_index
 public :: compute_xy_coords, compute_dual_xy_coords
 public :: make_setdir, make_rowdir
 public :: def_var_1d, def_var_1d_plain, def_time_1d, put_global_att
 public :: pad3, pad4 
 public :: c_usleep
 public :: io_fatal
 public :: is_world_leader, is_row_leader
 public :: precreate_rows_global, precreate_sets_for_this_row
 public :: global_barrier
 public :: compute_rec_index
 public :: default_time_tolerance

 !============== POSIX mkdir binding =============
 interface
   function c_mkdir(path, mode) bind(C, name="mkdir") result(rc)
     import :: c_int, c_char
     character(c_char), dimension(*), intent(in) :: path
     integer(c_int), value :: mode
     integer(c_int) :: rc
   end function c_mkdir
 end interface

 interface
   subroutine c_usleep(usec) bind(C, name="usleep")
     use iso_c_binding
     integer(c_int), value :: usec
   end subroutine c_usleep
 end interface

contains

 !=================================================
 ! Small string helpers
 !=================================================
 pure function pad3(i) result(s)
   integer, intent(in) :: i
   character(len=3) :: s
   write(s,'(I3.3)') i
 end function pad3

 pure function pad4(i) result(s)
   integer, intent(in) :: i
   character(len=4) :: s
   write(s,'(I4.4)') i
 end function pad4

 subroutine lowercase_inplace(s)
   character(len=*), intent(inout) :: s
   integer :: i, c, shift
   shift = iachar('a') - iachar('A')
   do i = 1, len_trim(s)
     c = iachar(s(i:i))
     if (c >= iachar('A') .and. c <= iachar('Z')) s(i:i) = achar(c + shift)
   end do
 end subroutine lowercase_inplace

 subroutine sanitise_registry_name(src, dst)
   character(len=*), intent(in) :: src
   character(len=*), intent(out) :: dst
   integer :: i, n
   dst = ""
   n = min(len_trim(src), len(dst))
   do i = 1, n
     if (src(i:i) == ' ') then
       dst(i:i) = '_'
     else
       dst(i:i) = src(i:i)
     end if
   end do
 end subroutine sanitise_registry_name

 !=================================================
 ! Simple path + directory utils
 !=================================================
 pure function join_path(a,b) result(p)
   character(len=*), intent(in) :: a,b
   character(len=:), allocatable :: p
   integer :: la
   la = len_trim(a)
   if (la == 0) then
     p = trim(b)
   elseif (a(la:la) == '/') then
     p = trim(a)//trim(b)
   else
     p = trim(a)//'/'//trim(b)
   end if
 end function join_path

 subroutine ensure_dir(path)
   character(len=*), intent(in) :: path
   call mkdir_p(trim(path))
 end subroutine ensure_dir

 !=================================================
 ! Parallel-safe mkdir_p with node-local leader
 !=================================================
 subroutine mkdir_p(fullpath)
   implicit none
   character(len=*), intent(in) :: fullpath
   integer(c_int), parameter :: MODE = int(o'775', kind=c_int)

   character(len=:), allocatable :: prefix
   integer :: i, n, pos
   character(kind=c_char), allocatable :: buf(:)
   integer(c_int) :: rc
   logical :: exists

   if (len_trim(fullpath) == 0) return

   n = len_trim(fullpath)
   prefix = ""

   do i = 1, n
     if (fullpath(i:i) == '/') then
       if (i > 1) then
         prefix = fullpath(1:i-1)
         call make_one(prefix)
       end if
     end if
   end do

   prefix = fullpath(1:n)
   call make_one(prefix)

contains

   subroutine make_one(path)
     character(len=*), intent(in) :: path
     integer :: j, m, attempt
     integer(c_int) :: delay_us
     logical :: exists
     integer :: ios

     if (len_trim(path) == 0) return

     ! Skip if already exists
     inquire(file=trim(path), exist=exists, iostat=ios)
     if (exists .and. ios == 0) return

     ! Build C-string buffer
     m = len_trim(path)
     allocate(buf(m+1))
     do j = 1, m
       buf(j) = char(iachar(path(j:j)), kind=c_char)
     end do
     buf(m+1) = char(0, kind=c_char)

     ! Try to mkdir; tolerate races from other ranks
     do attempt = 1, 3
       rc = c_mkdir(buf, MODE)
       if (rc == 0_c_int) exit

       ! On failure: check whether another process just created it
       inquire(file=trim(path), exist=exists)
       if (exists) exit

       ! Short backoff (max ~1ms)
       delay_us = int(200 * 2**(attempt-1), c_int)
       call c_usleep(delay_us)
     end do

     deallocate(buf)
   end subroutine make_one

 end subroutine mkdir_p

 !=================================================
 ! Filename helpers
 !=================================================
 pure function make_setdir(root, idx) result(s)
   character(len=*), intent(in) :: root
   integer, intent(in) :: idx
   character(len=:), allocatable :: s
   s = join_path(root, "set"//pad4(idx))
 end function make_setdir

 pure function make_rowdir(root, ry) result(rdir)
   character(len=*), intent(in) :: root
   integer,          intent(in) :: ry
   character(len=:), allocatable :: rdir
   rdir = join_path(root, 'row'//pad3(ry))
 end function make_rowdir

 logical function is_world_leader() result(is_leader)
   integer :: rank, ierr
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   is_leader = (rank == 0)
 end function is_world_leader

 logical function is_row_leader(rx) result(is_leader)
   integer, intent(in) :: rx
   is_leader = (rx == 0)
 end function is_row_leader

 ! Row leaders create <root>/rowYYY
 subroutine precreate_rows_global(root_dir, ry, rx)
   character(len=*), intent(in) :: root_dir
   integer, intent(in) :: ry, rx
   if (is_world_leader()) call ensure_dir(root_dir)                ! <-- world leader only
   call global_barrier()                                          ! <-- ensure root exists
   if (.not. is_row_leader(rx)) return
   call ensure_dir(make_rowdir(root_dir, ry))                     ! rowYYY per row leader
 end subroutine precreate_rows_global

 ! Row leaders create <root>/rowYYY/setNNNN for all windows on their row
 subroutine precreate_sets_for_this_row(root_dir, ry, rx, anchor, window, t_max)
   character(len=*), intent(in) :: root_dir
   integer,          intent(in) :: ry, rx
   real(real64),     intent(in) :: anchor, window, t_max
   integer :: i_max, iset
   character(len=:), allocatable :: rowdir, setdir

   if (.not. is_row_leader(rx)) return
   if (window <= 0.0_real64) return

   rowdir = make_rowdir(root_dir, ry)
   call ensure_dir(rowdir)

   i_max = compute_file_index(t_max, anchor, window)
   if (i_max < 0) i_max = 0

   do iset = 0, i_max
     setdir = make_setdir(rowdir, iset)
     call ensure_dir(setdir)
   end do
 end subroutine precreate_sets_for_this_row

 !=================================================
 ! micro_stagger_sleep
 !=================================================
 subroutine micro_stagger_sleep(max_us)
   real(real64), intent(in) :: max_us
   integer :: rank, size, ierr
   integer(c_int) :: delay_us
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   delay_us = int(min(max_us, max(1.0d0, max_us * ( &
                real(rank+1,real64) / max(1.0d0, real(size,real64)) ))), c_int)
   if (delay_us > 0_c_int) call c_usleep(delay_us)
 end subroutine micro_stagger_sleep

 !-------------------------------------------------
 ! io_jitter
 !-------------------------------------------------
 subroutine io_jitter()
   integer :: rank, ierr
   integer(c_int) :: delay_us
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   ! 0–5000 microseconds = 0–5 ms jitter
   delay_us = modulo(rank * 7919, 5000_c_int)
   call c_usleep(delay_us)
 end subroutine io_jitter


 !=================================================
 ! NetCDF helpers
 !=================================================
 subroutine inq_varid_if_present(ncid,name,vid)
   integer, intent(in) :: ncid
   character(len=*), intent(in) :: name
   integer, intent(out) :: vid
   integer :: ierr
   ierr = nf90_inq_varid(ncid, trim(name), vid)
   if (ierr /= NF90_NOERR) vid = -1
 end subroutine inq_varid_if_present

 integer function get_dim_length(ncid,dimid) result(n)
   integer, intent(in) :: ncid, dimid
   integer :: ierr
   n = 0
   ierr = nf90_inquire_dimension(ncid, dimid, len=n)
 end function get_dim_length

 integer function get_dim_length_by_name(ncid,dimname) result(n)
   integer, intent(in) :: ncid
   character(len=*), intent(in) :: dimname
   integer :: dimid, ierr
   n = 0
   ierr = nf90_inq_dimid(ncid, trim(dimname), dimid)
   if (ierr==NF90_NOERR) n = get_dim_length(ncid, dimid)
 end function get_dim_length_by_name


 !=================================================
 ! q metadata mapping
 !=================================================
 subroutine get_q_cf_metadata(monc_name, stdn, units)
   character(len=*), intent(in) :: monc_name
   character(len=*), intent(out) :: stdn, units
   character(len=128) :: tmp
   tmp = monc_name
   call lowercase_inplace(tmp)
   stdn  = "mass_fraction_of_unspecified_scalar_in_air"
   units = "kg kg-1"
   select case (trim(tmp))
   case ("vapour")
     stdn="specific_humidity"
   case ("cloud_liquid_mass")
     stdn="mass_fraction_of_cloud_liquid_water_in_air"
   case ("rain_mass")
     stdn="mass_fraction_of_rain_in_air"
   case ("snow_mass")
     stdn="mass_fraction_of_snow_in_air"
   case ("ice_mass")
     stdn="mass_fraction_of_cloud_ice_in_air"
   case ("graupel_mass")
     stdn="mass_fraction_of_graupel_in_air"
   case ("cloud_liquid_number")
     stdn="number_concentration_of_cloud_liquid_water_particles_in_air"
     units="m-3"
   case ("rain_number")
     stdn="number_concentration_of_rain_particles_in_air"
     units="m-3"
   case ("snow_number")
     stdn="number_concentration_of_snow_particles_in_air"
     units="m-3"
   case ("ice_number")
     stdn="number_concentration_of_cloud_ice_particles_in_air"
     units="m-3"
   case ("graupel_number")
     stdn="number_concentration_of_graupel_particles_in_air"
     units="m-3"
   end select
 end subroutine get_q_cf_metadata

 !=================================================
 ! Coordinate verification (eps tuned)
 !=================================================
 logical function verify_coords_2d(ncid,var_x,var_y,ix0,ix1,iy0,iy1,dxx,dyy)
   integer, intent(in) :: ncid,var_x,var_y
   integer, intent(in) :: ix0,ix1,iy0,iy1
   real(real64), intent(in) :: dxx,dyy
   integer :: nx,ny,ierr
   real(real64) :: x0,xN,y0,yN,x0f,xNf,y0f,yNf,eps
   real(real64) :: buf(1)

   verify_coords_2d = .false.
   eps = 5d-6

   nx = ix1-ix0+1
   ny = iy1-iy0+1

   x0 = (ix0-0.5d0)*dxx; xN = (ix1-0.5d0)*dxx
   y0 = (iy0-0.5d0)*dyy; yN = (iy1-0.5d0)*dyy

   ierr = nf90_get_var(ncid,var_x,buf,start=(/1/), count=(/1/));  if (ierr/=NF90_NOERR) return
   x0f = buf(1)
   ierr = nf90_get_var(ncid,var_x,buf,start=(/nx/),count=(/1/));  if (ierr/=NF90_NOERR) return
   xNf = buf(1)
   ierr = nf90_get_var(ncid,var_y,buf,start=(/1/), count=(/1/));  if (ierr/=NF90_NOERR) return
   y0f = buf(1)
   ierr = nf90_get_var(ncid,var_y,buf,start=(/ny/),count=(/1/));  if (ierr/=NF90_NOERR) return
   yNf = buf(1)

   verify_coords_2d = (abs(x0f-x0)<=eps*max(1d0,abs(x0)) .and. &
                       abs(xNf-xN)<=eps*max(1d0,abs(xN)) .and. &
                       abs(y0f-y0)<=eps*max(1d0,abs(y0)) .and. &
                       abs(yNf-yN)<=eps*max(1d0,abs(yN)))
 end function verify_coords_2d

 logical function verify_coords_3d(ncid,var_x,var_y,var_z,ix0,ix1,iy0,iy1,z0,zN,dxx,dyy)
   integer, intent(in) :: ncid,var_x,var_y,var_z
   integer, intent(in) :: ix0,ix1,iy0,iy1
   real(real64), intent(in) :: z0,zN,dxx,dyy
   real(real64) :: x0,xN,y0,yN,x0f,xNf,y0f,yNf,z0f,zNf,eps
   integer :: nx,ny,nz,ierr
   real(real64) :: buf(1)

   verify_coords_3d = .false.
   eps = 5d-6

   nx = ix1-ix0+1
   ny = iy1-iy0+1

   x0 = (ix0-0.5d0)*dxx; xN = (ix1-0.5d0)*dxx
   y0 = (iy0-0.5d0)*dyy; yN = (iy1-0.5d0)*dyy

   ierr = nf90_get_var(ncid,var_x,buf,start=(/1/), count=(/1/));  if (ierr/=NF90_NOERR) return
   x0f = buf(1)
   ierr = nf90_get_var(ncid,var_x,buf,start=(/nx/),count=(/1/));  if (ierr/=NF90_NOERR) return
   xNf = buf(1)
   ierr = nf90_get_var(ncid,var_y,buf,start=(/1/), count=(/1/));  if (ierr/=NF90_NOERR) return
   y0f = buf(1)
   ierr = nf90_get_var(ncid,var_y,buf,start=(/ny/),count=(/1/));  if (ierr/=NF90_NOERR) return
   yNf = buf(1)

   nz = get_dim_length_by_name(ncid,"z"); if (nz<=0) return
   ierr = nf90_get_var(ncid,var_z,buf,start=(/1/), count=(/1/));  if (ierr/=NF90_NOERR) return
   z0f = buf(1)
   ierr = nf90_get_var(ncid,var_z,buf,start=(/nz/),count=(/1/));  if (ierr/=NF90_NOERR) return
   zNf = buf(1)

   verify_coords_3d = (abs(x0f-x0)<=eps*max(1d0,abs(x0)) .and. &
                       abs(xNf-xN)<=eps*max(1d0,abs(xN)) .and. &
                       abs(y0f-y0)<=eps*max(1d0,abs(y0)) .and. &
                       abs(yNf-yN)<=eps*max(1d0,abs(yN)) .and. &
                       abs(z0f-z0)<=eps*max(1d0,abs(z0)) .and. &
                       abs(zNf-zN)<=eps*max(1d0,abs(zN)))
 end function verify_coords_3d

 !=================================================
 ! NetCDF var helpers
 !=================================================
 subroutine cf_define_2d_var(ncid, nm, dims, vid, stdname, units, fillv)
   integer, intent(in) :: ncid
   character(len=*), intent(in) :: nm, stdname, units
   integer, intent(in) :: dims(3)
   integer, intent(out) :: vid
   real(real64), intent(in) :: fillv
   integer :: ierr

   ierr = nf90_def_var(ncid, nm, NF90_DOUBLE, dims, vid); call chk(ierr,"def "//nm)

   if (len_trim(stdname)>0) then
     ierr = nf90_put_att(ncid,vid,"standard_name",trim(stdname)); call chk(ierr,"att standard_name "//nm)
   end if

   if (len_trim(units)>0) then
     ierr = nf90_put_att(ncid,vid,"units",trim(units));          call chk(ierr,"att units "//nm)
   end if

   ierr = nf90_put_att(ncid,vid,"_FillValue",fillv);             call chk(ierr,"att _FillValue "//nm)
   ierr = nf90_put_att(ncid,vid,"missing_value",fillv);          call chk(ierr,"att missing_value "//nm)
   ierr = nf90_put_att(ncid,vid,"coordinates","time y x");       call chk(ierr,"att coordinates "//nm)
 end subroutine cf_define_2d_var

 subroutine cf_define_3d_var(ncid, nm, dims, vid, stdname, units, coords, fillv)
   integer, intent(in) :: ncid
   character(len=*), intent(in) :: nm, stdname, units, coords
   integer, intent(in) :: dims(4)
   integer, intent(out) :: vid
   real(real64), intent(in) :: fillv
   integer :: ierr

   ierr = nf90_def_var(ncid, nm, NF90_DOUBLE, dims, vid); call chk(ierr,"def "//nm)

   if (len_trim(stdname)>0) then
     ierr = nf90_put_att(ncid,vid,"standard_name",trim(stdname)); call chk(ierr,"att standard_name "//nm)
   end if

   if (len_trim(units)>0) then
     ierr = nf90_put_att(ncid,vid,"units",trim(units));          call chk(ierr,"att units "//nm)
   end if

   ierr = nf90_put_att(ncid,vid,"_FillValue",fillv);             call chk(ierr,"att _FillValue "//nm)
   ierr = nf90_put_att(ncid,vid,"missing_value",fillv);          call chk(ierr,"att missing_value "//nm)
   ierr = nf90_put_att(ncid,vid,"coordinates",trim(coords));     call chk(ierr,"att coordinates "//nm)
 end subroutine cf_define_3d_var

 subroutine def_var_1d(ncid,name,varid,dimid,stdlong,units,axis,is_standard,long_name,positive)
   integer, intent(in) :: ncid, dimid
   character(len=*), intent(in) :: name, stdlong, units
   character(len=*), intent(in), optional :: axis
   logical, intent(in),  optional :: is_standard
   character(len=*), intent(in), optional :: long_name, positive
   integer, intent(out) :: varid
   integer :: ierr
   logical :: stdflag

   ierr = nf90_def_var(ncid,name,NF90_DOUBLE,(/dimid/),varid); call chk(ierr,'def '//trim(name))

   stdflag = .true.; if (present(is_standard)) stdflag = is_standard
   if (stdflag) then
     ierr = nf90_put_att(ncid,varid,"standard_name",trim(stdlong)); call chk(ierr,'att standard_name '//trim(name))
   else
     ierr = nf90_put_att(ncid,varid,"long_name",trim(stdlong));    call chk(ierr,'att long_name '//trim(name))
   end if
 
   ierr = nf90_put_att(ncid,varid,"units",trim(units)); call chk(ierr,'att units '//trim(name))

   if (present(axis)) then
     ierr = nf90_put_att(ncid,varid,"axis",trim(axis)); call chk(ierr,'att axis '//trim(name))
   end if

   if (present(positive)) then
     ierr = nf90_put_att(ncid,varid,"positive",trim(positive)); call chk(ierr,'att positive '//trim(name))
   end if

   if (present(long_name)) then
     ierr = nf90_put_att(ncid,varid,"long_name",trim(long_name)); call chk(ierr,'att long_name2 '//trim(name))
   end if
 end subroutine def_var_1d

 subroutine def_var_1d_plain(ncid,name,varid,dimid)
   integer, intent(in) :: ncid, dimid
   character(len=*), intent(in) :: name
   integer, intent(out) :: varid
   integer :: ierr

   ierr = nf90_def_var(ncid,name,NF90_DOUBLE,(/dimid/),varid); call chk(ierr,'def '//trim(name))
 end subroutine def_var_1d_plain

 subroutine def_time_1d(ncid,varid,dimid)
   integer, intent(in) :: ncid, dimid
   integer, intent(out) :: varid
   integer :: ierr

   ierr = nf90_def_var(ncid,"time",NF90_DOUBLE,(/dimid/),varid); call chk(ierr,'def time')
   ierr = nf90_put_att(ncid,varid,"standard_name","time");       call chk(ierr,'att time standard_name')
   ierr = nf90_put_att(ncid,varid,"units","seconds since 0");    call chk(ierr,'att time units')
   ierr = nf90_put_att(ncid,varid,"axis","T");                   call chk(ierr,'att time axis')
 end subroutine def_time_1d

 subroutine put_global_att(ncid,key,val)
   integer, intent(in) :: ncid
   character(len=*), intent(in):: key,val
   integer :: ierr

   ierr = nf90_put_att(ncid,NF90_GLOBAL,trim(key),trim(val)); call chk(ierr,'global att '//trim(key))
 end subroutine put_global_att

 !=================================================
 ! Scheduling helpers
 !=================================================
 pure function compute_file_index(t,anchor,window) result(idx)
   real(real64), intent(in) :: t,anchor,window
   integer :: idx
   if (window <= 0.0_real64) then
     idx = 0
   else
     idx = floor((t-anchor)/window)
     if (idx < 0) idx = 0
   end if
 end function compute_file_index

 !=================================================
 ! Coordinate constructors
 !=================================================
 subroutine compute_xy_coords(ix0,ix1,iy0,iy1,dxx,dyy,x,y)
   integer,     intent(in) :: ix0,ix1,iy0,iy1
   real(real64),intent(in) :: dxx,dyy
   real(real64),allocatable,intent(out) :: x(:),y(:)
   integer :: nx,ny,i
   nx = ix1-ix0+1
   ny = iy1-iy0+1
   allocate(x(nx),y(ny))
   do i = 1, nx
     x(i) = (ix0+i-1.5d0)*dxx
   end do
   do i = 1, ny
     y(i) = (iy0+i-1.5d0)*dyy
   end do
 end subroutine compute_xy_coords

 subroutine compute_dual_xy_coords(ix0,ix1,iy0,iy1,dxx,dyy,xn,yn)
   integer,     intent(in) :: ix0,ix1,iy0,iy1
   real(real64),intent(in) :: dxx,dyy
   real(real64),allocatable,intent(out) :: xn(:),yn(:)
   integer :: nx,ny,i
   nx = ix1-ix0+1
   ny = iy1-iy0+1
   allocate(xn(nx),yn(ny))
   do i = 1, nx
     xn(i) = (ix0+i-2.0d0)*dxx
   end do
   do i = 1, ny
     yn(i) = (iy0+i-2.0d0)*dyy
   end do
 end subroutine compute_dual_xy_coords

 !=================================================
 ! NetCDF error helper
 !=================================================
 subroutine chk(st,where)
   integer, intent(in) :: st
   character(len=*), intent(in) :: where
   if (st /= NF90_NOERR) call io_fatal(where, st)
 end subroutine chk

 pure function compute_rec_index(t, win_start, interval) result(rec)
   ! r = 1 + floor((t - win_start + eps) / interval), clamped to >= 1
   real(real64), intent(in) :: t, win_start, interval
   integer :: rec
   real(real64) :: rreal, eps

   eps = default_time_tolerance(interval)
   rreal = 1.0d0 + floor((t - win_start + eps) / interval)
   if (rreal < 1.0d0) rreal = 1.0d0
   rec = int(rreal)
 end function compute_rec_index

 pure function default_time_tolerance(interval) result(tol)
   real(real64), intent(in) :: interval
   real(real64) :: tol
   tol = 1.0d-9 * max(1.0d0, interval)
 end function default_time_tolerance

 subroutine io_fatal(where, st)
   character(len=*), intent(in) :: where
   integer,          intent(in) :: st
   integer :: ierr
   logical :: mpi_inited

   write(*,*) 'FATAL NetCDF error at ', trim(where), ': ', trim(nf90_strerror(st))
   call flush(6)

   call MPI_Initialized(mpi_inited, ierr)
   if (mpi_inited) then
     call MPI_Abort(MPI_COMM_WORLD, merge(st,1,st/=0), ierr)
   else
     error stop 'FATAL NetCDF error (MPI not initialized)'
   end if
 end subroutine io_fatal

 subroutine global_barrier()
   integer :: ierr
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
 end subroutine global_barrier

end module monc_simple_io_utils_mod
