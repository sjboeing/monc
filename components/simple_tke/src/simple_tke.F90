module simple_tke_mod
 ! The content of this file has been produced with the assistance of
 ! University of Leeds, m365 co-pilot
   use monc_component_mod,  only : component_descriptor_type
   use state_mod,           only : model_state_type, FORWARD_STEPPING
   use smagorinsky_mod,     only : &
        calculate_half_squared_strain_rate, &
        calculate_richardson_number
   use grids_mod,           only : z_index
   use datadefn_mod,        only : default_precision
   implicit none
   private

   public :: simple_tke_get_descriptor

   integer, parameter :: dp = default_precision

contains

!====================================================================
function simple_tke_get_descriptor() result(desc)
   type(component_descriptor_type) :: desc
   desc%name           = "simple_tke"
   desc%version        = 4.0_dp
   desc%initialisation => initialise
   desc%timestep       => timestep
   desc%finalisation   => finalise
end function simple_tke_get_descriptor

!====================================================================
subroutine initialise(current_state)
   type(model_state_type), target, intent(inout) :: current_state
   ! Nothing to do: fields allocated elsewhere
end subroutine initialise

!====================================================================
subroutine timestep(current_state)
   type(model_state_type), target, intent(inout) :: current_state

   integer :: nz, k
   integer :: ix, iy
   real(dp), allocatable :: ssq_col(:), Ri_col(:)
   real(dp) :: nut, kappa
   real(dp) :: diss, L2

   nz = current_state%local_grid%size(z_index)

   ix = current_state%column_local_x
   iy = current_state%column_local_y

   allocate(ssq_col(nz), Ri_col(nz))

   ! Note we only apply to interior, no halo copies
   ! That should be OK I think (as no gradients of output needed in other modules)
   if (.not. current_state%halo_column) then
      if (current_state%field_stepping == FORWARD_STEPPING) then
         ssq_col=calculate_half_squared_strain_rate(current_state, current_state%u, current_state%v, current_state%w)
         Ri_col=calculate_richardson_number(current_state, ssq_col, current_state%th, current_state%q)
      else
         ssq_col=calculate_half_squared_strain_rate(current_state, current_state%zu, current_state%zv, current_state%zw)
         Ri_col=calculate_richardson_number(current_state, ssq_col, current_state%zth, current_state%zq)
      end if
   end if

   do k = 1, nz

      nut   = current_state%vis_coefficient%data(k,iy,ix)
      kappa = current_state%diff_coefficient%data(k,iy,ix)

      if (nut > 0.0_dp .and. Ri_col(k) < 0.25_dp) then
         diss = ssq_col(k) * (nut - Ri_col(k)*kappa)
         diss = max(diss, 0.0_dp)
      else
         diss = 0.0_dp
      end if

      current_state%sgs_dissipation%data(k,iy,ix) = diss

      if (diss > 0.0_dp .and. nut > 0.0_dp) then
         L2 = nut / sqrt(max(1.0_dp - Ri_col(k)*(kappa/nut), 1.0e-12_dp))
         current_state%sgs_tke%data(k,iy,ix) = ((diss*diss*L2)**(1.0_dp/3.0_dp)) / 0.23_dp
      else
         current_state%sgs_tke%data(k,iy,ix) = 0.0_dp
      end if

   end do

   deallocate(ssq_col, Ri_col)

end subroutine timestep

!====================================================================
subroutine finalise(current_state)
   type(model_state_type), target, intent(inout) :: current_state
end subroutine finalise

end module simple_tke_mod
