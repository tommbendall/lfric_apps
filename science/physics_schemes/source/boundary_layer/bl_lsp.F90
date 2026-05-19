! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Convert temperature from liquid ice to liquid, and convert
!           the vapour+liquid+ice variable (Q) to vapour+liquid. This
!           subroutine is used if the mixed phase precipitation scheme
!           is selected and a full boundary layer treatment is not
!           performed.

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

module bl_lsp_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'BL_LSP_MOD'
contains

subroutine bl_lsp( bl_levels,qcf,q,t,qcl )

use atm_fields_bounds_mod, only: tdims
use planet_constants_mod, only: lsrcp, cpd => cp
use water_constants_mod, only: lc, lf, tm
use bl_cpm_mod, only: cpv_cpm, cl_cpm, ci_cpm
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

integer, intent(in) ::                                                         &
  bl_levels             ! in   Number of boundary layer levels

real(kind=real_umphys), intent(in out) ::                                      &
  qcf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                     &
      bl_levels),                                                              &
                                 ! INOUT Ice water content
  q(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                       &
    bl_levels),                                                                &
                                 ! INOUT
!                                  in    Vapour+liquid+ice content
!                                  out   Vapour+liquid content
    t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                     &
      bl_levels)                   ! INOUT
!                                  in    Liquid ice temperature
!                                  out   Liquid temperature
real(kind=real_umphys), intent(in) ::                                          &
  qcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                     &
      bl_levels)
                                 ! IN Cloud liquid water content
! Temporary Space
integer ::                                                                     &
        i,                                                                     &
                               ! Counter over points
        j,                                                                     &
                               ! Counter over points
        k                ! Counter over boundary layer levels
real(kind=real_umphys) :: newqcf              ! Temporary variable for QCF
real(kind=real_umphys) :: L_sub_val           ! Temp-dependent latent heat of sublimation
real(kind=real_umphys) :: cp_moist_val        ! Temp-dependent moist specific heat
real(kind=real_umphys) :: lsrcp_moist         ! Temp-dependent L_sub / cp_moist

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BL_LSP'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP          private(i,j,k,newqcf,L_sub_val,cp_moist_val,lsrcp_moist)        &
!$OMP          SHARED(bl_levels,tdims,q,qcf,t,qcl,lsrcp,cpd,cpv_cpm,cl_cpm,    &
!$OMP                 ci_cpm)
do k = 1, bl_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Convert Q (vapour+liquid+ice) to (vapour+liquid)
      q(i,j,k)=q(i,j,k)-qcf(i,j,k)
      ! Check that Q is not negative
      if (q(i,j,k)  <   0.0) then
        ! Evaporate ice to keep Q positive, but don't let ice go negative
        ! itself
        newqcf=max(qcf(i,j,k)+q(i,j,k),0.0)
        q(i,j,k)=q(i,j,k)+(qcf(i,j,k)-newqcf)
        qcf(i,j,k)=newqcf
      end if
      ! Adjust T from T liquid ice to T liquid
      L_sub_val    = (lc + lf) - (ci_cpm - cpv_cpm) * (t(i,j,k) - tm)
      cp_moist_val = cpd + (q(i,j,k)-qcl(i,j,k)) * cpv_cpm                    &
           + qcl(i,j,k) * cl_cpm + qcf(i,j,k) * ci_cpm
      lsrcp_moist  = L_sub_val / cp_moist_val
      t(i,j,k)=t(i,j,k)+lsrcp_moist*qcf(i,j,k)
    end do
  end do
end do
!$OMP end PARALLEL do
! End the subroutine
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bl_lsp
end module bl_lsp_mod
