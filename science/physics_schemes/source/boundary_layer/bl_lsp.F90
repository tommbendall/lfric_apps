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

subroutine bl_lsp( bl_levels,qcf,q,t,qcl,qrain,qgraupel )

use atm_fields_bounds_mod, only: tdims
use planet_constants_mod, only: cpd => cp
use bl_cpm_mod, only: cpv_cpm, cl_cpm, ci_cpm
use water_constants_mod, only: tm => tm_bl, lc => lc_bl, lf => lf_bl
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
      bl_levels),                                                              &
                                 ! IN Cloud liquid water content
  qrain(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                                 ! IN Rain mixing ratio
  qgraupel(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           bl_levels)
                                 ! IN Graupel mixing ratio
! Temporary Space
integer ::                                                                     &
        i,                                                                     &
                               ! Counter over points
        j,                                                                     &
                               ! Counter over points
        k                ! Counter over boundary layer levels
real(kind=real_umphys) :: newqcf              ! Temporary variable for QCF
real(kind=real_umphys) :: cpm                 ! Moist heat capacity before phase change (qcf as ice)
real(kind=real_umphys) :: cpm_dag             ! Moist heat capacity after phase change (qcf as vapour)
real(kind=real_umphys) :: lrs0                ! Reference latent heat of sublimation

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BL_LSP'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
lrs0 = (lc + lf) + (ci_cpm - cpv_cpm) * tm
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP          private(i,j,k,newqcf,cpm,cpm_dag)                               &
!$OMP          SHARED(bl_levels,tdims,q,qcf,t,qcl,qrain,qgraupel,              &
!$OMP                 cpd,cpv_cpm,cl_cpm,ci_cpm,lrs0)
do k = 1, bl_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Compute cpm BEFORE phase change for ice-liquid temperature
      cpm = cpd + cpv_cpm*(q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))                 &
                + cl_cpm*(qcl(i,j,k) + qrain(i,j,k))                           &
                + ci_cpm*(qcf(i,j,k) + qgraupel(i,j,k))
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
      ! Compute cpm_dag AFTER phase change for liquid temperature
      cpm_dag = cpd + cpv_cpm*(q(i,j,k) - qcl(i,j,k))                          &
                    + cl_cpm*(qcl(i,j,k) + qrain(i,j,k))                       &
                    + ci_cpm*(qcf(i,j,k) + qgraupel(i,j,k))
      ! Adjust T from T liquid ice to T liquid
      t(i,j,k) = (cpm/cpm_dag) * t(i,j,k) + (lrs0/cpm_dag)*qcf(i,j,k)
    end do
  end do
end do
!$OMP end PARALLEL do
! End the subroutine
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bl_lsp
end module bl_lsp_mod
