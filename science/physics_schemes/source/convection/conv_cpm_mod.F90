! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module conv_cpm_mod

use um_types,              only: real_umphys
use water_constants_mod,   only: hcapv, hcapw, hcapi
use planet_constants_mod,  only: cpd => cp
use convection_config_mod, only: conv_cp_none,                                 &
                                 conv_cp_dry,                                  &
                                 conv_cp_moist

implicit none

public :: cpv_cpm, cl_cpm, ci_cpm
public :: set_conv_cp_coeffs

real(kind=real_umphys) :: cpv_cpm = 0.0_real_umphys
real(kind=real_umphys) :: cl_cpm  = 0.0_real_umphys
real(kind=real_umphys) :: ci_cpm  = 0.0_real_umphys

contains

subroutine set_conv_cp_coeffs(conv_mode)

integer, intent(in) :: conv_mode

select case (conv_mode)
  case (conv_cp_none)
    cpv_cpm = 0.0_real_umphys
    cl_cpm  = 0.0_real_umphys
    ci_cpm  = 0.0_real_umphys
  case (conv_cp_dry)
    cpv_cpm = cpd
    cl_cpm  = cpd
    ci_cpm  = cpd
  case (conv_cp_moist)
    cpv_cpm = hcapv
    cl_cpm  = hcapw
    ci_cpm  = hcapi
  case default
    cpv_cpm = 0.0_real_umphys
    cl_cpm  = 0.0_real_umphys
    ci_cpm  = 0.0_real_umphys
end select

end subroutine set_conv_cp_coeffs

end module conv_cpm_mod