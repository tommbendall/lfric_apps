! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsc_cpm_mod

use um_types,             only: real_umphys
use water_constants_mod,  only: hcapv, hcapw, hcapi
use planet_constants_mod, only: cpd => cp
use cloud_config_mod,     only: lsc_cp_none,                       &
                                lsc_cp_dry,                        &
                                lsc_cp_moist

implicit none

public :: cpv_cpm, cl_cpm, ci_cpm
public :: set_lsc_cp_coeffs

real(kind=real_umphys) :: cpv_cpm = 0.0_real_umphys
real(kind=real_umphys) :: cl_cpm  = 0.0_real_umphys
real(kind=real_umphys) :: ci_cpm  = 0.0_real_umphys

contains

subroutine set_lsc_cp_coeffs(lsc_mode)

integer, intent(in) :: lsc_mode

select case (lsc_mode)
  case (lsc_cp_none)
    cpv_cpm = 0.0_real_umphys
    cl_cpm  = 0.0_real_umphys
    ci_cpm  = 0.0_real_umphys
  case (lsc_cp_dry)
    cpv_cpm = cpd
    cl_cpm  = cpd
    ci_cpm  = cpd
  case (lsc_cp_moist)
    cpv_cpm = hcapv
    cl_cpm  = hcapw
    ci_cpm  = hcapi
  case default
    cpv_cpm = 0.0_real_umphys
    cl_cpm  = 0.0_real_umphys
    ci_cpm  = 0.0_real_umphys
end select

end subroutine set_lsc_cp_coeffs

end module lsc_cpm_mod