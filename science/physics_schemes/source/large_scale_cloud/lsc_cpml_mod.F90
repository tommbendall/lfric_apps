! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsc_cpml_mod

use um_types,             only: real_umphys
use water_constants_mod,  only: hcapv, hcapw, hcapi
use planet_constants_mod, only: cpd => cp
use cloud_config_mod,     only: lsc_moist_heat_cap_none,                       &
                                lsc_moist_heat_cap_dry,                        &
                                lsc_moist_heat_cap_moist

implicit none

public :: cpv_cpml, cl_cpml, ci_cpml
public :: set_lsc_moist_heat_cap_coeffs

real(kind=real_umphys) :: cpv_cpml = 0.0_real_umphys
real(kind=real_umphys) :: cl_cpml  = 0.0_real_umphys
real(kind=real_umphys) :: ci_cpml  = 0.0_real_umphys

contains

subroutine set_lsc_moist_heat_cap_coeffs(lsc_mode)

integer, intent(in) :: lsc_mode

select case (lsc_mode)
  case (lsc_moist_heat_cap_none)
    cpv_cpml = 0.0_real_umphys
    cl_cpml  = 0.0_real_umphys
    ci_cpml  = 0.0_real_umphys
  case (lsc_moist_heat_cap_dry)
    cpv_cpml = cpd
    cl_cpml  = cpd
    ci_cpml  = cpd
  case (lsc_moist_heat_cap_moist)
    cpv_cpml = hcapv
    cl_cpml  = hcapw
    ci_cpml  = hcapi
  case default
    cpv_cpml = 0.0_real_umphys
    cl_cpml  = 0.0_real_umphys
    ci_cpml  = 0.0_real_umphys
end select

end subroutine set_lsc_moist_heat_cap_coeffs

end module lsc_cpml_mod