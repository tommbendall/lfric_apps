! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_cpml_mod

use um_types,             only: real_lsprec
use water_constants_mod,  only: hcapv, hcapw, hcapi
use planet_constants_mod,    only: cpd => cp
use microphysics_config_mod, only: lsp_moist_heat_cap_none,                    &
                                   lsp_moist_heat_cap_dry,                     &
                                   lsp_moist_heat_cap_moist

implicit none

public :: cpv_cpml, cl_cpml, ci_cpml
public :: set_lsp_moist_heat_cap_coeffs

real(kind=real_lsprec) :: cpv_cpml = 0.0_real_lsprec
real(kind=real_lsprec) :: cl_cpml  = 0.0_real_lsprec
real(kind=real_lsprec) :: ci_cpml  = 0.0_real_lsprec

contains

subroutine set_lsp_moist_heat_cap_coeffs(lsc_mode)

integer, intent(in) :: lsc_mode

select case (lsc_mode)
  case (lsp_moist_heat_cap_none)
    cpv_cpml = 0.0_real_lsprec
    cl_cpml  = 0.0_real_lsprec
    ci_cpml  = 0.0_real_lsprec
  case (lsp_moist_heat_cap_dry)
    cpv_cpml = real(cpd, kind=real_lsprec)
    cl_cpml  = real(cpd, kind=real_lsprec)
    ci_cpml  = real(cpd, kind=real_lsprec)
  case (lsp_moist_heat_cap_moist)
    cpv_cpml = real(hcapv, kind=real_lsprec)
    cl_cpml  = real(hcapw, kind=real_lsprec)
    ci_cpml  = real(hcapi, kind=real_lsprec)
  case default
    cpv_cpml = 0.0_real_lsprec
    cl_cpml  = 0.0_real_lsprec
    ci_cpml  = 0.0_real_lsprec
end select

end subroutine set_lsp_moist_heat_cap_coeffs

end module lsp_cpml_mod
