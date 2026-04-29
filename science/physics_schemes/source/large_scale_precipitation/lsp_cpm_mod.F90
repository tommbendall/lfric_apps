! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_cpm_mod

use um_types,             only: real_lsprec
use water_constants_mod,  only: hcapv, hcapw, hcapi
use planet_constants_mod,    only: cpd => cp
use microphysics_config_mod, only: lsp_cp_none,                    &
                                   lsp_cp_dry,                     &
                                   lsp_cp_moist

implicit none

public :: cpv_cpm, cl_cpm, ci_cpm
public :: set_lsp_cp_coeffs

real(kind=real_lsprec) :: cpv_cpm = 0.0_real_lsprec
real(kind=real_lsprec) :: cl_cpm  = 0.0_real_lsprec
real(kind=real_lsprec) :: ci_cpm  = 0.0_real_lsprec

contains

subroutine set_lsp_cp_coeffs(lsc_mode)

integer, intent(in) :: lsc_mode

select case (lsc_mode)
  case (lsp_cp_none)
    cpv_cpm = 0.0_real_lsprec
    cl_cpm  = 0.0_real_lsprec
    ci_cpm  = 0.0_real_lsprec
  case (lsp_cp_dry)
    cpv_cpm = real(cpd, kind=real_lsprec)
    cl_cpm  = real(cpd, kind=real_lsprec)
    ci_cpm  = real(cpd, kind=real_lsprec)
  case (lsp_cp_moist)
    cpv_cpm = real(hcapv, kind=real_lsprec)
    cl_cpm  = real(hcapw, kind=real_lsprec)
    ci_cpm  = real(hcapi, kind=real_lsprec)
  case default
    cpv_cpm = 0.0_real_lsprec
    cl_cpm  = 0.0_real_lsprec
    ci_cpm  = 0.0_real_lsprec
end select

end subroutine set_lsp_cp_coeffs

end module lsp_cpm_mod
