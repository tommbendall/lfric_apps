! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module bl_cpml_mod

use um_types,             only: real_umphys, r_bl
use iso_fortran_env,      only: real32
use water_constants_mod,  only: hcapv, hcapw, hcapi
use planet_constants_mod, only: cpd => cp, cp_bl
use blayer_config_mod,    only: bl_moist_heat_cap_none,                        &
                                bl_moist_heat_cap_dry,                         &
                                bl_moist_heat_cap_moist

implicit none

! Coefficients at real_umphys precision (for buoy_tq 64b)
public :: cpv_cpml, cl_cpml, ci_cpml
! Coefficients at r_bl precision (for bdy_impl3, kmkhz_9c)
public :: cpv_cpml_bl, cl_cpml_bl, ci_cpml_bl
! Coefficients at real32 precision (for buoy_tq 32b)
public :: cpv_cpml_32b, cl_cpml_32b, ci_cpml_32b
public :: set_bl_moist_heat_cap_coeffs

real(kind=real_umphys) :: cpv_cpml = 0.0_real_umphys
real(kind=real_umphys) :: cl_cpml  = 0.0_real_umphys
real(kind=real_umphys) :: ci_cpml  = 0.0_real_umphys

real(kind=r_bl) :: cpv_cpml_bl = 0.0_r_bl
real(kind=r_bl) :: cl_cpml_bl  = 0.0_r_bl
real(kind=r_bl) :: ci_cpml_bl  = 0.0_r_bl

real(real32) :: cpv_cpml_32b = 0.0
real(real32) :: cl_cpml_32b  = 0.0
real(real32) :: ci_cpml_32b  = 0.0

contains

subroutine set_bl_moist_heat_cap_coeffs(bl_mode)

integer, intent(in) :: bl_mode

select case (bl_mode)
  case (bl_moist_heat_cap_none)
    cpv_cpml     = 0.0_real_umphys
    cl_cpml      = 0.0_real_umphys
    ci_cpml      = 0.0_real_umphys
    cpv_cpml_bl  = 0.0_r_bl
    cl_cpml_bl   = 0.0_r_bl
    ci_cpml_bl   = 0.0_r_bl
    cpv_cpml_32b = 0.0
    cl_cpml_32b  = 0.0
    ci_cpml_32b  = 0.0
  case (bl_moist_heat_cap_dry)
    cpv_cpml     = cpd
    cl_cpml      = cpd
    ci_cpml      = cpd
    cpv_cpml_bl  = cp_bl
    cl_cpml_bl   = cp_bl
    ci_cpml_bl   = cp_bl
    cpv_cpml_32b = real(cpd, real32)
    cl_cpml_32b  = real(cpd, real32)
    ci_cpml_32b  = real(cpd, real32)
  case (bl_moist_heat_cap_moist)
    cpv_cpml     = hcapv
    cl_cpml      = hcapw
    ci_cpml      = hcapi
    cpv_cpml_bl  = real(hcapv, r_bl)
    cl_cpml_bl   = real(hcapw, r_bl)
    ci_cpml_bl   = real(hcapi, r_bl)
    cpv_cpml_32b = real(hcapv, real32)
    cl_cpml_32b  = real(hcapw, real32)
    ci_cpml_32b  = real(hcapi, real32)
  case default
    cpv_cpml     = 0.0_real_umphys
    cl_cpml      = 0.0_real_umphys
    ci_cpml      = 0.0_real_umphys
    cpv_cpml_bl  = 0.0_r_bl
    cl_cpml_bl   = 0.0_r_bl
    ci_cpml_bl   = 0.0_r_bl
    cpv_cpml_32b = 0.0
    cl_cpml_32b  = 0.0
    ci_cpml_32b  = 0.0
end select

end subroutine set_bl_moist_heat_cap_coeffs

end module bl_cpml_mod
