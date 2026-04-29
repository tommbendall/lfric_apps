! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module bl_cpm_mod

use um_types,             only: real_umphys, r_bl
use iso_fortran_env,      only: real32
use water_constants_mod,  only: hcapv, hcapw, hcapi
use planet_constants_mod, only: cpd => cp, cp_bl
use blayer_config_mod,    only: bl_cp_none,                        &
                                bl_cp_dry,                         &
                                bl_cp_moist

implicit none

! Coefficients at real_umphys precision (for buoy_tq 64b)
public :: cpv_cpm, cl_cpm, ci_cpm
! Coefficients at r_bl precision (for bdy_impl3, kmkhz_9c)
public :: cpv_cpm_bl, cl_cpm_bl, ci_cpm_bl
! Coefficients at real32 precision (for buoy_tq 32b)
public :: cpv_cpm_32b, cl_cpm_32b, ci_cpm_32b
public :: set_bl_cp_coeffs

real(kind=real_umphys) :: cpv_cpm = 0.0_real_umphys
real(kind=real_umphys) :: cl_cpm  = 0.0_real_umphys
real(kind=real_umphys) :: ci_cpm  = 0.0_real_umphys

real(kind=r_bl) :: cpv_cpm_bl = 0.0_r_bl
real(kind=r_bl) :: cl_cpm_bl  = 0.0_r_bl
real(kind=r_bl) :: ci_cpm_bl  = 0.0_r_bl

real(real32) :: cpv_cpm_32b = 0.0
real(real32) :: cl_cpm_32b  = 0.0
real(real32) :: ci_cpm_32b  = 0.0

contains

subroutine set_bl_cp_coeffs(bl_mode)

integer, intent(in) :: bl_mode

select case (bl_mode)
  case (bl_cp_none)
    cpv_cpm     = 0.0_real_umphys
    cl_cpm      = 0.0_real_umphys
    ci_cpm      = 0.0_real_umphys
    cpv_cpm_bl  = 0.0_r_bl
    cl_cpm_bl   = 0.0_r_bl
    ci_cpm_bl   = 0.0_r_bl
    cpv_cpm_32b = 0.0
    cl_cpm_32b  = 0.0
    ci_cpm_32b  = 0.0
  case (bl_cp_dry)
    cpv_cpm     = cpd
    cl_cpm      = cpd
    ci_cpm      = cpd
    cpv_cpm_bl  = cp_bl
    cl_cpm_bl   = cp_bl
    ci_cpm_bl   = cp_bl
    cpv_cpm_32b = real(cpd, real32)
    cl_cpm_32b  = real(cpd, real32)
    ci_cpm_32b  = real(cpd, real32)
  case (bl_cp_moist)
    cpv_cpm     = hcapv
    cl_cpm      = hcapw
    ci_cpm      = hcapi
    cpv_cpm_bl  = real(hcapv, r_bl)
    cl_cpm_bl   = real(hcapw, r_bl)
    ci_cpm_bl   = real(hcapi, r_bl)
    cpv_cpm_32b = real(hcapv, real32)
    cl_cpm_32b  = real(hcapw, real32)
    ci_cpm_32b  = real(hcapi, real32)
  case default
    cpv_cpm     = 0.0_real_umphys
    cl_cpm      = 0.0_real_umphys
    ci_cpm      = 0.0_real_umphys
    cpv_cpm_bl  = 0.0_r_bl
    cl_cpm_bl   = 0.0_r_bl
    ci_cpm_bl   = 0.0_r_bl
    cpv_cpm_32b = 0.0
    cl_cpm_32b  = 0.0
    ci_cpm_32b  = 0.0
end select

end subroutine set_bl_cp_coeffs

end module bl_cpm_mod
