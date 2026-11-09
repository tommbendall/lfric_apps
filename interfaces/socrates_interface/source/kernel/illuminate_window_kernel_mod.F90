!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------
!> @brief Mean solar angle and lit fraction over a window of the timestep

module illuminate_window_kernel_mod

use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_SCALAR,                              &
                             GH_REAL, GH_INTEGER,                              &
                             GH_READ, GH_WRITE,                                &
                             CELL_COLUMN,                                      &
                             ANY_DISCONTINUOUS_SPACE_1
use constants_mod,     only: r_def, i_def
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: illuminate_window_kernel_type
  private
  type(arg_type) :: meta_args(12) = (/                                         &
      arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! cos_zenith_angle_win
      arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! lit_fraction_win
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! sin_stellar_declination_rts
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! stellar_eqn_of_time_rts
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! latitude
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! longitude
      arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                & ! window_length
      arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                & ! window_start
      arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                & ! spin_option
      arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                & ! hour_angle_inc
      arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                & ! fixed_zenith_angle
      arg_type(GH_SCALAR, GH_REAL,    GH_READ)                                 & ! fixed_azimuth_angle
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: illuminate_window_code
end type

public :: illuminate_window_code

contains

!> @param[in]     nlayers                     Number of layers
!> @param[in,out] cos_zenith_angle_win        Cosine of the zenith angle
!> @param[in,out] lit_fraction_win            Lit fraction of the window
!> @param[in]     sin_stellar_declination_rts Stellar declination
!> @param[in]     stellar_eqn_of_time_rts     Stellar equation of time
!> @param[in]     latitude                    Latitude field
!> @param[in]     longitude                   Longitude field
!> @param[in]     window_length               Length of the window (s)
!> @param[in]     window_start                Seconds into the day at its start
!> @param[in]     spin_option                 Motion of the sun in the sky
!> @param[in]     hour_angle_inc              Rate of change of the hour angle
!> @param[in]     fixed_zenith_angle          Zenith angle for a fixed sun
!> @param[in]     fixed_azimuth_angle         Azimuth angle for a fixed sun
!> @param[in]     ndf_2d                      Num DoFs per cell for 2d space
!> @param[in]     undf_2d                     Num DoFs in this partition for 2d
!> @param[in]     map_2d                      Map of base cell DoFs for 2d
subroutine illuminate_window_code(                                             &
    nlayers,                                                                   &
    cos_zenith_angle_win, lit_fraction_win,                                    &
    sin_stellar_declination_rts, stellar_eqn_of_time_rts,                      &
    latitude, longitude,                                                       &
    window_length, window_start, spin_option, hour_angle_inc,                  &
    fixed_zenith_angle, fixed_azimuth_angle,                                   &
    ndf_2d, undf_2d, map_2d                                                    &
)

  use orbit_config_mod,    only: spin_user, spin_earth_day, spin_fixed_sun
  use socrates_illuminate, only: illuminate,                                   &
                                 ip_spin_user,                                 &
                                 ip_spin_earth_day,                            &
                                 ip_spin_fixed_sun

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_2d
  integer(kind=i_def), intent(in)    :: undf_2d
  integer(kind=i_def), intent(in)    :: spin_option

  real(kind=r_def),    intent(in)    :: window_length
  real(kind=r_def),    intent(in)    :: window_start
  real(kind=r_def),    intent(in)    :: hour_angle_inc
  real(kind=r_def),    intent(in)    :: fixed_zenith_angle
  real(kind=r_def),    intent(in)    :: fixed_azimuth_angle

  real(kind=r_def),    intent(inout) :: cos_zenith_angle_win(undf_2d)
  real(kind=r_def),    intent(inout) :: lit_fraction_win(undf_2d)
  real(kind=r_def),    intent(in)    :: sin_stellar_declination_rts(undf_2d)
  real(kind=r_def),    intent(in)    :: stellar_eqn_of_time_rts(undf_2d)
  real(kind=r_def),    intent(in)    :: latitude(undf_2d)
  real(kind=r_def),    intent(in)    :: longitude(undf_2d)

  integer(kind=i_def), intent(in)    :: map_2d(ndf_2d)

  ! Internal variables
  integer(kind=i_def), parameter :: n_profile = 1_i_def
  integer(kind=i_def)            :: i_spin
  real(kind=r_def)               :: sin_stellar_declination
  real(kind=r_def)               :: stellar_eqn_of_time

  ! Local copies, as illuminate has these as intent(inout)
  sin_stellar_declination = sin_stellar_declination_rts(map_2d(1))
  stellar_eqn_of_time = stellar_eqn_of_time_rts(map_2d(1))

  ! Set motion of sun across the sky
  select case (spin_option)
  case (spin_user)
    i_spin = ip_spin_user
  case (spin_earth_day)
    i_spin = ip_spin_earth_day
  case (spin_fixed_sun)
    i_spin = ip_spin_fixed_sun
  case default
    i_spin = ip_spin_earth_day
  end select

  call illuminate(                                                             &
      l_stellar_angle         = .true.,                                        &
      n_profile               = n_profile,                                     &
      i_spin                  = i_spin,                                        &
      second_of_day           = window_start,                                  &
      length_of_timestep      = window_length,                                 &
      hour_angle_inc          = hour_angle_inc,                                &
      fixed_zenith_angle      = fixed_zenith_angle,                            &
      fixed_azimuth_angle     = fixed_azimuth_angle,                           &
      latitude                = latitude(map_2d(1):map_2d(1)),                 &
      longitude               = longitude(map_2d(1):map_2d(1)),                &
      sin_stellar_declination = sin_stellar_declination,                       &
      stellar_eqn_of_time     = stellar_eqn_of_time,                           &
      cos_zenith_angle        = cos_zenith_angle_win(map_2d(1):map_2d(1)),     &
      lit_fraction            = lit_fraction_win(map_2d(1):map_2d(1))          &
  )

end subroutine illuminate_window_code

end module illuminate_window_kernel_mod
