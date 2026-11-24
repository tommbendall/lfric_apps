!-------------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief   Shifts longitudes to the east of the dateline by 2*pi.
!> @details Longitudes below a threshold value have 2*pi added to them, so that
!!          a limited area domain spanning the dateline has a continuous
!!          longitude coordinate. Other values are copied unchanged.
module shift_longitude_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, GH_REAL,         &
                               GH_READ, GH_WRITE, DOF, ANY_SPACE_1
  use constants_mod,     only: r_def, PI
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: shift_longitude_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_SPACE_1),                  &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1),                  &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                 &
    /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: shift_longitude_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: shift_longitude_code

contains

!> @brief Adds 2*pi to longitudes that are below a threshold
!> @param[in,out] chi_shifted  Adjusted longitude
!> @param[in]     chi          Original longitude
!> @param[in]     threshold    Longitudes below this value are shifted
subroutine shift_longitude_code( chi_shifted, chi, threshold )

  implicit none

  ! Arguments
  real(kind=r_def), intent(inout) :: chi_shifted
  real(kind=r_def), intent(in)    :: chi
  real(kind=r_def), intent(in)    :: threshold

  if ( chi < threshold ) then
    chi_shifted = chi + 2.0_r_def * PI
  else
    chi_shifted = chi
  end if

end subroutine shift_longitude_code

end module shift_longitude_kernel_mod
