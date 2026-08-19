!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes the heat capacity of moist air at constant volume

module compute_cv_kernel_mod

  use argument_mod,       only : arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL,                                      &
                                 GH_WRITE, GH_READ,                            &
                                 DOF
  use constants_mod,      only : r_def, i_def
  use fs_continuity_mod,  only : Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_cv_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),                        & ! cv_tot
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! mr_vap
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! mr_liq
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! mr_ice
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cvd
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cvv
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cl
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  & ! ci
    /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: compute_cv_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_cv_code

contains

  !> @brief Computes latent heats at the model surface
  !> @param[in]     nlayers        The number of layers
  !> @param[in,out] cv_tot         Heat capacity
  !> @param[in]     mr_vap         Mixing ratio of water vapour
  !> @param[in]     mr_liq         Mixing ratio of liquid water
  !> @param[in]     mr_ice         Mixing ratio of ice
  !> @param[in]     cvd            Heat capacity of dry air at const volume
  !> @param[in]     cvv            Heat capacity of water vapour at const volume
  !> @param[in]     cl             Heat capacity of liquid water
  !> @param[in]     ci             Heat capacity of ice
  subroutine compute_cv_code( cvm_tot,       &
                              mr_vap,        &
                              mr_liq,        &
                              mr_ice,        &
                              cvd,           &
                              cvv,           &
                              cl,            &
                              ci             )

    implicit none

    ! Arguments
    real(kind=r_def),    intent(inout) :: cvm_tot
    real(kind=r_def),    intent(in)    :: mr_vap
    real(kind=r_def),    intent(in)    :: mr_liq
    real(kind=r_def),    intent(in)    :: mr_ice
    real(kind=r_def),    intent(in)    :: cvd
    real(kind=r_def),    intent(in)    :: cvv
    real(kind=r_def),    intent(in)    :: cl
    real(kind=r_def),    intent(in)    :: ci

    ! Compute the heat capacity of moist air at constant volume
    cvm_tot = cvd + mr_vap * cvv + mr_liq * cl + mr_ice * ci

  end subroutine compute_cv_code

end module compute_cv_kernel_mod
