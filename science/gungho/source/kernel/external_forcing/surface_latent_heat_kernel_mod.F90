!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes the latent heat variables at the surface
!> @details Computes the latent heats of vaporisation and fusion at the surface
!!          of the model, using the potential temperature and Exner pressure.

module surface_latent_heat_kernel_mod

  use argument_mod,       only : arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL,                                      &
                                 GH_WRITE, GH_READ,                            &
                                 CELL_COLUMN,                                  &
                                 ANY_DISCONTINUOUS_SPACE_3
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
  type, public, extends(kernel_type) :: surface_latent_heat_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3),     & ! Lv
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3),     & ! Ls
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! exner_in_wth
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                        & ! theta
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! Lv0
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! Ls0
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cpv
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! cl
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 & ! ci
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  & ! Tm
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: surface_latent_heat_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: surface_latent_heat_code

contains

  !> @brief Computes latent heats at the model surface
  !> @param[in]     nlayers        The number of layers (in the 2D mesh)
  !> @param[in,out] Lv             Latent heat of vaporisation
  !> @param[in,out] Ls             Latent heat of sublimation
  !> @param[in]     exner_in_wth   Exner pressure at Wtheta points
  !> @param[in]     theta          Potential temperature at Wtheta points
  !> @param[in]     Lv0            Latent heat of vaporisation at Tm
  !> @param[in]     Ls0            Latent heat of sublimation at Tm
  !> @param[in]     cpv            Heat capacity of water vapour
  !> @param[in]     cl             Heat capacity of liquid water
  !> @param[in]     ci             Heat capacity of ice
  !> @param[in]     Tm             Reference temperature
  !> @param[in]     ndf_w3_2d      Number of DoFs for 2D W3 per cell
  !> @param[in]     undf_w3_2d     Number of DoFs in this partition for 2D W3
  !> @param[in]     map_w3_2d      Dofmap for 2D W3
  !> @param[in]     ndf_wtheta     Number of DoFs for Wtheta per cell
  !> @param[in]     undf_wtheta    Number of DoFs in this partition for Wtheta
  !> @param[in]     map_wtheta     Dofmap for Wtheta
  subroutine surface_latent_heat_code( nlayers,       &
                                       Lv,            &
                                       Ls,            &
                                       exner_in_wth,  &
                                       theta,         &
                                       Lv0,           &
                                       Ls0,           &
                                       cpv,           &
                                       cl,            &
                                       ci,            &
                                       Tm,            &
                                       ndf_w3_2d,     &
                                       undf_w3_2d,    &
                                       map_w3_2d,     &
                                       ndf_wtheta,    &
                                       undf_wtheta,   &
                                       map_wtheta )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_w3_2d
    integer(kind=i_def), intent(in)    :: undf_w3_2d
    integer(kind=i_def), intent(in)    :: ndf_wtheta
    integer(kind=i_def), intent(in)    :: undf_wtheta
    integer(kind=i_def), intent(in)    :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def), intent(in)    :: map_wtheta(ndf_wtheta)
    real(kind=r_def),    intent(inout) :: Lv(undf_w3_2d)
    real(kind=r_def),    intent(inout) :: Ls(undf_w3_2d)
    real(kind=r_def),    intent(in)    :: theta(undf_wtheta)
    real(kind=r_def),    intent(in)    :: exner_in_wth(undf_wtheta)
    real(kind=r_def),    intent(in)    :: Lv0
    real(kind=r_def),    intent(in)    :: Ls0
    real(kind=r_def),    intent(in)    :: cpv
    real(kind=r_def),    intent(in)    :: cl
    real(kind=r_def),    intent(in)    :: ci
    real(kind=r_def),    intent(in)    :: Tm

    ! Internal variables
    real(kind=r_def) :: T

    T = theta(map_wtheta(1))*exner_in_wth(map_wtheta(1))
    Lv(map_w3_2d(1)) = Lv0 - (cl - cpv)*(T - Tm)
    ! Add fusion latent heat to vaporisation to get sublimation
    Ls(map_w3_2d(1)) = Ls0 - (ci - cpv)*(T - Tm)

  end subroutine surface_latent_heat_code

end module surface_latent_heat_kernel_mod
