!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------

!> @brief   Finite difference predictor for the transport wind.
!> @details Increments a W2 field of fluxes by alpha_dt times the acceleration
!!          a = -(cp*theta_v*grad(exner) + grad(geopotential)) - 2 Omega x u
!!          for lowest-order elements with native coordinates.
!!          Each cell adds its own half of the increment to each of its faces,
!!          using only its own reconstructed gradients and wind, so that the
!!          two cells sharing a face together give the full face value. The
!!          kernel must therefore be called on owned cells and the depth 1
!!          halo cells, and sums into the W2 field (GH_INC).
!!          Faces on the model top and bottom boundaries are not changed, and
!!          nor are side faces with no neighbouring cell.
!!          The horizontal geometry comes from the 2D constants computed by
!!          calc_fd_predictor_geometry_kernel_mod, with lengths at unit radius.
!!          The DoFs of each face are given by sigma*(outward flux), where
!!          sigma is -1 for DoFs 1, 4, 5 and +1 for DoFs 2, 3, 6.

module fd_wind_predictor_kernel_mod

use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_SCALAR, GH_REAL, GH_INTEGER,         &
                             GH_LOGICAL, GH_READ, GH_INC,                      &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             STENCIL, CROSS2D, CELL_COLUMN
use constants_mod,     only: r_def, i_def, l_def
use fs_continuity_mod, only: W2, W3, Wtheta
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!! Psy layer.
type, public, extends(kernel_type) :: fd_wind_predictor_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/                                         &
      arg_type(GH_FIELD,  GH_REAL,    GH_INC,  W2),                            &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, W3, STENCIL(CROSS2D)),          &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, W3, STENCIL(CROSS2D)),          &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, Wtheta, STENCIL(CROSS2D)),      &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, W3, STENCIL(CROSS2D)),          &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, Wtheta, STENCIL(CROSS2D)),      &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, W2),                            &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, W2),                            &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_DISCONTINUOUS_SPACE_1),     &
      arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                &
      arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                &
      arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                &
      arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &
      arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                                &
      arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                                 &
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: fd_wind_predictor_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: fd_wind_predictor_code

contains

!> @brief Adds this cell's contribution to the predicted wind for one column
!! @param[in]     nlayers          Number of layers
!! @param[in,out] u_pred           Predicted wind (fluxes)
!! @param[in]     exner            Exner pressure
!! @param[in]     exner_smap_size  Size of the stencil map for each branch
!! @param[in]     exner_smap_max   Max size of the stencil map
!! @param[in]     exner_smap       W3 stencil map
!! @param[in]     phi              Geopotential
!! @param[in]     phi_smap_size    Size of the stencil map for each branch
!! @param[in]     phi_smap_max     Max size of the stencil map
!! @param[in]     phi_smap         W3 stencil map
!! @param[in]     theta_v          Virtual potential temperature
!! @param[in]     theta_smap_size  Size of the stencil map for each branch
!! @param[in]     theta_smap_max   Max size of the stencil map
!! @param[in]     theta_smap       Wtheta stencil map
!! @param[in]     height_w3        Height of the W3 DoFs
!! @param[in]     hw3_smap_size    Size of the stencil map for each branch
!! @param[in]     hw3_smap_max     Max size of the stencil map
!! @param[in]     hw3_smap         W3 stencil map
!! @param[in]     height_wth       Height of the Wtheta DoFs
!! @param[in]     hwt_smap_size    Size of the stencil map for each branch
!! @param[in]     hwt_smap_max     Max size of the stencil map
!! @param[in]     hwt_smap         Wtheta stencil map
!! @param[in]     u                Wind at the current time (fluxes)
!! @param[in]     dA_at_w2         Face areas at W2 DoFs
!! @param[in]     fd_geom          2D constants, with fd_geom_ndata per column
!! @param[in]     alpha_dt         Weight of the acceleration times time step
!! @param[in]     cp               Specific heat of dry air at constant pressure
!! @param[in]     scaled_radius    Planet radius
!! @param[in]     geometry         Integer indicating the domain geometry
!! @param[in]     shallow          Whether to use the shallow-atmosphere radius
!! @param[in]     rotating         Whether to include the Coriolis term
!! @param[in]     ndf_w2           Num DoFs per cell for W2
!! @param[in]     undf_w2          Num DoFs in this partition for W2
!! @param[in]     map_w2           Map of base cell DoFs for W2
!! @param[in]     ndf_w3           Num DoFs per cell for W3
!! @param[in]     undf_w3          Num DoFs in this partition for W3
!! @param[in]     map_w3           Map of base cell DoFs for W3
!! @param[in]     ndf_wt           Num DoFs per cell for Wtheta
!! @param[in]     undf_wt          Num DoFs in this partition for Wtheta
!! @param[in]     map_wt           Map of base cell DoFs for Wtheta
!! @param[in]     ndf_2d           Num DoFs per cell for the 2D constants
!! @param[in]     undf_2d          Num DoFs in this partition for 2D constants
!! @param[in]     map_2d           Map of base cell DoFs for 2D constants
subroutine fd_wind_predictor_code(                                             &
    nlayers,                                                                   &
    u_pred,                                                                    &
    exner, exner_smap_size, exner_smap_max, exner_smap,                        &
    phi, phi_smap_size, phi_smap_max, phi_smap,                                &
    theta_v, theta_smap_size, theta_smap_max, theta_smap,                      &
    height_w3, hw3_smap_size, hw3_smap_max, hw3_smap,                          &
    height_wth, hwt_smap_size, hwt_smap_max, hwt_smap,                         &
    u, dA_at_w2, fd_geom,                                                      &
    alpha_dt, cp, scaled_radius, geometry, shallow, rotating,                  &
    ndf_w2, undf_w2, map_w2,                                                   &
    ndf_w3, undf_w3, map_w3,                                                   &
    ndf_wt, undf_wt, map_wt,                                                   &
    ndf_2d, undf_2d, map_2d                                                    &
)

  use base_mesh_config_mod, only: geometry_spherical

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: map_wt(ndf_wt)
  integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

  integer(kind=i_def), intent(in) :: exner_smap_max, phi_smap_max
  integer(kind=i_def), intent(in) :: theta_smap_max
  integer(kind=i_def), intent(in) :: hw3_smap_max, hwt_smap_max
  integer(kind=i_def), intent(in) :: exner_smap_size(4), phi_smap_size(4)
  integer(kind=i_def), intent(in) :: theta_smap_size(4)
  integer(kind=i_def), intent(in) :: hw3_smap_size(4), hwt_smap_size(4)
  integer(kind=i_def), intent(in) :: exner_smap(ndf_w3, exner_smap_max, 4)
  integer(kind=i_def), intent(in) :: phi_smap(ndf_w3, phi_smap_max, 4)
  integer(kind=i_def), intent(in) :: theta_smap(ndf_wt, theta_smap_max, 4)
  integer(kind=i_def), intent(in) :: hw3_smap(ndf_w3, hw3_smap_max, 4)
  integer(kind=i_def), intent(in) :: hwt_smap(ndf_wt, hwt_smap_max, 4)

  real(kind=r_def), intent(inout) :: u_pred(undf_w2)
  real(kind=r_def), intent(in)    :: exner(undf_w3)
  real(kind=r_def), intent(in)    :: phi(undf_w3)
  real(kind=r_def), intent(in)    :: theta_v(undf_wt)
  real(kind=r_def), intent(in)    :: height_w3(undf_w3)
  real(kind=r_def), intent(in)    :: height_wth(undf_wt)
  real(kind=r_def), intent(in)    :: u(undf_w2)
  real(kind=r_def), intent(in)    :: dA_at_w2(undf_w2)
  real(kind=r_def), intent(in)    :: fd_geom(undf_2d)
  real(kind=r_def), intent(in)    :: alpha_dt
  real(kind=r_def), intent(in)    :: cp
  real(kind=r_def), intent(in)    :: scaled_radius
  integer(kind=i_def), intent(in) :: geometry
  logical(kind=l_def), intent(in) :: shallow
  logical(kind=l_def), intent(in) :: rotating

  ! Sign of each W2 DoF relative to the outward normal
  real(kind=r_def), parameter :: sigma(6) = (/ -1.0_r_def, 1.0_r_def,          &
      1.0_r_def, -1.0_r_def, -1.0_r_def, 1.0_r_def /)

  ! Internal variables
  integer(kind=i_def) :: df, k, base, nb_w3
  logical(kind=l_def) :: has_neighbour(4)

  real(kind=r_def) :: normal(2, 4), coeff(4), correction(2, 4)
  real(kind=r_def) :: w_grad(2, 4), w_wind(2, 4), omega(3)
  real(kind=r_def) :: rad(0:nlayers - 1), rad_lev(0:nlayers)
  real(kind=r_def) :: pz_face(nlayers), pz(0:nlayers - 1)
  real(kind=r_def) :: px(0:nlayers - 1), py(0:nlayers - 1)
  real(kind=r_def) :: zx(0:nlayers), zy(0:nlayers)
  real(kind=r_def) :: uc(0:nlayers - 1), vc(0:nlayers - 1), wc(0:nlayers - 1)
  real(kind=r_def) :: ax(0:nlayers - 1), ay(0:nlayers - 1), az(0:nlayers - 1)
  real(kind=r_def) :: x_bot(0:nlayers - 1), x_top(0:nlayers - 1)
  real(kind=r_def) :: theta_f, d_val, dz, un, w_bot, w_top, flux
  real(kind=r_def) :: dz_nb, dh, sum_x, sum_y

  ! Constants for this column
  do df = 1, 4
    base = map_2d(1) + (df - 1) * 9 - 1
    normal(1, df) = fd_geom(base + 1)
    normal(2, df) = fd_geom(base + 2)
    coeff(df) = fd_geom(base + 3)
    correction(1, df) = fd_geom(base + 4)
    correction(2, df) = fd_geom(base + 5)
    w_grad(1, df) = fd_geom(base + 6)
    w_grad(2, df) = fd_geom(base + 7)
    w_wind(1, df) = fd_geom(base + 8)
    w_wind(2, df) = fd_geom(base + 9)
    has_neighbour(df) = ( exner_smap_size(df) > 1 )
  end do
  base = map_2d(1) + 36 - 1
  omega(:) = fd_geom(base + 1 : base + 3)

  ! Radius for horizontal lengths. Heights are above the surface
  if ( geometry == geometry_spherical ) then
    if ( shallow ) then
      rad(:) = scaled_radius
      rad_lev(:) = scaled_radius
    else
      do k = 0, nlayers - 1
        rad(k) = scaled_radius + height_w3(map_w3(1) + k)
      end do
      do k = 0, nlayers
        rad_lev(k) = scaled_radius + height_wth(map_wt(1) + k)
      end do
    end if
  else
    rad(:) = 1.0_r_def
    rad_lev(:) = 1.0_r_def
  end if

  ! Pressure-gradient forcing terms X = cp*theta*exner + phi at the bottom and
  ! top face of each cell
  do k = 0, nlayers - 1
    x_bot(k) = cp * theta_v(map_wt(1) + k) * exner(map_w3(1) + k)              &
             + phi(map_w3(1) + k)
    x_top(k) = cp * theta_v(map_wt(2) + k) * exner(map_w3(1) + k)              &
             + phi(map_w3(1) + k)
  end do

  ! Vertical pressure gradient at interior faces, and the mean over each cell
  ! (one-sided at the boundaries)
  do k = 1, nlayers - 1
    theta_f = theta_v(map_wt(1) + k)
    pz_face(k) = ( cp * theta_f * ( exner(map_w3(1) + k)                       &
                                  - exner(map_w3(1) + k - 1) )                 &
                 + phi(map_w3(1) + k) - phi(map_w3(1) + k - 1) )               &
               / (height_w3(map_w3(1) + k) - height_w3(map_w3(1) + k - 1))
  end do
  if ( nlayers == 1 ) then
    pz(0) = 0.0_r_def
  else
    pz(0) = pz_face(1)
    pz(nlayers - 1) = pz_face(nlayers - 1)
    do k = 1, nlayers - 2
      pz(k) = 0.5_r_def * (pz_face(k) + pz_face(k + 1))
    end do
  end if

  ! Horizontal pressure gradient by least squares over the neighbours
  px(:) = 0.0_r_def
  py(:) = 0.0_r_def
  do df = 1, 4
    if ( has_neighbour(df) ) then
      do k = 0, nlayers - 1
        nb_w3 = exner_smap(1, 2, df) + k
        theta_f = 0.25_r_def * ( theta_v(map_wt(1) + k)                        &
                               + theta_v(map_wt(2) + k)                        &
                               + theta_v(theta_smap(1, 2, df) + k)             &
                               + theta_v(theta_smap(2, 2, df) + k) )
        dz = height_w3(hw3_smap(1, 2, df) + k) - height_w3(map_w3(1) + k)
        d_val = cp * theta_f * (exner(nb_w3) - exner(map_w3(1) + k))           &
              + phi(phi_smap(1, 2, df) + k) - phi(map_w3(1) + k)               &
              - dz * pz(k)
        px(k) = px(k) + w_grad(1, df) * d_val / rad(k)
        py(k) = py(k) + w_grad(2, df) * d_val / rad(k)
      end do
    end if
  end do

  ! Slopes of the levels at the bottom of each cell and at the top of the column
  zx(:) = 0.0_r_def
  zy(:) = 0.0_r_def
  do df = 1, 4
    if ( has_neighbour(df) ) then
      do k = 0, nlayers
        dh = height_wth(hwt_smap(1, 2, df) + k) - height_wth(map_wt(1) + k)
        zx(k) = zx(k) + w_grad(1, df) * dh / rad_lev(k)
        zy(k) = zy(k) + w_grad(2, df) * dh / rad_lev(k)
      end do
    end if
  end do

  ! Wind at the cell centre
  do k = 0, nlayers - 1
    sum_x = 0.0_r_def
    sum_y = 0.0_r_def
    do df = 1, 4
      un = sigma(df) * u(map_w2(df) + k) / dA_at_w2(map_w2(df) + k)
      sum_x = sum_x + w_wind(1, df) * un
      sum_y = sum_y + w_wind(2, df) * un
    end do
    uc(k) = sum_x
    vc(k) = sum_y
    w_bot = u(map_w2(5) + k) / dA_at_w2(map_w2(5) + k)                         &
          + zx(k) * uc(k) + zy(k) * vc(k)
    w_top = u(map_w2(6) + k) / dA_at_w2(map_w2(6) + k)                         &
          + zx(k + 1) * uc(k) + zy(k + 1) * vc(k)
    wc(k) = 0.5_r_def * (w_bot + w_top)
  end do

  ! Coriolis acceleration -2 Omega x u
  if ( rotating ) then
    do k = 0, nlayers - 1
      ax(k) = -(omega(2) * wc(k) - omega(3) * vc(k))
      ay(k) = -(omega(3) * uc(k) - omega(1) * wc(k))
      az(k) = -(omega(1) * vc(k) - omega(2) * uc(k))
    end do
  else
    ax(:) = 0.0_r_def
    ay(:) = 0.0_r_def
    az(:) = 0.0_r_def
  end if

  ! Side faces: compact pressure gradient, correction for the non-orthogonal
  ! and sloping neighbour, and Coriolis
  do df = 1, 4
    if ( has_neighbour(df) ) then
      do k = 0, nlayers - 1
        dz_nb = height_w3(hw3_smap(1, 2, df) + k) - height_w3(map_w3(1) + k)
        theta_f = 0.25_r_def * ( theta_v(map_wt(1) + k)                        &
                               + theta_v(map_wt(2) + k)                        &
                               + theta_v(theta_smap(1, 2, df) + k)             &
                               + theta_v(theta_smap(2, 2, df) + k) )
        d_val = cp * theta_f * exner(map_w3(1) + k) + phi(map_w3(1) + k)
        flux = coeff(df) / rad(k) * d_val                                      &
             - 0.5_r_def * ( correction(1, df) * px(k)                         &
                           + correction(2, df) * py(k)                         &
                           - dz_nb * coeff(df) / rad(k) * pz(k) )              &
             + 0.5_r_def * ( normal(1, df) * ax(k) + normal(2, df) * ay(k) )
        u_pred(map_w2(df) + k) = u_pred(map_w2(df) + k)                        &
            + alpha_dt * sigma(df) * dA_at_w2(map_w2(df) + k) * flux
      end do
    end if
  end do

  ! Vertical faces: compact pressure gradient, slope correction and Coriolis.
  ! The outward normal is up at the top face and down at the bottom face. The
  ! boundary faces are not changed
  do k = 1, nlayers - 1
    dz = height_w3(map_w3(1) + k) - height_w3(map_w3(1) + k - 1)
    flux = x_bot(k) / dz                                                       &
         - 0.5_r_def * (zx(k) * px(k) + zy(k) * py(k))                         &
         - 0.5_r_def * (az(k) - zx(k) * ax(k) - zy(k) * ay(k))
    u_pred(map_w2(5) + k) = u_pred(map_w2(5) + k)                              &
        + alpha_dt * sigma(5) * dA_at_w2(map_w2(5) + k) * flux
  end do
  do k = 0, nlayers - 2
    dz = height_w3(map_w3(1) + k + 1) - height_w3(map_w3(1) + k)
    flux = x_top(k) / dz                                                       &
         + 0.5_r_def * (zx(k + 1) * px(k) + zy(k + 1) * py(k))                 &
         + 0.5_r_def * (az(k) - zx(k + 1) * ax(k) - zy(k + 1) * ay(k))
    u_pred(map_w2(6) + k) = u_pred(map_w2(6) + k)                              &
        + alpha_dt * sigma(6) * dA_at_w2(map_w2(6) + k) * flux
  end do

end subroutine fd_wind_predictor_code

end module fd_wind_predictor_kernel_mod
