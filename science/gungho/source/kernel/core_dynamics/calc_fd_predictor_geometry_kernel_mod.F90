!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------

!> @brief   Computes horizontal geometry constants for the finite difference
!! wind predictor.
!> @details For lowest-order meshes with native coordinates, where the
!! horizontal coordinates do not vary with height. Everything is
!! expressed in a local orthonormal frame (east, north, up) at the
!! cell centre. Lengths are for a unit-radius sphere, or are physical
!! lengths for planar geometry. The values are stored in a 2D
!! multi-data field with fd_geom_ndata values per column, ordered as:
!! - for each face df=1,4 (W,S,E,N), a block of 9 values at
!! (df-1)*9 + (1..9):
!! 1-2  outward horizontal unit normal
!! 3    compact coefficient 1/(d e.n), where d is the distance to
!! the neighbouring cell centre and e the unit vector to it
!! 4-5  horizontal part of n - e/(e.n)
!! 6-7  least-squares weights for the horizontal gradient
!! 8-9  least-squares weights for the horizontal wind from the
!! outward normal wind components of the 4 faces
!! - 37-39 twice the planet rotation vector in the local frame
!! - 40    horizontal cell area
!! Faces with no neighbouring cell (e.g. limited-area edges) have
!! zeros for values 3-7.

module calc_fd_predictor_geometry_kernel_mod

use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_SCALAR, GH_REAL, GH_INTEGER,         &
                             GH_READ, GH_WRITE,                                &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_3,                        &
                             STENCIL, CROSS2D, CELL_COLUMN
use constants_mod,     only: r_def, i_def, l_def
use fs_continuity_mod, only: Wchi
use kernel_mod,        only: kernel_type

implicit none

private

!> Number of values stored in each column
integer(kind=i_def), parameter, public :: fd_geom_ndata = 40

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!! Psy layer.
type, public, extends(kernel_type) :: calc_fd_predictor_geometry_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                          &
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wchi, STENCIL(CROSS2D)),      &
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,    &
               STENCIL(CROSS2D)),                                              &
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                               &
      arg_type(GH_SCALAR,  GH_REAL,    GH_READ),                               &
      arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                                &
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: calc_fd_predictor_geometry_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: calc_fd_predictor_geometry_code

contains

!> @brief Computes the geometry constants for one column
!! @param[in]     nlayers        Number of layers
!! @param[in,out] fd_geom        Geometry constants, in a 2D multi-data field
!! @param[in]     chi_1          1st coordinate field
!! @param[in]     chi_2          2nd coordinate field
!! @param[in]     chi_3          3rd coordinate field
!! @param[in]     chi_smap_size  Size of the chi stencil map for each branch
!! @param[in]     chi_smap_max   Max size of the chi stencil map
!! @param[in]     chi_smap       Chi stencil map
!! @param[in]     panel_id       Field giving the ID for mesh panels
!! @param[in]     pid_smap_size  Size of the panel ID stencil map for each
!!                               branch
!! @param[in]     pid_smap_max   Max size of the panel ID stencil map
!! @param[in]     pid_smap       Panel ID stencil map
!! @param[in]     geometry       Integer indicating the domain geometry
!! @param[in]     omega          Planet angular velocity
!! @param[in]     f_lat          F-plane latitude
!! @param[in]     ndf_w3         Num DoFs per cell for W3
!! @param[in]     undf_w3        Num DoFs in this partition for W3
!! @param[in]     map_w3         Map of base cell DoFs for W3
!! @param[in]     ndf_chi        Num DoFs per cell for Wchi
!! @param[in]     undf_chi       Num DoFs in this partition for Wchi
!! @param[in]     map_chi        Map of base cell DoFs for Wchi
!! @param[in]     ndf_pid        Num DoFs per cell for panel ID
!! @param[in]     undf_pid       Num DoFs in this partition for panel ID
!! @param[in]     map_pid        Map of base cell DoFs for panel ID
subroutine calc_fd_predictor_geometry_code(                                    &
    nlayers,                                                                   &
    fd_geom,                                                                   &
    chi_1, chi_2, chi_3,                                                       &
    chi_smap_size,                                                             &
    chi_smap_max,                                                              &
    chi_smap,                                                                  &
    panel_id,                                                                  &
    pid_smap_size,                                                             &
    pid_smap_max,                                                              &
    pid_smap,                                                                  &
    geometry,                                                                  &
    omega,                                                                     &
    f_lat,                                                                     &
    ndf_w3, undf_w3, map_w3,                                                   &
    ndf_chi, undf_chi, map_chi,                                                &
    ndf_pid, undf_pid, map_pid                                                 &
)

  use base_mesh_config_mod,  only: geometry_spherical
  use sci_chi_transform_mod, only: chi2xyz

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: geometry
  integer(kind=i_def), intent(in) :: chi_smap_max
  integer(kind=i_def), intent(in) :: pid_smap_max

  real(kind=r_def), intent(inout) :: fd_geom(undf_w3)
  real(kind=r_def), intent(in)    :: chi_1(undf_chi)
  real(kind=r_def), intent(in)    :: chi_2(undf_chi)
  real(kind=r_def), intent(in)    :: chi_3(undf_chi)
  real(kind=r_def), intent(in)    :: panel_id(undf_pid)
  real(kind=r_def), intent(in)    :: omega
  real(kind=r_def), intent(in)    :: f_lat

  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in) :: chi_smap_size(4)
  integer(kind=i_def), intent(in) :: chi_smap(ndf_chi, chi_smap_max, 4)
  integer(kind=i_def), intent(in) :: pid_smap_size(4)
  integer(kind=i_def), intent(in) :: pid_smap(ndf_pid, pid_smap_max, 4)

  ! The corners of a cell at the base of the column are chi dofs 1 to 4 (SW,
  ! SE, NW, NE). Each face has two of them, and the neighbouring cell has the
  ! matching corner in position nb_corner_1
  integer(kind=i_def), parameter :: face_corner_1(4) = (/ 1, 1, 2, 3 /)
  integer(kind=i_def), parameter :: face_corner_2(4) = (/ 3, 2, 4, 4 /)
  integer(kind=i_def), parameter :: nb_corner_1(4)   = (/ 2, 3, 1, 1 /)
  real(kind=r_def),    parameter :: min_determinant  = 1.0e-12_r_def
  real(kind=r_def),    parameter :: min_alignment    = 1.0e-6_r_def

  ! Internal variables
  integer(kind=i_def) :: df, cell, corner_index, ipanel_cell, base
  integer(kind=i_def) :: cell_dofs(4)
  logical(kind=l_def) :: has_neighbour(4)

  real(kind=r_def) :: corner(3, 4, 0:4)
  real(kind=r_def) :: centre(3), centre_nb(3), up(3), east(3), north(3)
  real(kind=r_def) :: edge_vec(3), face_mid(3), separation(3)
  real(kind=r_def) :: omega_global(3), omega_local(3)
  real(kind=r_def) :: normal_h(2, 4), separation_h(2, 4), direction_h(2, 4)
  real(kind=r_def) :: correction_h(2, 4), distance(4), coefficient(4)
  real(kind=r_def) :: weights_gradient(2, 4), weights_wind(2, 4)
  real(kind=r_def) :: corner_h(2, 4), area
  real(kind=r_def) :: normal_matrix(2, 2), inverse_matrix(2, 2), determinant
  real(kind=r_def) :: edge_h(2), face_mid_h(2), length, dir_dot_normal
  real(kind=r_def) :: x, y, z, radius

  ! Corners of this cell (cell 0) and of each neighbouring cell. The height
  ! only scales a point, which is removed by normalising
  corner(:, :, :) = 0.0_r_def
  has_neighbour(:) = .false.
  do cell = 0, 4
    if ( cell == 0 ) then
      cell_dofs(:) = map_chi(1:4)
      ipanel_cell = int(panel_id(map_pid(1)), i_def)
    else if ( chi_smap_size(cell) > 1 ) then
      cell_dofs(:) = chi_smap(1:4, 2, cell)
      ipanel_cell = int(panel_id(pid_smap(1, 2, cell)), i_def)
      has_neighbour(cell) = .true.
    else
      cycle
    end if
    do corner_index = 1, 4
      if ( geometry == geometry_spherical ) then
        call chi2xyz( chi_1(cell_dofs(corner_index)),                          &
                      chi_2(cell_dofs(corner_index)),                          &
                      chi_3(cell_dofs(corner_index)),                          &
                      ipanel_cell, x, y, z )
        radius = sqrt(x * x + y * y + z * z)
        corner(:, corner_index, cell) = (/ x, y, z /) / radius
      else
        corner(:, corner_index, cell) = (/ chi_1(cell_dofs(corner_index)),     &
                                           chi_2(cell_dofs(corner_index)),     &
                                           0.0_r_def /)
      end if
    end do
  end do

  ! Cell centre and local frame
  centre(:) = 0.25_r_def * sum(corner(:, :, 0), dim=2)
  if ( geometry == geometry_spherical ) then
    up(:) = centre(:) / sqrt(dot_product(centre, centre))
    centre(:) = up(:)
  else
    up(:) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
  end if
  edge_vec(:) = 0.5_r_def * ( corner(:, 2, 0) + corner(:, 4, 0)                &
      - corner(:, 1, 0) - corner(:, 3, 0) )
  edge_vec(:) = edge_vec(:) - dot_product(edge_vec, up) * up(:)
  east(:) = edge_vec(:) / sqrt(dot_product(edge_vec, edge_vec))
  north(:) = (/ up(2) * east(3) - up(3) * east(2),                             &
                up(3) * east(1) - up(1) * east(3),                             &
                up(1) * east(2) - up(2) * east(1) /)

  ! Horizontal area from the corners projected into the local frame
  do df = 1, 4
    corner_h(1, df) = dot_product(corner(:, df, 0), east)
    corner_h(2, df) = dot_product(corner(:, df, 0), north)
  end do
  area = 0.5_r_def * abs( (corner_h(1, 4) - corner_h(1, 1))                    &
      * (corner_h(2, 3) - corner_h(2, 2))                                      &
      - (corner_h(2, 4) - corner_h(2, 1))                                      &
      * (corner_h(1, 3) - corner_h(1, 2)) )

  normal_h(:, :) = 0.0_r_def
  separation_h(:, :) = 0.0_r_def
  direction_h(:, :) = 0.0_r_def
  correction_h(:, :) = 0.0_r_def
  distance(:) = 0.0_r_def
  coefficient(:) = 0.0_r_def

  do df = 1, 4

    ! Outward horizontal unit normal of the face
    face_mid(:) = 0.5_r_def * ( corner(:, face_corner_1(df), 0)                &
        + corner(:, face_corner_2(df), 0) ) - centre(:)
    face_mid_h(1) = dot_product(face_mid, east)
    face_mid_h(2) = dot_product(face_mid, north)
    edge_vec(:) = corner(:, face_corner_2(df), 0)                              &
        - corner(:, face_corner_1(df), 0)
    edge_h(1) = dot_product(edge_vec, east)
    edge_h(2) = dot_product(edge_vec, north)
    normal_h(:, df) = (/ edge_h(2), -edge_h(1) /)
    if ( dot_product(normal_h(:, df), face_mid_h) < 0.0_r_def ) then
      normal_h(:, df) = -normal_h(:, df)
    end if
    length = sqrt(dot_product(normal_h(:, df), normal_h(:, df)))
    normal_h(:, df) = normal_h(:, df) / length

    ! Separation of the cell centres, in the local frame
    if ( has_neighbour(df) ) then
      centre_nb(:) = 0.25_r_def * sum(corner(:, :, df), dim=2)
      if ( geometry == geometry_spherical ) then
        centre_nb(:) = centre_nb(:) / sqrt(dot_product(centre_nb, centre_nb))
      else
        ! Planar meshes can be periodic, so move the neighbour to be
        ! adjacent to this cell
        centre_nb(:) = centre_nb(:) + corner(:, face_corner_1(df), 0)          &
            - corner(:, nb_corner_1(df), df)
      end if
      separation(:) = centre_nb(:) - centre(:)
      separation_h(1, df) = dot_product(separation, east)
      separation_h(2, df) = dot_product(separation, north)
      distance(df) = sqrt(dot_product(separation_h(:, df),                     &
                                      separation_h(:, df)))
      direction_h(:, df) = separation_h(:, df) / distance(df)
      dir_dot_normal = dot_product(direction_h(:, df), normal_h(:, df))
      has_neighbour(df) = ( dir_dot_normal > min_alignment )
      if ( has_neighbour(df) ) then
        coefficient(df) = 1.0_r_def / (distance(df) * dir_dot_normal)
        correction_h(:, df) = normal_h(:, df)                                  &
            - direction_h(:, df) / dir_dot_normal
      end if
    end if

  end do

  ! Least-squares weights for the horizontal gradient from differences with
  ! the available neighbours
  weights_gradient(:, :) = 0.0_r_def
  normal_matrix(:, :) = 0.0_r_def
  do df = 1, 4
    if ( has_neighbour(df) ) then
      normal_matrix(1, 1) = normal_matrix(1, 1)                                &
          + separation_h(1, df) * separation_h(1, df)
      normal_matrix(1, 2) = normal_matrix(1, 2)                                &
          + separation_h(1, df) * separation_h(2, df)
      normal_matrix(2, 2) = normal_matrix(2, 2)                                &
          + separation_h(2, df) * separation_h(2, df)
    end if
  end do
  normal_matrix(2, 1) = normal_matrix(1, 2)
  determinant = normal_matrix(1, 1) * normal_matrix(2, 2)                      &
      - normal_matrix(1, 2) * normal_matrix(2, 1)
  if ( determinant > min_determinant                                           &
      * (normal_matrix(1, 1) + normal_matrix(2, 2))**2 ) then
    inverse_matrix(1, 1) = normal_matrix(2, 2) / determinant
    inverse_matrix(1, 2) = -normal_matrix(1, 2) / determinant
    inverse_matrix(2, 1) = -normal_matrix(2, 1) / determinant
    inverse_matrix(2, 2) = normal_matrix(1, 1) / determinant
    do df = 1, 4
      if ( has_neighbour(df) ) then
        weights_gradient(:, df) = matmul(inverse_matrix, separation_h(:, df))
      end if
    end do
  end if

  ! Least-squares weights for the horizontal wind from the outward normal
  ! wind components of all 4 faces
  normal_matrix(:, :) = 0.0_r_def
  do df = 1, 4
    normal_matrix(1, 1) = normal_matrix(1, 1) +                                &
                          normal_h(1, df) * normal_h(1, df)
    normal_matrix(1, 2) = normal_matrix(1, 2) +                                &
                          normal_h(1, df) * normal_h(2, df)
    normal_matrix(2, 2) = normal_matrix(2, 2) +                                &
                          normal_h(2, df) * normal_h(2, df)
  end do
  normal_matrix(2, 1) = normal_matrix(1, 2)
  determinant = normal_matrix(1, 1) * normal_matrix(2, 2)                      &
              - normal_matrix(1, 2) * normal_matrix(2, 1)
  inverse_matrix(1, 1) = normal_matrix(2, 2) / determinant
  inverse_matrix(1, 2) = -normal_matrix(1, 2) / determinant
  inverse_matrix(2, 1) = -normal_matrix(2, 1) / determinant
  inverse_matrix(2, 2) = normal_matrix(1, 1) / determinant
  do df = 1, 4
    weights_wind(:, df) = matmul(inverse_matrix, normal_h(:, df))
  end do

  ! Twice the planet rotation vector, in the local frame
  if ( geometry == geometry_spherical ) then
    omega_global(:) = (/ 0.0_r_def, 0.0_r_def, 2.0_r_def * omega /)
  else
    omega_global(:) = (/ 0.0_r_def,                                            &
                         2.0_r_def * omega * cos(f_lat),                       &
                         2.0_r_def * omega * sin(f_lat) /)
  end if
  omega_local(1) = dot_product(omega_global, east)
  omega_local(2) = dot_product(omega_global, north)
  omega_local(3) = dot_product(omega_global, up)

  ! Store the values
  do df = 1, 4
    base = map_w3(1) + (df - 1) * 9 - 1
    fd_geom(base + 1) = normal_h(1, df)
    fd_geom(base + 2) = normal_h(2, df)
    fd_geom(base + 3) = coefficient(df)
    fd_geom(base + 4) = correction_h(1, df)
    fd_geom(base + 5) = correction_h(2, df)
    fd_geom(base + 6) = weights_gradient(1, df)
    fd_geom(base + 7) = weights_gradient(2, df)
    fd_geom(base + 8) = weights_wind(1, df)
    fd_geom(base + 9) = weights_wind(2, df)
  end do
  base = map_w3(1) + 36 - 1
  fd_geom(base + 1) = omega_local(1)
  fd_geom(base + 2) = omega_local(2)
  fd_geom(base + 3) = omega_local(3)
  fd_geom(base + 4) = area

end subroutine calc_fd_predictor_geometry_code

end module calc_fd_predictor_geometry_kernel_mod
