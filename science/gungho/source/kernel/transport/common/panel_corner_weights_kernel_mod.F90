!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Compute the interpolation weights and indices to remap a W3 or Wtheta
!!        scalar field from the standard cubed-sphere mesh onto the
!!        neighbouring panels around a cubed-sphere corner.
!!
!> @details This kernel extends the edge-remap logic of
!!          panel_edge_weights_kernel_mod to the case where a target cell is
!!          near both an x-edge and a y-edge simultaneously. At such a corner, a
!!          field must be remapped using the coordinate prolongation computed by
!!          panel_corner_coords_kernel_mod: the x-type correction is evaluated in
!!          the anti-clockwise corner coordinate pair (reprojected into the x-
!!          crossing neighbour), while the y-type correction is evaluated in the
!!          clockwise corner coordinate pair (reprojected into the y-crossing
!!          neighbour).
!!
!!          The remapping is a 2D bilinear interpolation using four cells of a
!!          REGION stencil. The target and the centres of all stencil cells are
!!          converted to the owned panel's alpha-beta frame. The four cells
!!          forming the smallest quadrilateral enclosing the target are chosen,
!!          and the weights come from the inverse bilinear map of that
!!          quadrilateral. Where no such quadrilateral exists (e.g. in the
!!          missing quadrant at a cube corner), three cells with barycentric
!!          weights are used, or else the nearest cell. The indices are
!!          positions in the REGION stencil.
module panel_corner_weights_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type, func_type,                          &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER,                          &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_9,                    &
                                 ANY_SPACE_9,                                  &
                                 GH_BASIS, GH_EVALUATOR,                       &
                                 CELL_COLUMN, STENCIL, REGION
use constants_mod,         only: r_tran, r_def, i_def, l_def, LARGE_REAL_POSITIVE
use reference_element_mod, only: W, S, N, E
use sci_chi_transform_mod, only: chi2xyz, get_to_rotate, get_to_stretch,         &
                                 get_inverse_mesh_rotation_matrix,                  &
                                 get_stretch_factor
use coord_transform_mod,   only: alphabetar2xyz, xyz2alphabetar,                 &
                                 inverse_schmidt_transform_xyz

! Configuration modules
use base_mesh_config_mod,      only: geometry, topology
use finite_element_config_mod, only: coord_system
use planet_config_mod,         only: scaled_radius

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_corner_weights_kernel_type
  private
  type(arg_type) :: meta_args(10) = (/                                         &
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
      arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
      arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9,                 &
                                                       STENCIL(REGION)),      &
      arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),  &
      arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),  &
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,   &
                                                       STENCIL(REGION)),      &
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                               &
  /)
  type(func_type) :: meta_funcs(2) = (/                                        &
      func_type(ANY_SPACE_9, GH_BASIS),                                        &
      func_type(ANY_DISCONTINUOUS_SPACE_9, GH_BASIS)                           &
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: panel_corner_weights_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_corner_weights_code

contains

!> @brief Compute the interpolation weights and indices to remap a W3 or Wtheta
!!        scalar field from the standard cubed-sphere mesh to the two panels
!!        meeting at a cubed-sphere corner.
!> @param[in]     nlayers                 Number of layers
!> @param[in,out] panel_corner_weights_x  Remapping interpolation weights for x
!> @param[in,out] panel_corner_weights_y  Remapping interpolation weights for y
!> @param[in,out] panel_corner_indices_x  Remapping interpolation indices for x
!> @param[in,out] panel_corner_indices_y  Remapping interpolation indices for y
!> @param[in]     alpha                   Alpha coordinate
!> @param[in]     beta                    Beta coordinate
!> @param[in]     height                  Height coordinate
!> @param[in]     wx_stencil_size         Size of the Wchi stencil (num cells)
!> @param[in]     wx_stencil_map          DoF map for the Wchi stencil
!> @param[in]     alpha_anti              2D alpha coord extended around corner
!> @param[in]     beta_anti               2D beta coord extended around corner
!> @param[in]     alpha_clock             2D alpha coord extended around corner
!> @param[in]     beta_clock              2D beta coord extended around corner
!> @param[in]     panel_id                ID of cubed sphere panel for each cell
!> @param[in]     pid_stencil_size        Size of panel ID stencil (num cells)
!> @param[in]     pid_stencil_map         DoF map for the panel id stencil
!> @param[in]     panel_edge_dist_W       Distance to the West panel edge
!> @param[in]     panel_edge_dist_E       Distance to the East panel edge
!> @param[in]     panel_edge_dist_S       Distance to the South panel edge
!> @param[in]     panel_edge_dist_N       Distance to the North panel edge
!> @param[in]     remap_depth             Depth to which the remapping may happen
!> @param[in]     ndf_ww                  Num DoFs per cell for remapping weights
!> @param[in]     undf_ww                 Num DoFs for this partition for weights
!> @param[in]     map_ww                  DoF map for remapping weights
!> @param[in]     ndf_wx                  Num DoFs per cell for coord fields
!> @param[in]     undf_wx                 Num DoFs for this partition for coords
!> @param[in]     map_wx                  DoF map for coord fields
!> @param[in]     basis_wx                Basis functions of coordinate fields,
!!                                      evaluated at nodes of scalar fields
!> @param[in]     ndf_wx_2d               Num DoFs per cell for 2D coord fields
!> @param[in]     undf_wx_2d              Num DoFs in this partition for 2D coords
!> @param[in]     map_wx_2d               DoF map for 2D coord fields
!> @param[in]     basis_wx_2d             Basis functions of 2D coordinate fields,
!!                                      evaluated at nodes of scalar fields
!> @param[in]     ndf_pid                 Num DoFs per cell for panel ID field
!> @param[in]     undf_pid                Num DoFs for this partition for panel ID
!> @param[in]     map_pid                 DoF map for panel ID field
subroutine panel_corner_weights_code( nlayers,                                   &
                                     panel_corner_weights_x,                    &
                                     panel_corner_weights_y,                    &
                                     panel_corner_indices_x,                    &
                                     panel_corner_indices_y,                    &
                                     alpha, beta, height,                      &
                                     wx_stencil_size,                          &
                                     wx_stencil_map,                           &
                                     alpha_anti, beta_anti,                    &
                                     alpha_clock, beta_clock,                  &
                                     panel_id,                                 &
                                     pid_stencil_size,                         &
                                     pid_stencil_map,                          &
                                     panel_edge_dist_W,                        &
                                     panel_edge_dist_E,                        &
                                     panel_edge_dist_S,                        &
                                     panel_edge_dist_N,                        &
                                     remap_depth,                              &
                                     ndf_ww, undf_ww, map_ww,                  &
                                     ndf_wx, undf_wx, map_wx, basis_wx,        &
                                     ndf_wx_2d, undf_wx_2d, map_wx_2d,         &
                                     basis_wx_2d,                              &
                                     ndf_pid, undf_pid, map_pid                &
  )

  use panel_edge_support_mod, only: panel_neighbour

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: wx_stencil_size
  integer(kind=i_def), intent(in)    :: pid_stencil_size
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: remap_depth

  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wx_stencil_map(ndf_wx, wx_stencil_size)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid,                 &
                                                      pid_stencil_size)
  real(kind=r_def),    intent(in)    :: basis_wx(1, ndf_wx, ndf_ww)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_ww)

  real(kind=r_tran),   intent(inout) :: panel_corner_weights_x(undf_ww)
  real(kind=r_tran),   intent(inout) :: panel_corner_weights_y(undf_ww)
  integer(kind=i_def), intent(inout) :: panel_corner_indices_x(undf_ww)
  integer(kind=i_def), intent(inout) :: panel_corner_indices_y(undf_ww)
  real(kind=r_def),    intent(in)    :: alpha(undf_wx)
  real(kind=r_def),    intent(in)    :: beta(undf_wx)
  real(kind=r_def),    intent(in)    :: height(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_anti(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_anti(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: alpha_clock(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_clock(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)

  ! Internal variables
  integer(kind=i_def) :: pid_idx, owned_panel, swapped_panel_x, swapped_panel_y
  real(kind=r_def)    :: cell_alpha(wx_stencil_size)
  real(kind=r_def)    :: cell_beta(wx_stencil_size)

  ! Panel id for this column
  pid_idx = map_pid(1)

  swapped_panel_x = 0_i_def
  swapped_panel_y = 0_i_def

  ! Determine which corner-neighbouring panels are relevant. This uses the same
  ! logic and the same depth (the mesh halo depth) as
  ! panel_corner_coords_kernel_mod, so that the corner coordinates used below
  ! were generated for the same pair of neighbouring panels.
  owned_panel = int(panel_id(pid_idx), i_def)

  if (ABS(panel_edge_dist_W(pid_idx)) <= remap_depth) then
    swapped_panel_x = panel_neighbour(owned_panel, W)
  end if
  if (ABS(panel_edge_dist_E(pid_idx)) <= remap_depth .and.                        &
      ABS(panel_edge_dist_E(pid_idx)) < ABS(panel_edge_dist_W(pid_idx))) then
    swapped_panel_x = panel_neighbour(owned_panel, E)
  end if

  if (ABS(panel_edge_dist_S(pid_idx)) <= remap_depth) then
    swapped_panel_y = panel_neighbour(owned_panel, S)
  end if
  if (ABS(panel_edge_dist_N(pid_idx)) <= remap_depth .and.                        &
      ABS(panel_edge_dist_N(pid_idx)) < ABS(panel_edge_dist_S(pid_idx))) then
    swapped_panel_y = panel_neighbour(owned_panel, N)
  end if

  if (swapped_panel_x /= 0_i_def .and. swapped_panel_y /= 0_i_def) then
    ! Centres of all stencil cells in the owned panel's alpha-beta frame
    call stencil_cell_positions(                                               &
            owned_panel, alpha, beta, panel_id,                                &
            wx_stencil_size, wx_stencil_map, pid_stencil_map,                  &
            basis_wx(1,:,1), cell_alpha, cell_beta,                            &
            ndf_wx, undf_wx, ndf_pid, undf_pid                                 &
    )

    ! x-type remap around the corner: evaluate the target in the
    ! anti-clockwise corner coordinates (reprojected into the x-crossing panel)
    call panel_corner_weights_2d(                                              &
            panel_corner_weights_x, panel_corner_indices_x,                    &
            owned_panel, swapped_panel_x,                                      &
            alpha_anti, beta_anti,                                             &
            wx_stencil_size, cell_alpha, cell_beta,                            &
            ndf_ww, undf_ww, map_ww,                                           &
            ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d                      &
    )

    ! y-type remap around the corner: evaluate the target in the clockwise
    ! corner coordinates (reprojected into the y-crossing panel)
    call panel_corner_weights_2d(                                              &
            panel_corner_weights_y, panel_corner_indices_y,                    &
            owned_panel, swapped_panel_y,                                      &
            alpha_clock, beta_clock,                                           &
            wx_stencil_size, cell_alpha, cell_beta,                            &
            ndf_ww, undf_ww, map_ww,                                           &
            ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d                      &
    )
  else
    ! Not near a panel corner: identity remap in both directions.
    panel_corner_weights_x(map_ww(1):map_ww(1)+ndf_ww-1) = 0.0_r_tran
    panel_corner_weights_y(map_ww(1):map_ww(1)+ndf_ww-1) = 0.0_r_tran
    panel_corner_indices_x(map_ww(1):map_ww(1)+ndf_ww-1) = 1
    panel_corner_indices_y(map_ww(1):map_ww(1)+ndf_ww-1) = 1
    panel_corner_weights_x(map_ww(1)) = 1.0_r_tran
    panel_corner_weights_y(map_ww(1)) = 1.0_r_tran
  end if

end subroutine panel_corner_weights_code

!> @brief Private routine to find the centres of all cells in the REGION
!!        stencil, expressed in the owned panel's alpha-beta frame.
!> @param[in]  owned_panel      Panel ID of this column
!> @param[in]  alpha            Alpha coordinate
!> @param[in]  beta             Beta coordinate
!> @param[in]  panel_id         ID of cubed sphere panel for each cell
!> @param[in]  stencil_size     Size of the stencil (num cells)
!> @param[in]  wx_stencil_map   DoF map for the Wchi stencil
!> @param[in]  pid_stencil_map  DoF map for the panel id stencil
!> @param[in]  basis            Basis functions of the coordinate field at the
!!                              cell centre
!> @param[out] cell_alpha       Alpha of each stencil cell centre
!> @param[out] cell_beta        Beta of each stencil cell centre
!> @param[in]  ndf_wx           Num DoFs per cell for coord fields
!> @param[in]  undf_wx          Num DoFs for this partition for coords
!> @param[in]  ndf_pid          Num DoFs per cell for panel ID field
!> @param[in]  undf_pid         Num DoFs for this partition for panel ID
subroutine stencil_cell_positions( owned_panel, alpha, beta, panel_id,         &
                                   stencil_size, wx_stencil_map,               &
                                   pid_stencil_map, basis,                     &
                                   cell_alpha, cell_beta,                      &
                                   ndf_wx, undf_wx, ndf_pid, undf_pid )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: owned_panel
  integer(kind=i_def), intent(in)    :: stencil_size
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: wx_stencil_map(ndf_wx, stencil_size)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid, stencil_size)
  real(kind=r_def),    intent(in)    :: basis(ndf_wx)
  real(kind=r_def),    intent(in)    :: alpha(undf_wx)
  real(kind=r_def),    intent(in)    :: beta(undf_wx)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(inout) :: cell_alpha(stencil_size)
  real(kind=r_def),    intent(inout) :: cell_beta(stencil_size)

  ! Internal variables
  real(kind=r_def), parameter :: unit_radius = 1.0_r_def
  integer(kind=i_def)         :: n, df, cell_panel
  real(kind=r_def)    :: xyz(3), h_dummy
  real(kind=r_def)    :: alpha_owned, beta_owned
  real(kind=r_def)    :: stretch_factor
  real(kind=r_def)    :: inverse_rot_matrix(3,3)
  logical(kind=l_def) :: to_rotate, to_stretch

  to_rotate = get_to_rotate()
  to_stretch = get_to_stretch()
  if (to_rotate) then
    inverse_rot_matrix = get_inverse_mesh_rotation_matrix()
  end if
  if (to_stretch) then
    stretch_factor = get_stretch_factor()
  end if

  do n = 1, stencil_size
    cell_alpha(n) = 0.0_r_def
    cell_beta(n) = 0.0_r_def

    if (pid_stencil_map(1, n) <= 0 .or. wx_stencil_map(1, n) <= 0) then
      ! Missing stencil cell: place it where it will never be selected
      cell_alpha(n) = LARGE_REAL_POSITIVE
      cell_beta(n) = LARGE_REAL_POSITIVE
    else
      cell_panel = int(panel_id(pid_stencil_map(1, n)), i_def)
      do df = 1, ndf_wx
        call chi2xyz(                                                          &
                alpha(wx_stencil_map(df, n)), beta(wx_stencil_map(df, n)),     &
                unit_radius, cell_panel, geometry, topology, coord_system,     &
                scaled_radius, xyz(1), xyz(2), xyz(3)                          &
        )
        if (to_rotate) then
          xyz = matmul(inverse_rot_matrix, xyz)
        end if
        if (to_stretch) then
          xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
        end if
        call xyz2alphabetar(                                                   &
                xyz(1), xyz(2), xyz(3), owned_panel,                           &
                alpha_owned, beta_owned, h_dummy                               &
        )
        cell_alpha(n) = cell_alpha(n) + alpha_owned * basis(df)
        cell_beta(n) = cell_beta(n) + beta_owned * basis(df)
      end do
    end if
  end do

end subroutine stencil_cell_positions

!> @brief Private routine for calculation of corner-remapping weights and
!!        indices for a single remap direction (x or y).
!> @param[in,out] remap_weights   Remapping interpolation weights
!> @param[in,out] remap_indices   Remapping interpolation indices (positions in
!!                                the REGION stencil)
!> @param[in]     owned_panel     Panel ID of this column
!> @param[in]     swapped_panel   Panel in which the target coordinates are
!!                                expressed
!> @param[in]     alpha_ext       2D alpha coord extended around corner
!> @param[in]     beta_ext        2D beta coord extended around corner
!> @param[in]     stencil_size    Size of the stencil (num cells)
!> @param[in]     cell_alpha      Alpha of each stencil cell centre, owned frame
!> @param[in]     cell_beta       Beta of each stencil cell centre, owned frame
!> @param[in]     ndf_ww          Num DoFs per cell for remapping weights
!> @param[in]     undf_ww         Num DoFs for this partition for weights
!> @param[in]     map_ww          DoF map for remapping weights
!> @param[in]     ndf_wx_2d       Num DoFs per cell for 2D coord fields
!> @param[in]     undf_wx_2d      Num DoFs in this partition for 2D coords
!> @param[in]     map_wx_2d       DoF map for 2D coord fields
!> @param[in]     basis_wx_2d     Basis functions of 2D coordinate fields,
!!                                evaluated at nodes of scalar fields
subroutine panel_corner_weights_2d( remap_weights, remap_indices,              &
                                    owned_panel, swapped_panel,                &
                                    alpha_ext, beta_ext,                       &
                                    stencil_size, cell_alpha, cell_beta,       &
                                    ndf_ww, undf_ww, map_ww,                   &
                                    ndf_wx_2d, undf_wx_2d, map_wx_2d,          &
                                    basis_wx_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: owned_panel
  integer(kind=i_def), intent(in)    :: swapped_panel
  integer(kind=i_def), intent(in)    :: stencil_size
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_ww)
  real(kind=r_def),    intent(in)    :: alpha_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: cell_alpha(stencil_size)
  real(kind=r_def),    intent(in)    :: cell_beta(stencil_size)
  real(kind=r_tran),   intent(inout) :: remap_weights(undf_ww)
  integer(kind=i_def), intent(inout) :: remap_indices(undf_ww)

  ! Internal variables
  real(kind=r_def), parameter :: unit_radius = 1.0_r_def
  integer(kind=i_def)         :: df, donor_ids(4)
  real(kind=r_def)            :: abh(2), xyz(3), h_dummy
  real(kind=r_def)            :: donor_weights(4)

  ! Target point in the swapped-panel frame
  abh = 0.0_r_def
  do df = 1, ndf_wx_2d
    abh(1) = abh(1) + alpha_ext(map_wx_2d(df)) * basis_wx_2d(1, df, 1)
    abh(2) = abh(2) + beta_ext(map_wx_2d(df)) * basis_wx_2d(1, df, 1)
  end do

  ! Convert the target from the swapped-panel frame to the owned-panel frame,
  ! in which the stencil cell centres are expressed
  call alphabetar2xyz(abh(1), abh(2), unit_radius,                             &
                     swapped_panel, xyz(1), xyz(2), xyz(3))
  call xyz2alphabetar(xyz(1), xyz(2), xyz(3), owned_panel,                     &
                      abh(1), abh(2), h_dummy)

  call bilinear_donors(abh(1), abh(2), stencil_size, cell_alpha, cell_beta,    &
                       donor_ids, donor_weights)

  remap_indices(map_ww(1):map_ww(1) + 3) = donor_ids
  remap_weights(map_ww(1):map_ww(1) + 3) = real(donor_weights, r_tran)

end subroutine panel_corner_weights_2d

!> @brief Private routine to choose four stencil cells and bilinear weights
!!        for interpolating to a target point.
!> @details The cells forming the smallest convex quadrilateral that encloses
!!          the target are used, with weights from the inverse bilinear map of
!!          that quadrilateral. If there is none, the smallest enclosing
!!          triangle gives barycentric weights (the fourth weight is zero). If
!!          there is none, the nearest cell is used. Only the nearest few cells
!!          are considered.
!> @param[in]  target_a       Alpha of the target point
!> @param[in]  target_b       Beta of the target point
!> @param[in]  n_cells        Number of stencil cells
!> @param[in]  cell_a         Alpha of each stencil cell centre
!> @param[in]  cell_b         Beta of each stencil cell centre
!> @param[out] ids            Stencil positions of the four donor cells
!> @param[out] wgts           Weights of the four donor cells
subroutine bilinear_donors(target_a, target_b, n_cells, cell_a, cell_b,        &
                           ids, wgts)

  implicit none

  integer(kind=i_def), intent(in)    :: n_cells
  real(kind=r_def),    intent(in)    :: target_a, target_b
  real(kind=r_def),    intent(in)    :: cell_a(n_cells), cell_b(n_cells)
  integer(kind=i_def), intent(inout) :: ids(4)
  real(kind=r_def),    intent(inout) :: wgts(4)

  integer(kind=i_def), parameter :: n_near_max = 10
  real(kind=r_def),    parameter :: tol = 1.0e-10_r_def

  integer(kind=i_def) :: n_near, m, i, j, k, l
  integer(kind=i_def) :: near(n_near_max), cand(4), ordered(4)
  real(kind=r_def)    :: work(n_cells)
  real(kind=r_def)    :: area, signed_area, best_area, s, t
  real(kind=r_def)    :: lambda(3)
  logical(kind=l_def) :: found, inside

  work = (cell_a - target_a)**2 + (cell_b - target_b)**2
  n_near = MIN(n_near_max, n_cells)
  do m = 1, n_near
    near(m) = MINLOC(work, 1)
    work(near(m)) = LARGE_REAL_POSITIVE
  end do

  ! Default: nearest cell only
  ids = near(1)
  wgts = 0.0_r_def
  wgts(1) = 1.0_r_def

  ! Smallest enclosing convex quadrilateral
  found = .false.
  best_area = LARGE_REAL_POSITIVE
  do i = 1, n_near - 3
    do j = i + 1, n_near - 2
      do k = j + 1, n_near - 1
        do l = k + 1, n_near
          cand = (/ near(i), near(j), near(k), near(l) /)
          call enclosing_quad(target_a, target_b, cell_a(cand), cell_b(cand),  &
                              ordered, area, inside)
          if (inside .and. area < best_area) then
            best_area = area
            ids = cand(ordered)
            found = .true.
          end if
        end do
      end do
    end do
  end do

  if (found) then
    call inverse_bilinear(target_a, target_b, cell_a(ids),                     &
                         cell_b(ids), s, t)
    wgts(1) = (1.0_r_def - s) * (1.0_r_def - t)
    wgts(2) = s * (1.0_r_def - t)
    wgts(3) = s * t
    wgts(4) = (1.0_r_def - s) * t
    return
  end if

  ! Smallest enclosing triangle
  best_area = LARGE_REAL_POSITIVE
  do i = 1, n_near - 2
    do j = i + 1, n_near - 1
      do k = j + 1, n_near
        cand(1:3) = (/ near(i), near(j), near(k) /)
        signed_area = 0.5_r_def * cross_2d(                                    &
            cell_a(cand(2)) - cell_a(cand(1)),                                 &
            cell_b(cand(2)) - cell_b(cand(1)),                                 &
            cell_a(cand(3)) - cell_a(cand(1)),                                 &
            cell_b(cand(3)) - cell_b(cand(1)))
        area = ABS(signed_area)
        if (area > tol .and. area < best_area) then
          lambda(1) = 0.5_r_def * cross_2d(                                    &
              cell_a(cand(2)) - target_a, cell_b(cand(2)) - target_b,          &
              cell_a(cand(3)) - target_a, cell_b(cand(3)) - target_b)          &
              / signed_area
          lambda(2) = 0.5_r_def * cross_2d(                                    &
              cell_a(cand(3)) - target_a, cell_b(cand(3)) - target_b,          &
              cell_a(cand(1)) - target_a, cell_b(cand(1)) - target_b)          &
              / signed_area
          lambda(3) = 1.0_r_def - lambda(1) - lambda(2)
          if (ALL(lambda >= -tol)) then
            best_area = area
            ids = (/ cand(1), cand(2), cand(3), cand(1) /)
            wgts = (/ lambda(1), lambda(2), lambda(3), 0.0_r_def /)
          end if
        end if
      end do
    end do
  end do

end subroutine bilinear_donors

!> @brief Private routine to order four points anticlockwise and check that
!!        they form a convex quadrilateral enclosing the target point.
!> @param[in]  target_a  Alpha of the target point
!> @param[in]  target_b  Beta of the target point
!> @param[in]  pa        Alpha of the four points
!> @param[in]  pb        Beta of the four points
!> @param[out] ordered   Indices of the points in anticlockwise order
!> @param[out] area      Area of the quadrilateral
!> @param[out] inside    Whether the quadrilateral is convex and contains the
!!                       target
subroutine enclosing_quad(target_a, target_b, pa, pb, ordered, area, inside)

  implicit none

  real(kind=r_def),    intent(in)    :: target_a, target_b
  real(kind=r_def),    intent(in)    :: pa(4), pb(4)
  integer(kind=i_def), intent(inout) :: ordered(4)
  real(kind=r_def),    intent(inout) :: area
  logical(kind=l_def), intent(inout) :: inside

  real(kind=r_def), parameter :: tol = 1.0e-10_r_def
  integer(kind=i_def)         :: i, j, e, next, after, tmp
  real(kind=r_def)            :: centre_a, centre_b, angle(4), c

  centre_a = SUM(pa) / 4.0_r_def
  centre_b = SUM(pb) / 4.0_r_def
  do i = 1, 4
    angle(i) = ATAN2(pb(i) - centre_b, pa(i) - centre_a)
    ordered(i) = i
  end do

  ! Sort by angle about the centroid (insertion sort of four items)
  do i = 2, 4
    tmp = ordered(i)
    j = i - 1
    do while (j >= 1)
      if (angle(ordered(j)) <= angle(tmp)) exit
      ordered(j + 1) = ordered(j)
      j = j - 1
    end do
    ordered(j + 1) = tmp
  end do

  inside = .true.
  area = 0.0_r_def
  do e = 1, 4
    next = ordered(MOD(e, 4) + 1)
    after = ordered(MOD(e + 1, 4) + 1)
    i = ordered(e)
    ! Convexity: every turn is anticlockwise
    c = cross_2d(pa(next) - pa(i), pb(next) - pb(i),                           &
                 pa(after) - pa(next), pb(after) - pb(next))
    if (c <= tol) inside = .false.
    ! The target lies to the left of every edge
    c = cross_2d(pa(next) - pa(i), pb(next) - pb(i),                           &
                 target_a - pa(i), target_b - pb(i))
    if (c < -tol) inside = .false.
    area = area + 0.5_r_def * cross_2d(pa(i), pb(i), pa(next), pb(next))
  end do

end subroutine enclosing_quad

!> @brief Private routine to solve for the local coordinates (s,t) of a point
!!        inside a quadrilateral with anticlockwise corners P1 to P4, where the
!!        point is P1 + s(P2-P1) + t(P4-P1) + st(P1-P2+P3-P4).
!> @param[in]  target_a  Alpha of the target point
!> @param[in]  target_b  Beta of the target point
!> @param[in]  pa        Alpha of the four corners
!> @param[in]  pb        Beta of the four corners
!> @param[out] s         Local coordinate along P1 to P2, limited to 0 to 1
!> @param[out] t         Local coordinate along P1 to P4, limited to 0 to 1
subroutine inverse_bilinear(target_a, target_b, pa, pb, s, t)

  implicit none

  real(kind=r_def),    intent(in)    :: target_a, target_b
  real(kind=r_def),    intent(in)    :: pa(4), pb(4)
  real(kind=r_def),    intent(inout) :: s, t

  integer(kind=i_def) :: iter
  real(kind=r_def)    :: ea, eb, fa, fb, dsa, dsb, dta, dtb, det, ds, dt

  ea = pa(1) - pa(2) + pa(3) - pa(4)
  eb = pb(1) - pb(2) + pb(3) - pb(4)

  s = 0.5_r_def
  t = 0.5_r_def
  do iter = 1, 30
    fa = pa(1) + s * (pa(2) - pa(1)) + t * (pa(4) - pa(1)) + s * t * ea &
      - target_a
    fb = pb(1) + s * (pb(2) - pb(1)) + t * (pb(4) - pb(1)) + s * t * eb &
      - target_b
    dsa = pa(2) - pa(1) + t * ea
    dsb = pb(2) - pb(1) + t * eb
    dta = pa(4) - pa(1) + s * ea
    dtb = pb(4) - pb(1) + s * eb
    det = dsa * dtb - dta * dsb
    if (ABS(det) < TINY(1.0_r_def)) exit
    ds = (fa * dtb - dta * fb) / det
    dt = (dsa * fb - dsb * fa) / det
    s = s - ds
    t = t - dt
    if (ABS(ds) + ABS(dt) < 1.0e-14_r_def) exit
  end do

  s = MIN(1.0_r_def, MAX(0.0_r_def, s))
  t = MIN(1.0_r_def, MAX(0.0_r_def, t))

end subroutine inverse_bilinear

!> @brief Private function giving the z-component of the cross product of two
!!        2D vectors
!> @param[in] ax  First component of the first vector
!> @param[in] ay  Second component of the first vector
!> @param[in] bx  First component of the second vector
!> @param[in] by  Second component of the second vector
function cross_2d(ax, ay, bx, by) result(c)

  implicit none

  real(kind=r_def), intent(in) :: ax, ay, bx, by
  real(kind=r_def)             :: c

  c = ax * by - ay * bx

end function cross_2d

end module panel_corner_weights_kernel_mod
