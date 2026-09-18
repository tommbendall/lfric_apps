!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
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
module panel_corner_weights_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type, func_type,                          &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_LOGICAL, GH_INTEGER,              &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_9,                    &
                                 ANY_SPACE_9,                                  &
                                 GH_BASIS, GH_EVALUATOR,                       &
                                 CELL_COLUMN, STENCIL, CROSS2D
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
use log_mod, only: log_event, LOG_LEVEL_DEBUG, log_scratch_space

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_corner_weights_kernel_type
  private
  type(arg_type) :: meta_args(11) = (/                                         &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9,                 &
                                                        STENCIL(CROSS2D)),     &
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),  &
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,   &
                                                        STENCIL(CROSS2D)),     &
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                              &
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
!> @param[in]     wx_stencil              DoF map for the Wchi-stencil
!> @param[in]     wx_max_length           Maximum stencil branch length
!> @param[in]     alpha_anti              2D alpha coord extended around corner
!> @param[in]     beta_anti               2D beta coord extended around corner
!> @param[in]     alpha_clock             2D alpha coord extended around corner
!> @param[in]     beta_clock              2D beta coord extended around corner
!> @param[in]     panel_id                ID of cubed sphere panel for each cell
!> @param[in]     pid_stencil_size        Size of panel ID stencil (num cells)
!> @param[in]     pid_stencil             DoF map for the panel id stencil
!> @param[in]     pid_max_length          Maximum stencil branch length
!> @param[in]     panel_edge_dist_W       Distance to the West panel edge
!> @param[in]     panel_edge_dist_E       Distance to the East panel edge
!> @param[in]     panel_edge_dist_S       Distance to the South panel edge
!> @param[in]     panel_edge_dist_N       Distance to the North panel edge
!> @param[in]     cubic_remap             Whether to use cubic remapping. If false
!!                                      linear remapping is used
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
                                     wx_max_length,                            &
                                     wx_stencil_map,                           &
                                     alpha_anti, beta_anti,                    &
                                     alpha_clock, beta_clock,                  &
                                     panel_id,                                 &
                                     pid_stencil_size,                         &
                                     pid_max_length,                           &
                                     pid_stencil_map,                          &
                                     panel_edge_dist_W,                        &
                                     panel_edge_dist_E,                        &
                                     panel_edge_dist_S,                        &
                                     panel_edge_dist_N,                        &
                                     cubic_remap,                              &
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
  integer(kind=i_def), intent(in)    :: wx_stencil_size(4)
  integer(kind=i_def), intent(in)    :: pid_stencil_size(4)
  integer(kind=i_def), intent(in)    :: wx_max_length
  integer(kind=i_def), intent(in)    :: pid_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: remap_depth
  logical(kind=l_def), intent(in)    :: cubic_remap

  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wx_stencil_map(ndf_wx, wx_max_length, 4)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid, pid_max_length, 4)
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
    ! x-type remap around the corner: evaluate the S/N donor stencil in the
    ! anti-clockwise corner coordinates (reprojected into the x-crossing panel)
    call panel_corner_weights_1d(                                            &
            panel_corner_weights_x, panel_corner_indices_x,                    &
            owned_panel, swapped_panel_x,                                      &
            alpha, beta, height, wx_max_length,                                &
            wx_stencil_size(S), wx_stencil_map(:,:,S),                         &
            wx_stencil_size(N), wx_stencil_map(:,:,N),                         &
            alpha_anti, beta_anti,                                             &
            panel_id, pid_max_length,                                          &
            pid_stencil_size(S), pid_stencil_map(:,:,S),                       &
            pid_stencil_size(N), pid_stencil_map(:,:,N),                       &
            cubic_remap,                                                       &
            ndf_ww, undf_ww, map_ww,                                           &
            ndf_wx, undf_wx, map_wx, basis_wx,                                 &
            ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,                     &
            ndf_pid, undf_pid, map_pid                                         &
    )

    ! y-type remap around the corner: evaluate the W/E donor stencil in the
    ! clockwise corner coordinates (reprojected into the y-crossing panel)
    call panel_corner_weights_1d(                                            &
            panel_corner_weights_y, panel_corner_indices_y,                    &
            owned_panel, swapped_panel_y,                                      &
            alpha, beta, height, wx_max_length,                                &
            wx_stencil_size(W), wx_stencil_map(:,:,W),                         &
            wx_stencil_size(E), wx_stencil_map(:,:,E),                         &
            alpha_clock, beta_clock,                                           &
            panel_id, pid_max_length,                                          &
            pid_stencil_size(W), pid_stencil_map(:,:,W),                       &
            pid_stencil_size(E), pid_stencil_map(:,:,E),                       &
            cubic_remap,                                                       &
            ndf_ww, undf_ww, map_ww,                                           &
            ndf_wx, undf_wx, map_wx, basis_wx,                                 &
            ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,                     &
            ndf_pid, undf_pid, map_pid                                         &
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

!> @brief Private routine for calculation of corner-remapping weights and
!!        indices for a single remap direction (x or y).
subroutine panel_corner_weights_1d( remap_weights, remap_indices,                &
                                  owned_panel, swapped_panel,                  &
                                  alpha, beta, height, wx_max_length,          &
                                  wx_stencil_size_l, wx_stencil_l,             &
                                  wx_stencil_size_r, wx_stencil_r,             &
                                  alpha_ext, beta_ext,                         &
                                  panel_id, pid_max_length,                    &
                                  pid_stencil_size_l, pid_stencil_l,           &
                                  pid_stencil_size_r, pid_stencil_r,           &
                                  cubic_remap,                                 &
                                  ndf_ww, undf_ww, map_ww,                     &
                                  ndf_wx, undf_wx, map_wx,                     &
                                  basis_wx,                                    &
                                  ndf_wx_2d, undf_wx_2d, map_wx_2d,            &
                                  basis_wx_2d,                                 &
                                  ndf_pid, undf_pid, map_pid                   &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: wx_stencil_size_l
  integer(kind=i_def), intent(in)    :: wx_stencil_size_r
  integer(kind=i_def), intent(in)    :: wx_max_length
  integer(kind=i_def), intent(in)    :: pid_stencil_size_l
  integer(kind=i_def), intent(in)    :: pid_stencil_size_r
  integer(kind=i_def), intent(in)    :: pid_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: owned_panel
  integer(kind=i_def), intent(in)    :: swapped_panel
  logical(kind=l_def), intent(in)    :: cubic_remap

  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wx_stencil_l(ndf_wx, wx_stencil_size_l)
  integer(kind=i_def), intent(in)    :: wx_stencil_r(ndf_wx, wx_stencil_size_r)
  integer(kind=i_def), intent(in)    :: pid_stencil_l(ndf_pid, pid_stencil_size_l)
  integer(kind=i_def), intent(in)    :: pid_stencil_r(ndf_pid, pid_stencil_size_r)
  real(kind=r_def),    intent(in)    :: basis_wx(1, ndf_wx, ndf_ww)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_ww)
  real(kind=r_def),    intent(in)    :: alpha(undf_wx)
  real(kind=r_def),    intent(in)    :: beta(undf_wx)
  real(kind=r_def),    intent(in)    :: height(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_tran),   intent(inout) :: remap_weights(undf_ww)
  integer(kind=i_def), intent(inout) :: remap_indices(undf_ww)

  ! Internal variables
  integer(kind=i_def), parameter   :: interp_dir_alpha = 1
  integer(kind=i_def), parameter   :: interp_dir_beta  = 2
  real(kind=r_def),    parameter   :: unit_radius = 1.0_r_def
  real(kind=r_def),    allocatable :: x1(:), dx(:)

  integer(kind=i_def) :: panel_edge
  integer(kind=i_def) :: panel
  integer(kind=i_def) :: df, id1, id2, id3, id4, n
  integer(kind=i_def) :: interp_dir
  integer(kind=i_def) :: ncells_in_stencil
  integer(kind=i_def) :: wx_stencil_1d(ndf_wx, 2*wx_max_length-1)
  integer(kind=i_def) :: pid_stencil_1d(2*pid_max_length-1)
  real(kind=r_def)    :: abh(3)
  real(kind=r_def)    :: x0, w1, w2, w3, w4
  real(kind=r_def)    :: xyz(3), h_dummy
  real(kind=r_def)    :: alpha_owned, beta_owned
  real(kind=r_def)    :: alpha_swapped, beta_swapped
  real(kind=r_def)    :: stretch_factor
  real(kind=r_def)    :: inverse_rot_matrix(3,3)
  logical(kind=l_def) :: to_rotate, to_stretch

  ! -------------------------------------------------------------------------- !
  ! Determine direction of interpolation (alpha or beta)
  ! -------------------------------------------------------------------------- !
  panel_edge = 10*owned_panel + swapped_panel
  select case (panel_edge)
  case (15, 26, 34, 43, 51, 62)
    interp_dir = interp_dir_alpha
  case (12, 21, 36, 45, 54, 63)
    interp_dir = interp_dir_beta
  case (16, 41, 32, 25, 64, 53)
    interp_dir = interp_dir_alpha
  case (61, 14, 23, 52, 46, 35)
    interp_dir = interp_dir_beta
  end select

  ! -------------------------------------------------------------------------- !
  ! Create a 1D stencil
  ! -------------------------------------------------------------------------- !
  ncells_in_stencil = wx_stencil_size_l + wx_stencil_size_r - 1

  wx_stencil_1d = 0
  do n = 1, wx_stencil_size_l
    wx_stencil_1d(:,n) = wx_stencil_l(:,n)
  end do
  do n = 1, wx_stencil_size_r - 1
    wx_stencil_1d(:,n+wx_stencil_size_l) = wx_stencil_r(:,n+1)
  end do

  pid_stencil_1d = 0
  do n = 1, pid_stencil_size_l
    pid_stencil_1d(n) = pid_stencil_l(1,n)
  end do
  do n = 1, pid_stencil_size_r - 1
    pid_stencil_1d(n+pid_stencil_size_l) = pid_stencil_r(1,n+1)
  end do

  ! -------------------------------------------------------------------------- !
  ! Compute 1D locations of points used in interpolation
  ! -------------------------------------------------------------------------- !
  abh = 0.0_r_def
  do df = 1, ndf_wx_2d
    abh(1) = abh(1) + alpha_ext(map_wx_2d(df))*basis_wx_2d(1,df,1)
    abh(2) = abh(2) + beta_ext(map_wx_2d(df))*basis_wx_2d(1,df,1)
    abh(3) = abh(3) + height(map_wx(df))*basis_wx_2d(1,df,1)
  end do

  alpha_swapped = abh(1)
  beta_swapped = abh(2)

  ! Convert the interpolation point from the swapped-panel frame back to the
  ! owned-panel frame before comparing against donor stencil points.
  call alphabetar2xyz(abh(1), abh(2), unit_radius, swapped_panel, xyz(1), xyz(2), xyz(3))
  call xyz2alphabetar(xyz(1), xyz(2), xyz(3), owned_panel, abh(1), abh(2), h_dummy)

  if (interp_dir == interp_dir_alpha) then
    x0 = abh(1)
  else
    x0 = abh(2)
  end if

  ! TODO: to remove
  ! write(log_scratch_space,*)                                              &
  !     'CORNER_TARGET_FRAME: owned_panel=', owned_panel,                    &
  !     ' swapped_panel=', swapped_panel,                                    &
  !     ' panel_edge=', panel_edge,                                          &
  !     ' interp_dir=', interp_dir,                                          &
  !     ' alpha_swapped=', alpha_swapped,                                    &
  !     ' beta_swapped=', beta_swapped,                                      &
  !     ' alpha_owned=', abh(1),                                             &
  !     ' beta_owned=', abh(2),                                              &
  !     ' x0=', x0
  ! call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  if (interp_dir == interp_dir_alpha) then
    x0 = abh(1)
  else
    x0 = abh(2)
  end if

  to_rotate = get_to_rotate()
  to_stretch = get_to_stretch()
  if (to_rotate) then
    inverse_rot_matrix = get_inverse_mesh_rotation_matrix()
  end if
  if (to_stretch) then
    stretch_factor = get_stretch_factor()
  end if

  allocate(x1(ncells_in_stencil), dx(ncells_in_stencil))
  do n = 1, ncells_in_stencil
    x1(n) = 0.0_r_def
    dx(n) = 0.0_r_def

    panel = int(panel_id(pid_stencil_1d(n)), i_def)
    if (owned_panel /= panel) then
      x1(n) = LARGE_REAL_POSITIVE
    else
      abh = 0.0_r_def
      do df = 1, ndf_wx
        call chi2xyz(                                                        &
                alpha(wx_stencil_1d(df,n)), beta(wx_stencil_1d(df,n)),       &
                unit_radius, owned_panel, geometry, topology, coord_system,    &
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
        abh(1) = abh(1) + alpha_owned * basis_wx(1,df,1)
        abh(2) = abh(2) + beta_owned * basis_wx(1,df,1)
      end do

      if (interp_dir == interp_dir_alpha) then
        x1(n) = abh(1)
      else
        x1(n) = abh(2)
      end if
    end if

    dx(n) = ABS(x1(n) - x0)

    ! TODO: to remove
    ! write(log_scratch_space,*)                             &
    !     'CORNER_WEIGHTS_1D: n=', n,                                            &
    !     ' owned_panel=', owned_panel,                                          &
    !     ' donor_panel_id=', int(panel_id(pid_stencil_1d(n)), i_def),           &
    !     ' x1=', x1(n),                                                        &
    !     ' x0=', x0
    ! call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
  end do

  id1 = MINLOC(dx(1:ncells_in_stencil), 1)
  dx(id1) = LARGE_REAL_POSITIVE
  id2 = MINLOC(dx(1:ncells_in_stencil), 1)
  dx(id2) = LARGE_REAL_POSITIVE
  id3 = MINLOC(dx(1:ncells_in_stencil), 1)
  dx(id3) = LARGE_REAL_POSITIVE
  id4 = MINLOC(dx(1:ncells_in_stencil), 1)

  ! TODO: to remove
  ! write(log_scratch_space, *) 'CORNER_WEIGHTS_1D: ids=', id1, id2, id3, id4, &
  !     ' x1(ids)=', x1(id1), x1(id2), x1(id3), x1(id4)
  ! call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  if (ABS(x1(id1) - x1(id2)) <= epsilon(1.0_r_def) .or.                          &
      ABS(x1(id1) - x1(id3)) <= epsilon(1.0_r_def) .or.                          &
      ABS(x1(id1) - x1(id4)) <= epsilon(1.0_r_def) .or.                          &
      ABS(x1(id2) - x1(id3)) <= epsilon(1.0_r_def) .or.                          &
      ABS(x1(id2) - x1(id4)) <= epsilon(1.0_r_def) .or.                          &
      ABS(x1(id3) - x1(id4)) <= epsilon(1.0_r_def)) then
    remap_weights(map_ww(1):map_ww(1)+3) = 0.0_r_tran
    remap_weights(map_ww(1)) = 1.0_r_tran
    remap_indices(map_ww(1):map_ww(1)+3) = id1
  else
    ! The remapping is cubic if requested; otherwise it falls back to linear
    if (cubic_remap) then
      w1 = (x0 - x1(id2))/(x1(id1) - x1(id2)) * (x0 - x1(id3))/(x1(id1) - x1(id3)) * &
           (x0 - x1(id4))/(x1(id1) - x1(id4))
      w2 = (x0 - x1(id1))/(x1(id2) - x1(id1)) * (x0 - x1(id3))/(x1(id2) - x1(id3)) * &
           (x0 - x1(id4))/(x1(id2) - x1(id4))
      w3 = (x0 - x1(id1))/(x1(id3) - x1(id1)) * (x0 - x1(id2))/(x1(id3) - x1(id2)) * &
           (x0 - x1(id4))/(x1(id3) - x1(id4))
      w4 = (x0 - x1(id1))/(x1(id4) - x1(id1)) * (x0 - x1(id2))/(x1(id4) - x1(id2)) * &
           (x0 - x1(id3))/(x1(id4) - x1(id3))
      ! Fill the output arrays with the four weights and their associated indices.
      remap_weights(map_ww(1)) = w1
      remap_weights(map_ww(1)+1) = w2
      remap_weights(map_ww(1)+2) = w3
      remap_weights(map_ww(1)+3) = w4
      remap_indices(map_ww(1)) = id1
      remap_indices(map_ww(1)+1) = id2
      remap_indices(map_ww(1)+2) = id3
      remap_indices(map_ww(1)+3) = id4
    else
      w1 = (x0 - x1(id2))/(x1(id1) - x1(id2))
      w2 = (x0 - x1(id1))/(x1(id2) - x1(id1))
      remap_weights(map_ww(1)) = w1
      remap_weights(map_ww(1)+1) = w2
      remap_indices(map_ww(1)) = id1
      remap_indices(map_ww(1)+1) = id2
    end if
  end if

  deallocate(x1, dx)

end subroutine panel_corner_weights_1d

end module panel_corner_weights_kernel_mod
