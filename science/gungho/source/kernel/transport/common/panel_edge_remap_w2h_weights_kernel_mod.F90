!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Compute the interpolation weights/indices and rotation factors
!!        needed to remap a W2H flux/wind field from the standard
!!        cubed-sphere mesh onto the "extended"/eave mesh at the edges of
!!        the mesh panels.
!!
!> @details This kernel generalises the donor-stencil, Lagrange-interpolation
!!          approach of panel_edge_weights_kernel_mod (used for W3/Wtheta
!!          scalar fields) to a W2H field. Because a W2H field carries four
!!          face values per cell (W, S, E, N), remapping it across a panel
!!          edge additionally requires accounting for the change of basis
!!          vectors between the native and eave coordinate frames: a face
!!          value that was "parallel" to the native frame picks up a
!!          "perpendicular" contribution once expressed in the eave frame.
!!
!!          Each selected face is computed for each adjacent panel. The panel
!!          axis label selects the storage slot, and donor indices store signed
!!          offsets from the face rather than positions in one cell's stencil.
!!
!!          Interpolation and rotation are kept separate, unlike the previous
!!          version of this kernel:
!!            - interp_w2h_weights / interp_w2h_indices are PURE Lagrange
!!              donor-stencil weights for interpolating a "parl_w2" value,
!!            - parl_w2_rot / perp_w2_rot are the cos/sin
!!              rotation factors between the native and eave coordinate
!!              frames, applied separately (by panel_edge_remap_w2h_kernel_mod)
!!              to the "parl_w2" and "perp_w2" staged values respectively.
!!          This phase only produces a non-identity "parl_w2" interpolation
!!          (the W2H donors share the target's own face orientation); the
!!          "perp_w2" contribution (from the W1V/W3-staged values computed by
!!          later phases) is represented here only by its rotation factor,
!!          perp_w2_rot, which is correct but has no effect until a later
!!          phase supplies a non-zero perp_w2 staged value.
!!
!!          interp_perp_weights / interp_perp_indices are the
!!          weights/indices for the "perp_w2" interpolation. Each parl_p face
!!          slot (S/N for x, W/E for y) holds two sets of ndata entries, one
!!          per face vertex. For x the donors are the W DoFs (west vertex) and
!!          E DoFs (east vertex); for y the S DoFs (south vertex) and N DoFs
!!          (north vertex). The target points are the eave vertex positions,
!!          evaluated via the gh_evaluator_targets(2) ("wv") basis arrays.
!!
!!          interp_w3_weights_{x,y} / interp_w3_indices_{x,y} are the
!!          weights/indices for the "parl_p donors -> W3 centre targets"
!!          staging interpolation (Phase 3 work, not yet done). In this
!!          phase they are written as identity/zero placeholders only, so
!!          that the fields exist with stable meta_args ahead of a later
!!          phase filling them in.
!!
!!          Away from any relevant panel edge, interp_w2h_weights defaults to
!!          the identity (donor = target's own DoF, weight = 1), parl_w2_rot
!!          defaults to 1 and perp_w2_rot defaults to 0, so that applying the
!!          remap leaves the field unchanged there.
!!
!!          NOTE: the exact combination of the Lagrange donor weight and the
!!          rotation factor implemented below is a new derivation and has
!!          not yet been validated. The intended acceptance test is that
!!          remapping an analytic zonal wind u = U0*cos(lat) across a panel
!!          edge should reproduce the expected value accurately.
module panel_edge_remap_w2h_weights_kernel_mod

use kernel_mod,             only: kernel_type
use argument_mod,           only: arg_type, func_type,                         &
                                  GH_FIELD, GH_SCALAR,                         &
                                  GH_REAL, GH_LOGICAL, GH_INTEGER,             &
                                  GH_READ, GH_WRITE,                           &
                                  GH_BASIS, GH_EVALUATOR,                      &
                                  ANY_DISCONTINUOUS_SPACE_3,                   &
                                  ANY_DISCONTINUOUS_SPACE_5,                   &
                                  ANY_DISCONTINUOUS_SPACE_7,                   &
                                  ANY_DISCONTINUOUS_SPACE_9,                   &
                                  ANY_DISCONTINUOUS_SPACE_10,                  &
                                  ANY_SPACE_9,                                 &
                                  CELL_COLUMN, STENCIL, CROSS2D
use fs_continuity_mod,       only: W2h, W2, W3
use constants_mod,           only: r_tran, r_def, i_def, l_def,               &
                                   LARGE_REAL_POSITIVE
use reference_element_mod,   only: W, S, N, E, SW, SE, NE, NW
use panel_edge_support_mod,  only: panel_neighbour
use sci_face_selector_support_mod, only: face_from_face_selector
use coord_transform_mod,     only: PANEL_ROT_MATRIX
use panel_edge_support_mod,   only: panel_axis_label

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!! PSy layer
type, public, extends(kernel_type) :: panel_edge_remap_w2h_weights_kernel_type
  private
    type(arg_type) :: meta_args(30) = (/                                         &
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! interp_w2h_weights
      arg_type(GH_FIELD*3, GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! interp_w2h_indices
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, W2h),                       & ! parl_w2_rot
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, W2h),                       & ! perp_w2_rot
      ! The perp weights/indices hold, per W2H face slot, two consecutive
      ! sets of ndata entries for the face's two vertices
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_10),& ! interp_perp_weights
      arg_type(GH_FIELD*3, GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_10),& ! interp_perp_indices
       ! NOTE: W1V is not yet a valid function space for GH_FIELD kernel
       ! metadata in PSyclone's LFRic API, so ANY_DISCONTINUOUS_SPACE_7 is
       ! used for this W1V field. It is never read: its only purpose is to
       ! provide the W1V nodes at which the vertex basis functions are
       ! evaluated
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7), & ! w1v_nodes (dummy)
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                        & ! interp_w3_weights_x (placeholder)
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, W3),                        & ! interp_w3_indices_x (placeholder)
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                        & ! interp_w3_weights_y (placeholder)
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, W3),                        & ! interp_w3_indices_y (placeholder)
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! eave_lon
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! eave_lat
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! correct_flux
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9,                &
                                                        STENCIL(CROSS2D)),     & ! chi
      arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9, &
                    STENCIL(CROSS2D)),     & ! alpha_x, beta_x
      arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9, &
                    STENCIL(CROSS2D)),     & ! alpha_y, beta_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2),                        & ! dA_at_w2 (diagnostic)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,  &
                                                        STENCIL(CROSS2D)),     & ! panel_id
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3, &
                    STENCIL(CROSS2D)),     & ! panel_edge_dist_W/E/S/N
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ew
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ns
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                             & ! cubic_remap
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! remap_depth
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! geometry
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! topology
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! coord_system
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),                             & ! scaled_radius
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),                             & ! u0 (diagnostic)
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              & ! dt (diagnostic)
  /)
  type(func_type) :: meta_funcs(2) = (/                                        &
      func_type(ANY_SPACE_9, GH_BASIS),                                        &
      func_type(ANY_DISCONTINUOUS_SPACE_9, GH_BASIS)                           &
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
  ! Three evaluator targets: ANY_DISCONTINUOUS_SPACE_5 ("ww", W2H-shaped face
  ! storage, set 1), ANY_DISCONTINUOUS_SPACE_7 ("wv", the W1V nodes of the
  ! dummy W1V field, at which the face-vertex positions are evaluated) and W3
  ! ("w3", W3 cell-centre storage, set 3). PSyclone generates one basis array
  ! per (meta_funcs, gh_evaluator_targets) pair - see
  ! sci_w3_to_w2_displacement_kernel_mod.F90
  integer :: gh_evaluator_targets(3) = (/ ANY_DISCONTINUOUS_SPACE_5,          &
                                           ANY_DISCONTINUOUS_SPACE_7,         &
                                           W3 /)
contains
  procedure, nopass :: panel_edge_remap_w2h_weights_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_remap_w2h_weights_code

contains

!> @brief Compute the interpolation weights/indices and rotation factors to
!!        remap a W2H field across cubed-sphere panel edges.
!> @param[in]     nlayers               Number of layers
!> @param[in,out] interp_w2h_weights_1 Weights for panel-axis label 1
!> @param[in,out] interp_w2h_weights_2 Weights for panel-axis label 2
!> @param[in,out] interp_w2h_weights_3 Weights for panel-axis label 3
!> @param[in,out] interp_w2h_indices_1 Signed donor offsets for label 1
!> @param[in,out] interp_w2h_indices_2 Signed donor offsets for label 2
!> @param[in,out] interp_w2h_indices_3 Signed donor offsets for label 3
!> @param[in,out] parl_w2_rot_1 Parallel rotation factor for label 1
!> @param[in,out] parl_w2_rot_2 Parallel rotation factor for label 2
!> @param[in,out] parl_w2_rot_3 Parallel rotation factor for label 3
!> @param[in,out] perp_w2_rot_1 Perpendicular rotation factor for label 1
!> @param[in,out] perp_w2_rot_2 Perpendicular rotation factor for label 2
!> @param[in,out] perp_w2_rot_3 Perpendicular rotation factor for label 3
!> @param[in,out] interp_perp_weights_1 Perpendicular weights for label 1
!> @param[in,out] interp_perp_weights_2 Perpendicular weights for label 2
!> @param[in,out] interp_perp_weights_3 Perpendicular weights for label 3
!> @param[in,out] interp_perp_indices_1 Signed vertex donor offsets for label 1
!> @param[in,out] interp_perp_indices_2 Signed vertex donor offsets for label 2
!> @param[in,out] interp_perp_indices_3 Signed vertex donor offsets for label 3
!> @param[in]     w1v_nodes             Dummy W1V field, only used to evaluate
!!                                      basis functions at the W1V nodes
!> @param[in,out] interp_w3_weights_x   Placeholder (identity/zero) for the
!!                                      "parl_p -> W3" staging weights (x)
!> @param[in,out] interp_w3_indices_x   Placeholder donor indices (x)
!> @param[in,out] interp_w3_weights_y   Placeholder (identity/zero) for the
!!                                      "parl_p -> W3" staging weights (y)
!> @param[in,out] interp_w3_indices_y   Placeholder donor indices (y)
!> @param[in,out] eave_lon_1            Diagnostic eave longitude for label 1
!> @param[in,out] eave_lon_2            Diagnostic eave longitude for label 2
!> @param[in,out] eave_lon_3            Diagnostic eave longitude for label 3
!> @param[in,out] eave_lat_1            Diagnostic eave latitude for label 1
!> @param[in,out] eave_lat_2            Diagnostic eave latitude for label 2
!> @param[in,out] eave_lat_3            Diagnostic eave latitude for label 3
!> @param[in,out] correct_flux_1        Diagnostic analytic flux for label 1
!> @param[in,out] correct_flux_2        Diagnostic analytic flux for label 2
!> @param[in,out] correct_flux_3        Diagnostic analytic flux for label 3
!> @param[in]     chi_1                 Native alpha coordinate field
!> @param[in]     chi_2                 Native beta coordinate field
!> @param[in]     chi_3                 Height coordinate field
!> @param[in]     wx_stencil_size       Size of the chi CROSS2D stencil branches
!> @param[in]     wx_max_length         Maximum stencil branch length
!> @param[in]     wx_stencil_map        DoF map for the chi stencil
!> @param[in]     alpha_x               2D alpha coord for mesh extended in x
!> @param[in]     beta_x                2D beta coord for mesh extended in x
!> @param[in]     alpha_y               2D alpha coord for mesh extended in y
!> @param[in]     beta_y                2D beta coord for mesh extended in y
!> @param[in]     dA_at_w2              W2 face-area field (diagnostic)
!> @param[in]     panel_id              Field containing IDs of mesh panels
!> @param[in]     pid_stencil_size      Size of panel ID CROSS2D stencil branches
!> @param[in]     pid_max_length        Maximum stencil branch length
!> @param[in]     pid_stencil_map       DoF map for the panel ID stencil
!> @param[in]     panel_edge_dist_W     Distance of this column to the W panel edge
!> @param[in]     panel_edge_dist_E     Distance of this column to the E panel edge
!> @param[in]     panel_edge_dist_S     Distance of this column to the S panel edge
!> @param[in]     panel_edge_dist_N     Distance of this column to the N panel edge
!> @param[in]     face_selector_ew      East-West face selector, used to ensure
!!                                      each W/E face is only ever written by
!!                                      one of its two neighbouring columns
!> @param[in]     face_selector_ns      North-South face selector, as
!!                                      face_selector_ew for S/N faces
!> @param[in]     cubic_remap           Whether to use cubic (else linear) remapping
!> @param[in]     remap_depth           Depth (in columns) to which remapping applies
!> @param[in]     geometry              Mesh geometry identifier
!> @param[in]     topology              Mesh topology identifier
!> @param[in]     coord_system          Coordinate system identifier
!> @param[in]     scaled_radius         Planet radius, appropriately scaled
!> @param[in]     u0                    Solid-body-rotation zonal wind
!!                                      amplitude (diagnostic)
!> @param[in]     dt                    Model timestep (diagnostic)
!> @param[in]     ndf_ww                Num DoFs per cell for the W2H-shaped
!!                                      interpolation/diagnostic storage
!> @param[in]     undf_ww               Num DoFs for this partition for the
!!                                      W2H-shaped interpolation/diagnostic
!!                                      storage
!> @param[in]     map_ww                DoF map for the W2H-shaped
!!                                      interpolation/diagnostic storage
!> @param[in]     ndf_wr                Num DoFs per cell for the rotation factors
!> @param[in]     undf_wr               Num DoFs for this partition for the
!!                                      rotation factors
!> @param[in]     map_wr                DoF map for the rotation factors
!> @param[in]     ndf_wp                Num DoFs per cell for the perp
!!                                      interpolation weights/indices
!> @param[in]     undf_wp               Num DoFs for this partition for the
!!                                      perp interpolation weights/indices
!> @param[in]     map_wp                DoF map for the perp interpolation
!!                                      weights/indices
!> @param[in]     ndf_wv                Num DoFs per cell for the dummy W1V field
!> @param[in]     undf_wv               Num DoFs for this partition for the
!!                                      dummy W1V field
!> @param[in]     map_wv                DoF map for the dummy W1V field
!> @param[in]     ndf_w3                Num DoFs per cell for the W3 staging
!!                                      placeholders
!> @param[in]     undf_w3               Num DoFs for this partition for the
!!                                      W3 staging placeholders
!> @param[in]     map_w3                DoF map for the W3 staging placeholders
!> @param[in]     ndf_wx                Num DoFs per cell for coordinate fields
!> @param[in]     undf_wx               Num DoFs for this partition for
!!                                      coordinate fields
!> @param[in]     map_wx                DoF map for coordinate fields
!> @param[in]     basis_wx_ww           Basis functions of chi, evaluated at
!!                                      W2H ("ww") nodes
!> @param[in]     basis_wx_wv           Basis functions of chi, evaluated at
!!                                      W1V ("wv") vertex nodes
!> @param[in]     basis_wx_w3           Basis functions of chi, evaluated at
!!                                      W3 ("w3") centre nodes
!> @param[in]     ndf_wx_2d             Num DoFs per cell for 2D coord fields
!> @param[in]     undf_wx_2d            Num 2D coord DoFs for this partition
!> @param[in]     map_wx_2d             DoF map for 2D coord fields
!> @param[in]     basis_wx_2d_ww        Basis functions of extended coords, at
!!                                      W2H ("ww") nodes
!> @param[in]     basis_wx_2d_wv        Basis functions of extended coords, at
!!                                      W1V ("wv") vertex nodes
!> @param[in]     basis_wx_2d_w3        Basis functions of extended coords, at
!!                                      W3 ("w3") centre nodes
!> @param[in]     ndf_w2                Num DoFs per cell for dA_at_w2
!> @param[in]     undf_w2               Num DoFs for this partition for dA_at_w2
!> @param[in]     map_w2                DoF map for dA_at_w2
!> @param[in]     ndf_pid               Num DoFs per cell for panel ID field
!> @param[in]     undf_pid              Num DoFs for this partition for panel ID
!> @param[in]     map_pid               DoF map for panel ID field
subroutine panel_edge_remap_w2h_weights_code(                                  &
    nlayers,                                                                   &
  interp_w2h_weights_1, interp_w2h_weights_2, interp_w2h_weights_3,          &
  interp_w2h_indices_1, interp_w2h_indices_2, interp_w2h_indices_3,          &
  parl_w2_rot_1, parl_w2_rot_2, parl_w2_rot_3,                              &
  perp_w2_rot_1, perp_w2_rot_2, perp_w2_rot_3,                              &
  interp_perp_weights_1, interp_perp_weights_2, interp_perp_weights_3,        &
  interp_perp_indices_1, interp_perp_indices_2, interp_perp_indices_3,        &
    w1v_nodes,                                                                 &
    interp_w3_weights_x, interp_w3_indices_x,                                  &
    interp_w3_weights_y, interp_w3_indices_y,                                  &
    eave_lon_1, eave_lon_2, eave_lon_3,                                        &
    eave_lat_1, eave_lat_2, eave_lat_3,                                        &
    correct_flux_1, correct_flux_2, correct_flux_3,                            &
    chi_1, chi_2, chi_3,                                                       &
    wx_stencil_size, wx_max_length, wx_stencil_map,                            &
    alpha_x, beta_x, coord_x_stencil_size, coord_x_max_length,                 &
    coord_x_stencil_map,                                                       &
    alpha_y, beta_y, coord_y_stencil_size, coord_y_max_length,                 &
    coord_y_stencil_map,                                                       &
    dA_at_w2,                                                                  &
    panel_id,                                                                  &
    pid_stencil_size, pid_max_length, pid_stencil_map,                         &
    panel_edge_dist_W, panel_edge_dist_E,                                      &
    panel_edge_dist_S, panel_edge_dist_N,                                      &
    edge_stencil_size, edge_max_length, edge_stencil_map,                      &
    face_selector_ew, face_selector_ns,                                        &
    cubic_remap, remap_depth,                                                  &
    geometry, topology, coord_system, scaled_radius,                           &
    u0, dt,                                                                    &
    ndf_ww, undf_ww, map_ww,                                                   &
    ndf_wr, undf_wr, map_wr,                                                   &
    ndf_wp, undf_wp, map_wp,                                                   &
    ndf_wv, undf_wv, map_wv,                                                   &
    ndf_w3, undf_w3, map_w3,                                                   &
    ndf_wx, undf_wx, map_wx, basis_wx_ww, basis_wx_wv, basis_wx_w3,           &
    ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d_ww, basis_wx_2d_wv,          &
    basis_wx_2d_w3,                                                           &
    ndf_w2, undf_w2, map_w2,                                                   &
    ndf_pid, undf_pid, map_pid                                                 &
)

  use sci_chi_transform_mod, only: chi2xyz, get_to_rotate, get_to_stretch,     &
                                   get_inverse_mesh_rotation_matrix,           &
                                   get_stretch_factor
  use coord_transform_mod,   only: alphabetar2xyz, xyz2alphabetar, xyz2ll,     &
                                   inverse_schmidt_transform_xyz,             &
                                   alphabetar2xyz_vector, cart2sphere_vector

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: wx_stencil_size(4)
  integer(kind=i_def), intent(in)    :: pid_stencil_size(4)
  integer(kind=i_def), intent(in)    :: wx_max_length
  integer(kind=i_def), intent(in)    :: pid_max_length
  integer(kind=i_def), intent(in)    :: coord_x_stencil_size(4)
  integer(kind=i_def), intent(in)    :: coord_x_max_length
  integer(kind=i_def), intent(in)    :: coord_y_stencil_size(4)
  integer(kind=i_def), intent(in)    :: coord_y_max_length
  integer(kind=i_def), intent(in)    :: edge_stencil_size(4)
  integer(kind=i_def), intent(in)    :: edge_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_wp, undf_wp
  integer(kind=i_def), intent(in)    :: ndf_wv, undf_wv
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: remap_depth
  logical(kind=l_def), intent(in)    :: cubic_remap
  integer(kind=i_def), intent(in)    :: geometry
  integer(kind=i_def), intent(in)    :: topology
  integer(kind=i_def), intent(in)    :: coord_system
  real(kind=r_def),    intent(in)    :: scaled_radius
  real(kind=r_def),    intent(in)    :: u0
  real(kind=r_def),    intent(in)    :: dt

  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_wp(ndf_wp)
  integer(kind=i_def), intent(in)    :: map_wv(ndf_wv)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wx_stencil_map(ndf_wx, wx_max_length, 4)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid, pid_max_length, 4)
    integer(kind=i_def), intent(in)    :: coord_x_stencil_map(                  &
      ndf_wx_2d, coord_x_max_length, 4)
    integer(kind=i_def), intent(in)    :: coord_y_stencil_map(                  &
      ndf_wx_2d, coord_y_max_length, 4)
    integer(kind=i_def), intent(in)    :: edge_stencil_map(                     &
      ndf_pid, edge_max_length, 4)
  real(kind=r_def),    intent(in)    :: basis_wx_ww(1, ndf_wx, ndf_ww)
  ! basis_wx_wv is generated because gh_evaluator_targets(2) applies to the
  ! chi func_type too, but is not needed/used: donor node positions are
  ! always sampled at the source's own "ww" face node (see
  ! compute_lagrange_weights), never at a "wv" vertex node
  real(kind=r_def),    intent(in)    :: basis_wx_wv(1, ndf_wx, ndf_wv)
  ! basis_wx_w3 is generated for the same reason as basis_wx_wv, and is
  ! likewise unused: donor nodes are always sampled at the source's own
  ! "ww" face node
  real(kind=r_def),    intent(in)    :: basis_wx_w3(1, ndf_wx, ndf_w3)
  real(kind=r_def),    intent(in)    :: basis_wx_2d_ww(1, ndf_wx_2d, ndf_ww)
  real(kind=r_def),    intent(in)    :: basis_wx_2d_wv(1, ndf_wx_2d, ndf_wv)
  real(kind=r_def),    intent(in)    :: basis_wx_2d_w3(1, ndf_wx_2d, ndf_w3)

  real(kind=r_tran),   intent(inout) :: interp_w2h_weights_1(undf_ww)
  real(kind=r_tran),   intent(inout) :: interp_w2h_weights_2(undf_ww)
  real(kind=r_tran),   intent(inout) :: interp_w2h_weights_3(undf_ww)
  integer(kind=i_def), intent(inout) :: interp_w2h_indices_1(undf_ww)
  integer(kind=i_def), intent(inout) :: interp_w2h_indices_2(undf_ww)
  integer(kind=i_def), intent(inout) :: interp_w2h_indices_3(undf_ww)
  real(kind=r_tran),   intent(inout) :: parl_w2_rot_1(undf_wr)
  real(kind=r_tran),   intent(inout) :: parl_w2_rot_2(undf_wr)
  real(kind=r_tran),   intent(inout) :: parl_w2_rot_3(undf_wr)
  real(kind=r_tran),   intent(inout) :: perp_w2_rot_1(undf_wr)
  real(kind=r_tran),   intent(inout) :: perp_w2_rot_2(undf_wr)
  real(kind=r_tran),   intent(inout) :: perp_w2_rot_3(undf_wr)
  real(kind=r_tran),   intent(inout) :: interp_perp_weights_1(undf_wp)
  real(kind=r_tran),   intent(inout) :: interp_perp_weights_2(undf_wp)
  real(kind=r_tran),   intent(inout) :: interp_perp_weights_3(undf_wp)
  integer(kind=i_def), intent(inout) :: interp_perp_indices_1(undf_wp)
  integer(kind=i_def), intent(inout) :: interp_perp_indices_2(undf_wp)
  integer(kind=i_def), intent(inout) :: interp_perp_indices_3(undf_wp)
  real(kind=r_tran),   intent(in)    :: w1v_nodes(undf_wv)
  real(kind=r_tran),   intent(inout) :: interp_w3_weights_x(undf_w3)
  integer(kind=i_def), intent(inout) :: interp_w3_indices_x(undf_w3)
  real(kind=r_tran),   intent(inout) :: interp_w3_weights_y(undf_w3)
  integer(kind=i_def), intent(inout) :: interp_w3_indices_y(undf_w3)
  real(kind=r_tran),   intent(inout) :: eave_lon_1(undf_ww)
  real(kind=r_tran),   intent(inout) :: eave_lon_2(undf_ww)
  real(kind=r_tran),   intent(inout) :: eave_lon_3(undf_ww)
  real(kind=r_tran),   intent(inout) :: eave_lat_1(undf_ww)
  real(kind=r_tran),   intent(inout) :: eave_lat_2(undf_ww)
  real(kind=r_tran),   intent(inout) :: eave_lat_3(undf_ww)
  real(kind=r_tran),   intent(inout) :: correct_flux_1(undf_ww)
  real(kind=r_tran),   intent(inout) :: correct_flux_2(undf_ww)
  real(kind=r_tran),   intent(inout) :: correct_flux_3(undf_ww)

  real(kind=r_def),    intent(in)    :: chi_1(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_2(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_3(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: alpha_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: dA_at_w2(undf_w2)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)

  ! Internal variables
  integer(kind=i_def), parameter   :: interp_dir_alpha = 1
  integer(kind=i_def), parameter   :: interp_dir_beta  = 2
  real(kind=r_def),    parameter   :: unit_radius = 1.0_r_def

  integer(kind=i_def) :: df, dfc, i, w2h_idx, wr_idx, wp_idx, w3_idx, pid_idx
  integer(kind=i_def) :: ew_sel, ns_sel, j
  integer(kind=i_def) :: ndata
  integer(kind=i_def) :: edge_dist(4)
  integer(kind=i_def) :: owned_panel, neighbour_panel
  integer(kind=i_def) :: other_panel, other_df, host_sign, label
  integer(kind=i_def) :: edge_direction_x, edge_direction_y
  real(kind=r_def) :: label_weights(4), label_perp_weights(8)
  real(kind=r_def) :: label_rotation_parallel, label_rotation_perp
  real(kind=r_def) :: label_eave_lon, label_eave_lat, label_correct_flux
  integer(kind=i_def) :: label_offsets(4), label_perp_offsets(8)
  logical(kind=l_def) :: near_panel_edge_x, near_panel_edge_y
  real(kind=r_def)    :: stretch_factor
  real(kind=r_def)    :: inverse_rot_matrix(3,3)
  logical(kind=l_def) :: to_rotate, to_stretch

  ! Lagrange donor-stencil interpolation results, for the currently active
  ! direction (x or y)
  real(kind=r_def)    :: lagrange_w(4)
  integer(kind=i_def) :: lagrange_id(4)

  ! Perp interpolation results for the two vertices of a face
  real(kind=r_def)    :: perp_w_a(4), perp_w_b(4)
  integer(kind=i_def) :: perp_id_a(4), perp_id_b(4)

  to_rotate = get_to_rotate()
  to_stretch = get_to_stretch()
  if (to_rotate) then
    inverse_rot_matrix = get_inverse_mesh_rotation_matrix()
  end if
  if (to_stretch) then
    stretch_factor = get_stretch_factor()
  end if

  ndata = merge(4, 2, cubic_remap)

  pid_idx = map_pid(1)
  owned_panel = int(panel_id(pid_idx), i_def)

  ew_sel = face_selector_ew(pid_idx)
  ns_sel = face_selector_ns(pid_idx)
  edge_dist = (/ panel_edge_dist_W(pid_idx), panel_edge_dist_E(pid_idx),      &
                 panel_edge_dist_S(pid_idx), panel_edge_dist_N(pid_idx) /)

  ! ======================================================================== !
  ! Set 3: "parl_p donors -> W3 centre targets" interpolation weights/
  ! indices (Route B of the instructions doc, Sec 2.4). There is a single W3
  ! DoF per cell column (unlike W2h/W1V, which have one per horizontal
  ! face), so this is done once, not per face. Computed to remap_depth + 1
  ! (Sec 2.5): the averaging step done later by the average kernel needs
  ! the W3 centre of the cell across the target face, which may be one
  ! column further out than remap_depth. The donor is always the native S
  ! DoF for the x remap and the native W DoF for the y remap (different
  ! columns use different donor-stencil indices, not a different source DoF)
  w3_idx = map_w3(1)

  near_panel_edge_x = (                                                       &
      ABS(panel_edge_dist_W(pid_idx)) <= remap_depth + 1                     &
      .or. ABS(panel_edge_dist_E(pid_idx)) <= remap_depth + 1                &
  )
  if (ABS(panel_edge_dist_W(pid_idx)) < ABS(panel_edge_dist_E(pid_idx))) then
    edge_direction_x = W
  else
    edge_direction_x = E
  end if

  if (near_panel_edge_x) then
    neighbour_panel = panel_neighbour(owned_panel, edge_direction_x)
    call compute_lagrange_weights(1, S, .true., neighbour_panel,              &
                                   lagrange_w, lagrange_id, is_w3_target=.true.)
    do i = 1, ndata
      interp_w3_weights_x(w3_idx+i-1) = REAL(lagrange_w(i), r_tran)
      interp_w3_indices_x(w3_idx+i-1) = lagrange_id(i)
    end do
  else
    interp_w3_weights_x(w3_idx) = 1.0_r_tran
    interp_w3_indices_x(w3_idx) = 1
    do i = 2, ndata
      interp_w3_weights_x(w3_idx+i-1) = 0.0_r_tran
      interp_w3_indices_x(w3_idx+i-1) = 1
    end do
  end if

  near_panel_edge_y = (                                                       &
      ABS(panel_edge_dist_S(pid_idx)) <= remap_depth + 1                     &
      .or. ABS(panel_edge_dist_N(pid_idx)) <= remap_depth + 1                &
  )
  if (ABS(panel_edge_dist_S(pid_idx)) < ABS(panel_edge_dist_N(pid_idx))) then
    edge_direction_y = S
  else
    edge_direction_y = N
  end if

  if (near_panel_edge_y) then
    neighbour_panel = panel_neighbour(owned_panel, edge_direction_y)
    call compute_lagrange_weights(1, W, .false., neighbour_panel,             &
                                   lagrange_w, lagrange_id, is_w3_target=.true.)
    do i = 1, ndata
      interp_w3_weights_y(w3_idx+i-1) = REAL(lagrange_w(i), r_tran)
      interp_w3_indices_y(w3_idx+i-1) = lagrange_id(i)
    end do
  else
    interp_w3_weights_y(w3_idx) = 1.0_r_tran
    interp_w3_indices_y(w3_idx) = 1
    do i = 2, ndata
      interp_w3_weights_y(w3_idx+i-1) = 0.0_r_tran
      interp_w3_indices_y(w3_idx+i-1) = 1
    end do
  end if

  ! Compute each selected face from both adjacent panels so either selector
  ! owner produces the same labelled values.
  do j = 1, ABS(ew_sel) + ABS(ns_sel)
    df = face_from_face_selector(j, ew_sel, ns_sel)

    w2h_idx = map_ww(df)
    wr_idx = map_wr(df)
    wp_idx = map_wp(df)

    other_panel = owned_panel
    other_df = opposite_face(df)
    if (pid_stencil_size(df) >= 2) then
      other_panel = int(panel_id(pid_stencil_map(1, 2, df)), i_def)
      if (other_panel /= owned_panel) then
        other_df = face_direction_to_panel(other_panel, owned_panel)
      end if
    end if

    if (other_panel /= owned_panel) then
      if (owned_panel < other_panel) then
        host_sign = 1_i_def
      else
        host_sign = -1_i_def
      end if
    else
      if (df == S .or. df == W) then
        host_sign = 1_i_def
      else
        host_sign = -1_i_def
      end if
    end if

    do label = 1, 3
      label_weights = 0.0_r_def
      label_offsets = 0_i_def
      label_perp_weights = 0.0_r_def
      label_perp_offsets = 0_i_def
      label_rotation_parallel = 1.0_r_def
      label_rotation_perp = 0.0_r_def
      label_eave_lon = 0.0_r_def
      label_eave_lat = 0.0_r_def
      label_correct_flux = 0.0_r_def

      if (label == panel_axis_label(owned_panel)) then
        call compute_labelled_side(                                           &
            df, owned_panel, .false., host_sign, label_weights, label_offsets,&
            label_perp_weights, label_perp_offsets,                          &
            label_rotation_parallel, label_rotation_perp,                    &
            label_eave_lon, label_eave_lat, label_correct_flux               &
        )
      else if (other_panel /= owned_panel .and.                               &
               label == panel_axis_label(other_panel)) then
        call compute_labelled_side(                                           &
            other_df, other_panel, .true., host_sign, label_weights,          &
            label_offsets, label_perp_weights, label_perp_offsets,            &
            label_rotation_parallel, label_rotation_perp,                    &
            label_eave_lon, label_eave_lat, label_correct_flux               &
        )
      end if

      call store_labelled_fields(label, w2h_idx, wr_idx, wp_idx,              &
          label_weights, label_offsets, label_perp_weights,                   &
          label_perp_offsets, label_rotation_parallel, label_rotation_perp,   &
          label_eave_lon, label_eave_lat, label_correct_flux)
    end do

  end do  ! faces

contains

  subroutine compute_labelled_side(target_df, target_panel, target_across,    &
                                   host_sign, weights, offsets,               &
                                   perp_weights, perp_offsets,                &
                                   parallel_factor, perp_factor,              &
                                   eave_lon, eave_lat, correct_flux)

    implicit none

    integer(kind=i_def), intent(in) :: target_df, target_panel, host_sign
    logical(kind=l_def), intent(in) :: target_across
    real(kind=r_def), intent(out) :: weights(4), perp_weights(8)
    integer(kind=i_def), intent(out) :: offsets(4), perp_offsets(8)
    real(kind=r_def), intent(out) :: parallel_factor, perp_factor
    real(kind=r_def), intent(out) :: eave_lon, eave_lat, correct_flux

    integer(kind=i_def) :: side_edge_dist(4), side_edge_idx
    integer(kind=i_def) :: edge_direction, swapped_panel
    logical(kind=l_def) :: is_x_direction, near_panel_edge

    weights = 0.0_r_def
    offsets = 0_i_def
    perp_weights = 0.0_r_def
    perp_offsets = 0_i_def
    parallel_factor = 1.0_r_def
    perp_factor = 0.0_r_def
    eave_lon = 0.0_r_def
    eave_lat = 0.0_r_def
    correct_flux = 0.0_r_def

    is_x_direction = (target_df == S .or. target_df == N)
    if (target_across) then
      if (edge_stencil_size(df) < 2) return
      if (is_x_direction) then
        if (coord_x_stencil_size(df) < 2) return
      else
        if (coord_y_stencil_size(df) < 2) return
      end if
      side_edge_idx = edge_stencil_map(1, 2, df)
      side_edge_dist = (/ panel_edge_dist_W(side_edge_idx),                  &
                          panel_edge_dist_E(side_edge_idx),                  &
                          panel_edge_dist_S(side_edge_idx),                  &
                          panel_edge_dist_N(side_edge_idx) /)
    else
      side_edge_dist = edge_dist
    end if

    if (is_x_direction) then
      if (ABS(side_edge_dist(1)) < ABS(side_edge_dist(2))) then
        edge_direction = W
      else
        edge_direction = E
      end if
      near_panel_edge = MIN(ABS(side_edge_dist(1)), ABS(side_edge_dist(2)))    &
                        <= remap_depth
    else
      if (ABS(side_edge_dist(3)) < ABS(side_edge_dist(4))) then
        edge_direction = S
      else
        edge_direction = N
      end if
      near_panel_edge = MIN(ABS(side_edge_dist(3)), ABS(side_edge_dist(4)))    &
                        <= remap_depth
    end if

    if (.not. near_panel_edge) return

    swapped_panel = panel_neighbour(target_panel, edge_direction)
    call compute_lagrange_weights(target_df, target_df, is_x_direction,        &
        swapped_panel, lagrange_w, lagrange_id,                                &
        target_panel_arg=target_panel, host_face=df, host_sign=host_sign,      &
        target_across=target_across)
    if (ANY(lagrange_id(1:ndata) == 0_i_def)) return
    weights(1:ndata) = lagrange_w(1:ndata)
    offsets(1:ndata) = lagrange_id(1:ndata)

    call compute_rotation_factors(target_df, is_x_direction, swapped_panel,   &
        target_panel, df, target_across, parallel_factor, perp_factor,        &
        eave_lon, eave_lat, correct_flux)

    call compute_perp_weights(target_df, is_x_direction, target_panel,        &
      swapped_panel, edge_direction, host_sign, target_across)
    perp_weights(1:ndata) = perp_w_a(1:ndata)
    perp_weights(ndata+1:2*ndata) = perp_w_b(1:ndata)
    perp_offsets(1:ndata) = perp_id_a(1:ndata)
    perp_offsets(ndata+1:2*ndata) = perp_id_b(1:ndata)

  end subroutine compute_labelled_side

  subroutine compute_perp_weights(target_df, is_x_direction, target_panel,    &
                                  swapped_panel, edge_direction, host_sign,   &
                                  target_across)

    implicit none

    integer(kind=i_def), intent(in) :: target_df, target_panel
    integer(kind=i_def), intent(in) :: swapped_panel, edge_direction, host_sign
    logical(kind=l_def), intent(in) :: is_x_direction, target_across

    integer(kind=i_def) :: corner_a, corner_b, source_a, source_b

    perp_w_a = 0.0_r_def
    perp_w_b = 0.0_r_def
    perp_id_a = 0_i_def
    perp_id_b = 0_i_def

    if (is_x_direction) then
      source_a = W
      source_b = E
      if (target_df == S) then
        corner_a = SW
        corner_b = SE
      else
        corner_a = NW
        corner_b = NE
      end if
    else
      source_a = S
      source_b = N
      if (target_df == W) then
        corner_a = SW
        corner_b = NW
      else
        corner_a = SE
        corner_b = NE
      end if
    end if

    call compute_lagrange_weights(corner_a - SW + 1, source_a,                &
        is_x_direction, swapped_panel, perp_w_a, perp_id_a,                    &
        is_corner_target=.true., target_panel_arg=target_panel,                &
        host_face=df, host_sign=host_sign, target_across=target_across)
    call compute_lagrange_weights(corner_b - SW + 1, source_b,                &
        is_x_direction, swapped_panel, perp_w_b, perp_id_b,                    &
        is_corner_target=.true., target_panel_arg=target_panel,                &
        host_face=df, host_sign=host_sign, target_across=target_across)
  end subroutine compute_perp_weights

  subroutine store_labelled_fields(label, weights_idx, rotation_idx, perp_idx,&
                                   weights, offsets, perp_weights,            &
                                   perp_offsets, parallel_factor, perp_factor,&
                                   eave_lon, eave_lat, correct_flux)

    implicit none

    integer(kind=i_def), intent(in) :: label, weights_idx, rotation_idx, perp_idx
    real(kind=r_def), intent(in) :: weights(4), perp_weights(8)
    integer(kind=i_def), intent(in) :: offsets(4), perp_offsets(8)
    real(kind=r_def), intent(in) :: parallel_factor, perp_factor
    real(kind=r_def), intent(in) :: eave_lon, eave_lat, correct_flux

    select case (label)
    case (1)
      interp_w2h_weights_1(weights_idx:weights_idx+ndata-1) =                  &
          REAL(weights(1:ndata), r_tran)
      interp_w2h_indices_1(weights_idx:weights_idx+ndata-1) =                  &
          offsets(1:ndata)
      parl_w2_rot_1(rotation_idx) = REAL(parallel_factor, r_tran)
      perp_w2_rot_1(rotation_idx) = REAL(perp_factor, r_tran)
      interp_perp_weights_1(perp_idx:perp_idx+2*ndata-1) =                    &
          REAL(perp_weights(1:2*ndata), r_tran)
      interp_perp_indices_1(perp_idx:perp_idx+2*ndata-1) =                    &
          perp_offsets(1:2*ndata)
      eave_lon_1(weights_idx) = REAL(eave_lon, r_tran)
      eave_lat_1(weights_idx) = REAL(eave_lat, r_tran)
      correct_flux_1(weights_idx) = REAL(correct_flux, r_tran)
    case (2)
      interp_w2h_weights_2(weights_idx:weights_idx+ndata-1) =                  &
          REAL(weights(1:ndata), r_tran)
      interp_w2h_indices_2(weights_idx:weights_idx+ndata-1) =                  &
          offsets(1:ndata)
      parl_w2_rot_2(rotation_idx) = REAL(parallel_factor, r_tran)
      perp_w2_rot_2(rotation_idx) = REAL(perp_factor, r_tran)
      interp_perp_weights_2(perp_idx:perp_idx+2*ndata-1) =                    &
          REAL(perp_weights(1:2*ndata), r_tran)
      interp_perp_indices_2(perp_idx:perp_idx+2*ndata-1) =                    &
          perp_offsets(1:2*ndata)
      eave_lon_2(weights_idx) = REAL(eave_lon, r_tran)
      eave_lat_2(weights_idx) = REAL(eave_lat, r_tran)
      correct_flux_2(weights_idx) = REAL(correct_flux, r_tran)
    case (3)
      interp_w2h_weights_3(weights_idx:weights_idx+ndata-1) =                  &
          REAL(weights(1:ndata), r_tran)
      interp_w2h_indices_3(weights_idx:weights_idx+ndata-1) =                  &
          offsets(1:ndata)
      parl_w2_rot_3(rotation_idx) = REAL(parallel_factor, r_tran)
      perp_w2_rot_3(rotation_idx) = REAL(perp_factor, r_tran)
      interp_perp_weights_3(perp_idx:perp_idx+2*ndata-1) =                    &
          REAL(perp_weights(1:2*ndata), r_tran)
      interp_perp_indices_3(perp_idx:perp_idx+2*ndata-1) =                    &
          perp_offsets(1:2*ndata)
      eave_lon_3(weights_idx) = REAL(eave_lon, r_tran)
      eave_lat_3(weights_idx) = REAL(eave_lat, r_tran)
      correct_flux_3(weights_idx) = REAL(correct_flux, r_tran)
    end select

  end subroutine store_labelled_fields

  integer(kind=i_def) function face_direction_to_panel(panel, target_panel)    &
      result(direction)

    implicit none

    integer(kind=i_def), intent(in) :: panel, target_panel
    integer(kind=i_def) :: directions(4), k

    directions = (/ W, S, E, N /)
    direction = 0_i_def
    do k = 1, 4
      if (panel_neighbour(panel, directions(k)) == target_panel) then
        direction = directions(k)
        return
      end if
    end do

  end function face_direction_to_panel

  integer(kind=i_def) function opposite_face(direction) result(opposite)

    implicit none

    integer(kind=i_def), intent(in) :: direction

    select case (direction)
    case (W)
      opposite = E
    case (E)
      opposite = W
    case (S)
      opposite = N
    case (N)
      opposite = S
    end select

  end function opposite_face

  !-----------------------------------------------------------------------------
  !> @brief Compute the Lagrange donor-stencil interpolation weights/indices
  !!        tangential to the panel edge, for the target W2H DoF "target_df".
  !!        This mirrors panel_edge_weights_kernel_mod's panel_edge_weights_1d,
  !!        but evaluated at a W2H face node instead of a W3/Wtheta cell node.
  !!        "source_df" is the face DoF whose values will actually be
  !!        interpolated (i.e. it determines the donor node positions sampled
  !!        in each donor cell), kept separate from "target_df" (which only
  !!        determines the interpolation point x0). For this phase's
  !!        "parl_w2" interpolation, source_df == target_df; a later phase's
  !!        "perp_w2" interpolation will pass a different source_df.
  !!        "is_corner_target", if present and .true., indicates "target_df"
  !!        is a W1V vertex (corner) index rather than a W2H face index, so
  !!        the interpolation point x0 is evaluated using the "wv" basis
  !!        arrays instead of the "ww" ones. The donor node positions are
  !!        always sampled at the source's own W2H face node ("ww" basis),
  !!        regardless of the target.
  subroutine compute_lagrange_weights(target_df, source_df, is_x_direction,    &
                                       swapped_panel, wgt, id, is_corner_target,&
                                       is_w3_target, target_panel_arg,         &
                                       host_face, host_sign, target_across)

    implicit none

    integer(kind=i_def), intent(in)  :: target_df
    integer(kind=i_def), intent(in)  :: source_df
    logical(kind=l_def), intent(in)  :: is_x_direction
    integer(kind=i_def), intent(in)  :: swapped_panel
    real(kind=r_def),    intent(out) :: wgt(4)
    integer(kind=i_def), intent(out) :: id(4)
    logical(kind=l_def), intent(in), optional :: is_corner_target
    logical(kind=l_def), intent(in), optional :: is_w3_target
    integer(kind=i_def), intent(in), optional :: target_panel_arg
    integer(kind=i_def), intent(in), optional :: host_face, host_sign
    logical(kind=l_def), intent(in), optional :: target_across

    integer(kind=i_def) :: branch_l, branch_r
    integer(kind=i_def) :: i, panel_for_target
    integer(kind=i_def) :: wx_size_l, wx_size_r
    integer(kind=i_def) :: pid_size_l, pid_size_r
    integer(kind=i_def) :: ncells_in_stencil
    integer(kind=i_def) :: stencil_extent, other_branch, stencil_branch
    integer(kind=i_def) :: stencil_position, offset, candidate_count
    integer(kind=i_def) :: coord_idx
    integer(kind=i_def) :: wx_stencil_1d(ndf_wx, 2*wx_max_length-1)
    integer(kind=i_def) :: pid_stencil_1d(2*pid_max_length-1)
    integer(kind=i_def) :: donor_offsets(2*wx_max_length-1)
    integer(kind=i_def) :: selected_index(4)
    real(kind=r_def)    :: x1(2*wx_max_length-1), dx(2*wx_max_length-1)
    real(kind=r_def)    :: abh(3), xyz(3), h_dummy
    real(kind=r_def)    :: alpha_owned, beta_owned
    real(kind=r_def)    :: x0
    integer(kind=i_def) :: interp_dir, panel_edge, panel_of_donor
    integer(kind=i_def) :: id1, id2, id3, id4

    ! Tangential-to-edge stencil branches (S/N for x remap, W/E for y remap)
    if (is_x_direction) then
      branch_l = S
      branch_r = N
    else
      branch_l = W
      branch_r = E
    end if

    wx_size_l = wx_stencil_size(branch_l)
    wx_size_r = wx_stencil_size(branch_r)
    pid_size_l = pid_stencil_size(branch_l)
    pid_size_r = pid_stencil_size(branch_r)

    ncells_in_stencil = wx_size_l + wx_size_r - 1

    wx_stencil_1d = 0
    do i = 1, wx_size_l
      wx_stencil_1d(:,i) = wx_stencil_map(:,i,branch_l)
    end do
    do i = 1, wx_size_r - 1
      wx_stencil_1d(:,i+wx_size_l) = wx_stencil_map(:,i+1,branch_r)
    end do

    pid_stencil_1d = 0
    do i = 1, pid_size_l
      pid_stencil_1d(i) = pid_stencil_map(1,i,branch_l)
    end do
    do i = 1, pid_size_r - 1
      pid_stencil_1d(i+pid_size_l) = pid_stencil_map(1,i+1,branch_r)
    end do

    panel_for_target = owned_panel
    if (present(target_panel_arg)) panel_for_target = target_panel_arg

    if (present(host_face)) then
      other_branch = opposite_face(host_face)
      stencil_extent = MIN(                                                &
          MIN(wx_stencil_size(host_face), wx_stencil_size(other_branch)),  &
          MIN(pid_stencil_size(host_face), pid_stencil_size(other_branch)) &
      ) - 1
      candidate_count = 0
      wx_stencil_1d = 0_i_def
      pid_stencil_1d = 0_i_def
      donor_offsets = 0_i_def
      do offset = -stencil_extent, stencil_extent
        if (offset == 0) cycle
        if (offset == host_sign) then
          stencil_branch = 0_i_def
          stencil_position = 1_i_def
        else if (offset*host_sign < 0_i_def) then
          stencil_branch = host_face
          stencil_position = ABS(offset) + 1
        else
          stencil_branch = other_branch
          stencil_position = ABS(offset)
        end if
        if (stencil_branch == 0_i_def) then
          panel_of_donor = int(panel_id(map_pid(1)), i_def)
          if (panel_of_donor /= panel_for_target) cycle
          candidate_count = candidate_count + 1
          wx_stencil_1d(:, candidate_count) = map_wx
          pid_stencil_1d(candidate_count) = map_pid(1)
        else
          if (stencil_position > wx_stencil_size(stencil_branch) .or.        &
              stencil_position > pid_stencil_size(stencil_branch)) cycle
          panel_of_donor = int(                                             &
              panel_id(pid_stencil_map(1, stencil_position, stencil_branch)), &
              i_def)
          if (panel_of_donor /= panel_for_target) cycle
          candidate_count = candidate_count + 1
          wx_stencil_1d(:, candidate_count) =                                &
              wx_stencil_map(:, stencil_position, stencil_branch)
          pid_stencil_1d(candidate_count) =                                  &
              pid_stencil_map(1, stencil_position, stencil_branch)
        end if
        donor_offsets(candidate_count) = offset
      end do
      ncells_in_stencil = candidate_count
      wgt = 0.0_r_def
      id = 0_i_def
      if (candidate_count < ndata) return
    end if

    ! Direction of interpolation (alpha or beta), based on owned/swapped panels
    panel_edge = 10*panel_for_target + swapped_panel
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

    ! Interpolation point x0: the target DoF's own extended-mesh coords,
    ! converted back into the owned panel's native coordinates. A "face"
    ! target (parl_w2, set 1) is evaluated at the "ww" nodes; a "corner"
    ! target (perp_w2, set 2) is evaluated at the "wv" (W1V vertex) nodes
    abh = 0.0_r_def
    if (present(is_w3_target) .and. is_w3_target) then
      if (is_x_direction) then
        do dfc = 1, ndf_wx_2d
          coord_idx = map_wx_2d(dfc)
          abh(1) = abh(1) + alpha_x(coord_idx)*basis_wx_2d_w3(1,dfc,target_df)
          abh(2) = abh(2) + beta_x(coord_idx)*basis_wx_2d_w3(1,dfc,target_df)
        end do
      else
        do dfc = 1, ndf_wx_2d
          coord_idx = map_wx_2d(dfc)
          abh(1) = abh(1) + alpha_y(coord_idx)*basis_wx_2d_w3(1,dfc,target_df)
          abh(2) = abh(2) + beta_y(coord_idx)*basis_wx_2d_w3(1,dfc,target_df)
        end do
      end if
    else if (present(is_corner_target) .and. is_corner_target) then
      if (is_x_direction) then
        do dfc = 1, ndf_wx_2d
          coord_idx = map_wx_2d(dfc)
          if (present(target_across)) then
            if (target_across) coord_idx =                                  &
                coord_x_stencil_map(dfc, 2, host_face)
          end if
          abh(1) = abh(1) + alpha_x(coord_idx)*basis_wx_2d_wv(1,dfc,target_df)
          abh(2) = abh(2) + beta_x(coord_idx)*basis_wx_2d_wv(1,dfc,target_df)
        end do
      else
        do dfc = 1, ndf_wx_2d
          coord_idx = map_wx_2d(dfc)
          if (present(target_across)) then
            if (target_across) coord_idx =                                  &
                coord_y_stencil_map(dfc, 2, host_face)
          end if
          abh(1) = abh(1) + alpha_y(coord_idx)*basis_wx_2d_wv(1,dfc,target_df)
          abh(2) = abh(2) + beta_y(coord_idx)*basis_wx_2d_wv(1,dfc,target_df)
        end do
      end if
    else
      if (is_x_direction) then
        do dfc = 1, ndf_wx_2d
          coord_idx = map_wx_2d(dfc)
          if (present(target_across)) then
            if (target_across) coord_idx =                                  &
                coord_x_stencil_map(dfc, 2, host_face)
          end if
          abh(1) = abh(1) + alpha_x(coord_idx)*basis_wx_2d_ww(1,dfc,target_df)
          abh(2) = abh(2) + beta_x(coord_idx)*basis_wx_2d_ww(1,dfc,target_df)
        end do
      else
        do dfc = 1, ndf_wx_2d
          coord_idx = map_wx_2d(dfc)
          if (present(target_across)) then
            if (target_across) coord_idx =                                  &
                coord_y_stencil_map(dfc, 2, host_face)
          end if
          abh(1) = abh(1) + alpha_y(coord_idx)*basis_wx_2d_ww(1,dfc,target_df)
          abh(2) = abh(2) + beta_y(coord_idx)*basis_wx_2d_ww(1,dfc,target_df)
        end do
      end if
    end if

    call alphabetar2xyz(                                                       &
            abh(1), abh(2), unit_radius, swapped_panel, xyz(1), xyz(2), xyz(3) &
    )
    call xyz2alphabetar(                                                       &
            xyz(1), xyz(2), xyz(3), panel_for_target, abh(1), abh(2), h_dummy  &
    )
    x0 = abh(interp_dir)

    ! Locations of donor points, in owned-panel native coordinates. Donor
    ! node positions are sampled at the SOURCE DoF's own face node of each
    ! donor cell (see the subroutine's docstring)
    do i = 1, ncells_in_stencil

      panel_of_donor = int(panel_id(pid_stencil_1d(i)), i_def)

      if ( panel_for_target /= panel_of_donor ) then
        x1(i) = LARGE_REAL_POSITIVE
      else
        abh = 0.0_r_def
        do dfc = 1, ndf_wx
          call chi2xyz(                                                        &
                  chi_1(wx_stencil_1d(dfc,i)),                                 &
                  chi_2(wx_stencil_1d(dfc,i)),                                 &
                  unit_radius, panel_for_target, geometry, topology,           &
                  coord_system,                                                &
                  scaled_radius, xyz(1), xyz(2), xyz(3)                       &
          )
          if (to_rotate) then
            xyz = matmul(inverse_rot_matrix, xyz)
          end if
          if (to_stretch) then
            xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
          end if
          call xyz2alphabetar(                                                 &
                  xyz(1), xyz(2), xyz(3), panel_for_target,                   &
                  alpha_owned, beta_owned, h_dummy                            &
          )
          abh(1) = abh(1) + alpha_owned * basis_wx_ww(1,dfc,source_df)
          abh(2) = abh(2) + beta_owned * basis_wx_ww(1,dfc,source_df)
        end do
        x1(i) = abh(interp_dir)
      end if
    end do

    ! Lagrange interpolation weights and donor indices
    dx = LARGE_REAL_POSITIVE
    dx(1:ncells_in_stencil) = ABS(x1(1:ncells_in_stencil) - x0)
    id1 = MINLOC(dx(1:ncells_in_stencil), 1)
    dx(id1) = LARGE_REAL_POSITIVE
    id2 = MINLOC(dx(1:ncells_in_stencil), 1)

    wgt = 0.0_r_def
    id = 1
    selected_index = 1_i_def

    if (cubic_remap) then
      dx(id2) = LARGE_REAL_POSITIVE
      id3 = MINLOC(dx(1:ncells_in_stencil), 1)
      dx(id3) = LARGE_REAL_POSITIVE
      id4 = MINLOC(dx(1:ncells_in_stencil), 1)

      wgt(1) = (x0 - x1(id2))/(x1(id1) - x1(id2))                              &
            * (x0 - x1(id3))/(x1(id1) - x1(id3))                               &
            * (x0 - x1(id4))/(x1(id1) - x1(id4))
      wgt(2) = (x0 - x1(id1))/(x1(id2) - x1(id1))                              &
            * (x0 - x1(id3))/(x1(id2) - x1(id3))                               &
            * (x0 - x1(id4))/(x1(id2) - x1(id4))
      wgt(3) = (x0 - x1(id1))/(x1(id3) - x1(id1))                              &
            * (x0 - x1(id2))/(x1(id3) - x1(id2))                               &
            * (x0 - x1(id4))/(x1(id3) - x1(id4))
      wgt(4) = (x0 - x1(id1))/(x1(id4) - x1(id1))                              &
            * (x0 - x1(id2))/(x1(id4) - x1(id2))                               &
            * (x0 - x1(id3))/(x1(id4) - x1(id3))

      selected_index = (/ id1, id2, id3, id4 /)
    else
      wgt(1) = (x0 - x1(id2))/(x1(id1) - x1(id2))
      wgt(2) = (x0 - x1(id1))/(x1(id2) - x1(id1))

      selected_index(1:2) = (/ id1, id2 /)
    end if

    if (present(host_face)) then
      id = 0_i_def
      id(1:ndata) = donor_offsets(selected_index(1:ndata))
    else
      id = selected_index
    end if

  end subroutine compute_lagrange_weights

  !-----------------------------------------------------------------------------
  !> @brief Unit normals to the alpha=const and beta=const grid lines of a
  !!        panel at (alpha, beta), as global Cartesian vectors. The signs are
  !!        such that each normal has a positive component along the
  !!        corresponding alpha or beta tangent.
  subroutine compute_face_normals(alpha, beta, panel, n_alpha, n_beta)

    implicit none

    real(kind=r_def),    intent(in)  :: alpha, beta
    integer(kind=i_def), intent(in)  :: panel
    real(kind=r_def),    intent(out) :: n_alpha(3), n_beta(3)

    real(kind=r_def) :: abr(3), t_alpha(3), t_beta(3), r_hat(3)

    abr = (/ alpha, beta, unit_radius /)
    t_alpha = alphabetar2xyz_vector((/ 1.0_r_def, 0.0_r_def, 0.0_r_def /),     &
                                    abr, panel)
    t_beta = alphabetar2xyz_vector((/ 0.0_r_def, 1.0_r_def, 0.0_r_def /),      &
                                   abr, panel)
    r_hat = alphabetar2xyz_vector((/ 0.0_r_def, 0.0_r_def, 1.0_r_def /),       &
                                  abr, panel)

    n_alpha = cross_product(t_beta, r_hat)
    n_alpha = SIGN(1.0_r_def, DOT_PRODUCT(n_alpha, t_alpha)) * n_alpha         &
              / SQRT(SUM(n_alpha**2))
    n_beta = cross_product(r_hat, t_alpha)
    n_beta = SIGN(1.0_r_def, DOT_PRODUCT(n_beta, t_beta)) * n_beta             &
             / SQRT(SUM(n_beta**2))

  end subroutine compute_face_normals

  pure function cross_product(a, b) result(c)

    implicit none

    real(kind=r_def), intent(in) :: a(3), b(3)
    real(kind=r_def)             :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_product

  !-----------------------------------------------------------------------------
  !> @brief Compute the factors that convert the "parl_w2" and "perp_w2"
  !!        values (P-panel face values of the slot's own and of the
  !!        cross orientation) into the value on the neighbour (eave) panel's
  !!        extended-mesh face stored in this DoF slot. The eave face normal is
  !!        decomposed in the two P-panel face normals, and the stored sign
  !!        convention (W/E and alpha-type faces positive, S/N and beta-type
  !!        faces negative) is applied. On rotated panel edges the slot's eave
  !!        face is of the other (alpha/beta) type. Also computes the analytic
  !!        solid-body-rotation flux through the eave face, for diagnostics.
  subroutine compute_rotation_factors(target_df, is_x_direction,               &
                                       neighbour_panel_id, target_panel,       &
                                       host_face, target_across,              &
                                       parallel_factor, perp_factor,           &
                                       eave_lon_out, eave_lat_out,             &
                                       correct_flux_out)

    implicit none

    integer(kind=i_def), intent(in)  :: target_df
    logical(kind=l_def), intent(in)  :: is_x_direction
    integer(kind=i_def), intent(in)  :: neighbour_panel_id
    integer(kind=i_def), intent(in)  :: target_panel, host_face
    logical(kind=l_def), intent(in)  :: target_across
    real(kind=r_def),    intent(out) :: parallel_factor, perp_factor
    real(kind=r_def),    intent(out) :: eave_lon_out, eave_lat_out
    real(kind=r_def),    intent(out) :: correct_flux_out

    integer(kind=i_def) :: dfc, idx
    real(kind=r_def)    :: chi_1_df, chi_2_df, chi_3_df
    real(kind=r_def)    :: alpha_native, beta_native
    real(kind=r_def)    :: alpha_eave, beta_eave
    real(kind=r_def)    :: xyz(3), eave_xyz_true(3), r
    real(kind=r_def)    :: alpha_own, beta_own
    real(kind=r_def)    :: n_alpha_own(3), n_beta_own(3)
    real(kind=r_def)    :: n_alpha_nbr(3), n_beta_nbr(3)
    real(kind=r_def)    :: n_par(3), n_perp(3), n_q(3), normal_llr(3)
    real(kind=r_def)    :: gram, rhs_par, rhs_perp, c_par, c_perp, s_par, s_q, s_axis
    logical(kind=l_def) :: slot_is_ns, rotated_edge, q_is_alpha

    chi_1_df = 0.0_r_def
    chi_2_df = 0.0_r_def
    chi_3_df = 0.0_r_def
    do dfc = 1, ndf_wx
      if (target_across) then
        idx = wx_stencil_map(dfc, 2, host_face)
      else
        idx = map_wx(dfc)
      end if
      chi_1_df = chi_1_df + chi_1(idx) * basis_wx_ww(1, dfc, target_df)
      chi_2_df = chi_2_df + chi_2(idx) * basis_wx_ww(1, dfc, target_df)
      chi_3_df = chi_3_df + chi_3(idx) * basis_wx_ww(1, dfc, target_df)
    end do

    call chi2xyz(                                                              &
        chi_1_df, chi_2_df, chi_3_df, target_panel,                            &
        geometry, topology, coord_system, scaled_radius,                       &
        xyz(1), xyz(2), xyz(3)                                                 &
    )
    if (to_rotate) then
      xyz = matmul(inverse_rot_matrix, xyz)
    end if
    if (to_stretch) then
      xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
    end if
    call xyz2alphabetar(                                                       &
      xyz(1), xyz(2), xyz(3), target_panel, alpha_native, beta_native, r     &
    )

    alpha_eave = 0.0_r_def
    beta_eave = 0.0_r_def
    if (is_x_direction) then
      do dfc = 1, ndf_wx_2d
        if (target_across) then
          idx = coord_x_stencil_map(dfc, 2, host_face)
        else
          idx = map_wx_2d(dfc)
        end if
        alpha_eave = alpha_eave + alpha_x(idx) * basis_wx_2d_ww(1, dfc, target_df)
        beta_eave = beta_eave + beta_x(idx) * basis_wx_2d_ww(1, dfc, target_df)
      end do
    else
      do dfc = 1, ndf_wx_2d
        if (target_across) then
          idx = coord_y_stencil_map(dfc, 2, host_face)
        else
          idx = map_wx_2d(dfc)
        end if
        alpha_eave = alpha_eave + alpha_y(idx) * basis_wx_2d_ww(1, dfc, target_df)
        beta_eave = beta_eave + beta_y(idx) * basis_wx_2d_ww(1, dfc, target_df)
      end do
    end if
    ! Position of the eave face DoF on the sphere, and its coordinates in the
    ! owned panel's frame (the frame of the donor values being interpolated)
    call alphabetar2xyz(                                                       &
        alpha_eave, beta_eave, 1.0_r_def, neighbour_panel_id,                  &
        xyz(1), xyz(2), xyz(3)                                                 &
    )
    call xyz2alphabetar(                                                       &
      xyz(1), xyz(2), xyz(3), target_panel, alpha_own, beta_own, r           &
    )

    ! True eave lon/lat, obtained by applying the same rotate/stretch
    ! corrections used above for the target DoF's own native position
    eave_xyz_true = xyz
    if (to_rotate) then
      eave_xyz_true = matmul(inverse_rot_matrix, eave_xyz_true)
    end if
    if (to_stretch) then
      eave_xyz_true = inverse_schmidt_transform_xyz(eave_xyz_true, stretch_factor)
    end if
    call xyz2ll(eave_xyz_true(1), eave_xyz_true(2), eave_xyz_true(3),          &
                eave_lon_out, eave_lat_out)

    call compute_face_normals(alpha_own, beta_own, target_panel,               &
                              n_alpha_own, n_beta_own)
    call compute_face_normals(alpha_eave, beta_eave, neighbour_panel_id,       &
                              n_alpha_nbr, n_beta_nbr)

    ! On rotated panel edges the extended mesh is rotated, so the eave face
    ! stored in an S/N slot is an alpha-type face (and W/E slots beta-type)
    select case (10*target_panel + neighbour_panel_id)
    case (14, 23, 35, 46, 52, 61, 16, 25, 32, 41, 53, 64)
      rotated_edge = .true.
    case default
      rotated_edge = .false.
    end select

    slot_is_ns = (target_df == S .or. target_df == N)
    if (slot_is_ns) then
      n_par = n_beta_own
      n_perp = n_alpha_own
      s_par = -1.0_r_def
    else
      n_par = n_alpha_own
      n_perp = n_beta_own
      s_par = 1.0_r_def
    end if

    q_is_alpha = (slot_is_ns .eqv. rotated_edge)
    if (q_is_alpha) then
      n_q = n_alpha_nbr
      s_q = 1.0_r_def
    else
      n_q = n_beta_nbr
      s_q = -1.0_r_def
    end if

    ! Decompose the eave face normal in the two (non-orthogonal) owned-panel
    ! face normals. The cross-orientation face has the opposite stored sign
    gram = DOT_PRODUCT(n_par, n_perp)
    rhs_par = DOT_PRODUCT(n_par, n_q)
    rhs_perp = DOT_PRODUCT(n_perp, n_q)
    c_par = (rhs_par - gram*rhs_perp) / (1.0_r_def - gram*gram)
    c_perp = (rhs_perp - gram*rhs_par) / (1.0_r_def - gram*gram)
    ! Re-express the stored value in the host slot's sign convention: the
    ! consumers read it as positive along the host line, not the eave axis
    s_axis = s_par * s_q * SIGN(1.0_r_def, c_par)
    parallel_factor = s_axis * c_par * s_q * s_par
    perp_factor = -s_axis * c_perp * s_q * s_par

    ! Diagnostic "correct" flux: the analytic solid-body-rotation wind
    ! (u0*cos(lat), v=0) through the eave face, in the stored sign convention.
    ! Mesh stretching is not accounted for in this vector projection (only in
    ! the position above)
    normal_llr = cart2sphere_vector(eave_xyz_true, n_q)
    correct_flux_out = s_axis * s_q * u0 * COS(eave_lat_out) * normal_llr(1) * dt &
              * dA_at_w2(map_w2(host_face))

  end subroutine compute_rotation_factors

end subroutine panel_edge_remap_w2h_weights_code

end module panel_edge_remap_w2h_weights_kernel_mod
