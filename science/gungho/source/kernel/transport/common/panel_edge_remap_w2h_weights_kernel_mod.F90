!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute the interpolation weights and indices needed to remap a W2H
!!        flux/wind field from the standard cubed-sphere mesh onto the
!!        "extended"/eave mesh at the edges of the mesh panels.
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
!!          Two independent remaps are computed:
!!            - the "x" remap, triggered by proximity to a W or E panel edge,
!!              which corrects the S/N (v-type) face values of a cell,
!!            - the "y" remap, triggered by proximity to a S or N panel edge,
!!              which corrects the W/E (u-type) face values of a cell.
!!          This mirrors the split into panel_edge_weights_x/y used for
!!          scalar fields.
!!
!!          For each remap direction, two donor-weight fields are produced:
!!            - the "parallel" weight, which interpolates donor values taken
!!              from the SAME face orientation as the target DoF (e.g. donor
!!              S/N values contribute to a target S/N value), using a
!!              Lagrange interpolation across the donor stencil tangential
!!              to the edge (exactly as in panel_edge_weights_kernel_mod),
!!              additionally scaled by the cos(eave)/cos(native) rotation
!!              factor used previously in peregrin_weights_kernel_mod,
!!            - the "perp" weight, which interpolates donor values taken from
!!              the OPPOSITE face orientation (e.g. donor W/E values
!!              contribute to a target S/N value), using the SAME donor
!!              stencil/indices as the parallel weight, scaled by the
!!              -sin(eave)/cos(native) rotation factor.
!!          Only ONE index field per direction (indices_x, indices_y) is
!!          needed, since the donor CELL is the same for both the parallel
!!          and perp weights of that direction; only the choice of which
!!          face DoF within that donor cell is sampled differs (this is
!!          determined by the remap kernel, not stored here).
!!
!!          Away from any relevant panel edge, all weights default to the
!!          identity (parallel = 1 at the target's own DoF, perp = 0), so
!!          that applying the remap leaves the field unchanged there.
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
                                  ANY_DISCONTINUOUS_SPACE_9,                   &
                                  ANY_SPACE_9,                                 &
                                  CELL_COLUMN, STENCIL, CROSS2D
use fs_continuity_mod,       only: W2h
use constants_mod,           only: r_tran, r_def, i_def, l_def,               &
                                   LARGE_REAL_POSITIVE
use reference_element_mod,   only: W, S, N, E
use panel_edge_support_mod,  only: panel_neighbour
use sci_face_selector_support_mod, only: face_from_face_selector
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_DEBUG

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!! PSy layer
type, public, extends(kernel_type) :: panel_edge_remap_w2h_weights_kernel_type
  private
  type(arg_type) :: meta_args(19) = (/                                         &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! weights_x_parallel
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! weights_x_perp
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! weights_y_parallel
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! weights_y_perp
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! indices_x
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! indices_y
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9,                &
                                                        STENCIL(CROSS2D)),     & ! chi
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9), & ! alpha_x, beta_x
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9), & ! alpha_y, beta_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,  &
                                                        STENCIL(CROSS2D)),     & ! panel_id
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_edge_dist_W/E/S/N
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ew
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ns
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                             & ! cubic_remap
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! remap_depth
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! geometry
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! topology
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! coord_system
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              & ! scaled_radius
  /)
  type(func_type) :: meta_funcs(2) = (/                                        &
      func_type(ANY_SPACE_9, GH_BASIS),                                        &
      func_type(ANY_DISCONTINUOUS_SPACE_9, GH_BASIS)                           &
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: panel_edge_remap_w2h_weights_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_remap_w2h_weights_code

contains

!> @brief Compute the interpolation weights/indices to remap a W2H field
!!        across cubed-sphere panel edges.
!> @param[in]     nlayers               Number of layers
!> @param[in,out] weights_x_parallel    Same-orientation donor weights for x remap
!> @param[in,out] weights_x_perp        Cross-orientation donor weights for x remap
!> @param[in,out] weights_y_parallel    Same-orientation donor weights for y remap
!> @param[in,out] weights_y_perp        Cross-orientation donor weights for y remap
!> @param[in,out] indices_x             Donor cell indices for the x remap
!> @param[in,out] indices_y             Donor cell indices for the y remap
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
!> @param[in]     ndf_ww                Num DoFs per cell for W2H
!> @param[in]     undf_ww               Num W2H DoFs for this partition
!> @param[in]     map_ww                DoF map for W2H
!> @param[in]     ndf_wx                Num DoFs per cell for coordinate fields
!> @param[in]     undf_wx               Num coordinate DoFs for this partition
!> @param[in]     map_wx                DoF map for coordinate fields
!> @param[in]     basis_wx              Basis functions of chi, evaluated at W2H nodes
!> @param[in]     ndf_wx_2d             Num DoFs per cell for 2D coord fields
!> @param[in]     undf_wx_2d            Num 2D coord DoFs for this partition
!> @param[in]     map_wx_2d             DoF map for 2D coord fields
!> @param[in]     basis_wx_2d           Basis functions of extended coords, at W2H nodes
!> @param[in]     ndf_pid               Num DoFs per cell for panel ID field
!> @param[in]     undf_pid              Num DoFs for this partition for panel ID
!> @param[in]     map_pid               DoF map for panel ID field
subroutine panel_edge_remap_w2h_weights_code(                                  &
    nlayers,                                                                   &
    weights_x_parallel, weights_x_perp,                                        &
    weights_y_parallel, weights_y_perp,                                        &
    indices_x, indices_y,                                                      &
    chi_1, chi_2, chi_3,                                                       &
    wx_stencil_size, wx_max_length, wx_stencil_map,                            &
    alpha_x, beta_x, alpha_y, beta_y,                                          &
    panel_id,                                                                  &
    pid_stencil_size, pid_max_length, pid_stencil_map,                         &
    panel_edge_dist_W, panel_edge_dist_E,                                      &
    panel_edge_dist_S, panel_edge_dist_N,                                      &
    face_selector_ew, face_selector_ns,                                        &
    cubic_remap, remap_depth,                                                  &
    geometry, topology, coord_system, scaled_radius,                           &
    ndf_ww, undf_ww, map_ww,                                                   &
    ndf_wx, undf_wx, map_wx, basis_wx,                                         &
    ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,                             &
    ndf_pid, undf_pid, map_pid                                                 &
)

  use sci_chi_transform_mod, only: chi2xyz, get_to_rotate, get_to_stretch,     &
                                   get_inverse_mesh_rotation_matrix,           &
                                   get_stretch_factor
  use coord_transform_mod,   only: alphabetar2xyz, xyz2alphabetar,             &
                                   inverse_schmidt_transform_xyz

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
  integer(kind=i_def), intent(in)    :: geometry
  integer(kind=i_def), intent(in)    :: topology
  integer(kind=i_def), intent(in)    :: coord_system
  real(kind=r_def),    intent(in)    :: scaled_radius

  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wx_stencil_map(ndf_wx, wx_max_length, 4)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid, pid_max_length, 4)
  real(kind=r_def),    intent(in)    :: basis_wx(1, ndf_wx, ndf_ww)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_ww)

  real(kind=r_tran),   intent(out)   :: weights_x_parallel(undf_ww)
  real(kind=r_tran),   intent(out)   :: weights_x_perp(undf_ww)
  real(kind=r_tran),   intent(out)   :: weights_y_parallel(undf_ww)
  real(kind=r_tran),   intent(out)   :: weights_y_perp(undf_ww)
  integer(kind=i_def), intent(out)   :: indices_x(undf_ww)
  integer(kind=i_def), intent(out)   :: indices_y(undf_ww)

  real(kind=r_def),    intent(in)    :: chi_1(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_2(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_3(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: alpha_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_y(undf_wx_2d)
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

  integer(kind=i_def) :: df, dfc, i, w2h_idx, pid_idx
  integer(kind=i_def) :: ew_sel, ns_sel, j
  integer(kind=i_def) :: ndata
  integer(kind=i_def) :: owned_panel, neighbour_panel, edge_direction
  logical(kind=l_def) :: x_direction, near_panel_edge, end_of_eave
  real(kind=r_def)    :: stretch_factor
  real(kind=r_def)    :: inverse_rot_matrix(3,3)
  logical(kind=l_def) :: to_rotate, to_stretch

  ! Lagrange donor-stencil interpolation results, for the currently active
  ! direction (x or y)
  real(kind=r_def)    :: lagrange_w(4)
  integer(kind=i_def) :: lagrange_id(4)

  ! Rotation-factor results, for the currently active direction (x or y)
  real(kind=r_def)    :: rotation_parallel, rotation_perp

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

  ! Loop only over the faces this column owns (per face_selector_ew/ns), so
  ! that each shared W2H face is written by exactly one of its two
  ! neighbouring columns
  do j = 1, ABS(ew_sel) + ABS(ns_sel)
    df = face_from_face_selector(j, ew_sel, ns_sel)

    w2h_idx = map_ww(df)

    ! ======================================================================== !
    ! Determine things depending on DoF: is DoF at the end of an eave? At such
    ! DoFs the eave value is shared between neighbouring eaves and is not used.
    select case (df)
    case (S)
      end_of_eave = panel_edge_dist_S(pid_idx) == 1
    case (N)
      end_of_eave = panel_edge_dist_N(pid_idx) == 1
    case (W)
      end_of_eave = panel_edge_dist_W(pid_idx) == 1
    case (E)
      end_of_eave = panel_edge_dist_E(pid_idx) == 1
    end select

    ! Determine which remap direction (x or y) applies to this DoF, and
    ! whether it is near the relevant type of panel edge
    select case (df)
    case (S, N)
      x_direction = .true.
      near_panel_edge = (                                                     &
          ABS(panel_edge_dist_W(pid_idx)) <= remap_depth                      &
          .or. ABS(panel_edge_dist_E(pid_idx)) <= remap_depth                 &
      )
      if (ABS(panel_edge_dist_W(pid_idx)) < ABS(panel_edge_dist_E(pid_idx))) then
        edge_direction = W
      else
        edge_direction = E
      end if
    case (W, E)
      x_direction = .false.
      near_panel_edge = (                                                     &
          ABS(panel_edge_dist_S(pid_idx)) <= remap_depth                      &
          .or. ABS(panel_edge_dist_N(pid_idx)) <= remap_depth                 &
      )
      if (ABS(panel_edge_dist_S(pid_idx)) < ABS(panel_edge_dist_N(pid_idx))) then
        edge_direction = S
      else
        edge_direction = N
      end if
    end select

    ! ======================================================================== !
    if (near_panel_edge .and. .not. end_of_eave) then

      neighbour_panel = panel_neighbour(owned_panel, edge_direction)

      call compute_lagrange_weights(df, x_direction, neighbour_panel,          &
                                     lagrange_w, lagrange_id)
      call compute_rotation_factors(df, x_direction, neighbour_panel,          &
                                     rotation_parallel, rotation_perp)

      write(log_scratch_space, *) &
        'W2H EDGE target_df=', df, ' panel=', owned_panel,    &
        ' neigh=', neighbour_panel, ' edge=', edge_direction,                  &
        ' dist_W=', panel_edge_dist_W(pid_idx), ' dist_E=', panel_edge_dist_E(pid_idx), &
        ' dist_S=', panel_edge_dist_S(pid_idx), ' dist_N=', panel_edge_dist_N(pid_idx), &
        ' rot_parallel=', rotation_parallel, ' rot_perp=', rotation_perp
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

      if (x_direction) then
        do i = 1, ndata
          weights_x_parallel(w2h_idx+i-1) =                                    &
            REAL(lagrange_w(i) * rotation_parallel, r_tran)
          weights_x_perp(w2h_idx+i-1) =                                        &
            REAL(lagrange_w(i) * rotation_perp, r_tran)
          indices_x(w2h_idx+i-1) = lagrange_id(i)
        end do
      else
        do i = 1, ndata
          weights_y_parallel(w2h_idx+i-1) =                                    &
            REAL(lagrange_w(i) * rotation_parallel, r_tran)
          weights_y_perp(w2h_idx+i-1) =                                        &
            REAL(lagrange_w(i) * rotation_perp, r_tran)
          indices_y(w2h_idx+i-1) = lagrange_id(i)
        end do
      end if

    else

      ! Not near a relevant panel edge, or at the end of an eave: default to
      ! the identity remap (own DoF, unit weight, zero cross-contribution)
      if (x_direction) then
        weights_x_parallel(w2h_idx) = 1.0_r_tran
        weights_x_perp(w2h_idx) = 0.0_r_tran
        indices_x(w2h_idx) = 1
        do i = 2, ndata
          weights_x_parallel(w2h_idx+i-1) = 0.0_r_tran
          weights_x_perp(w2h_idx+i-1) = 0.0_r_tran
          indices_x(w2h_idx+i-1) = 1
        end do
      else
        weights_y_parallel(w2h_idx) = 1.0_r_tran
        weights_y_perp(w2h_idx) = 0.0_r_tran
        indices_y(w2h_idx) = 1
        do i = 2, ndata
          weights_y_parallel(w2h_idx+i-1) = 0.0_r_tran
          weights_y_perp(w2h_idx+i-1) = 0.0_r_tran
          indices_y(w2h_idx+i-1) = 1
        end do
      end if

    end if

    ! The direction NOT active at this DoF is given the (harmless) identity
    ! default, since every W2H field must have a value at every DoF
    if (x_direction) then
      weights_y_parallel(w2h_idx) = 1.0_r_tran
      weights_y_perp(w2h_idx) = 0.0_r_tran
      indices_y(w2h_idx) = 1
      do i = 2, ndata
        weights_y_parallel(w2h_idx+i-1) = 0.0_r_tran
        weights_y_perp(w2h_idx+i-1) = 0.0_r_tran
        indices_y(w2h_idx+i-1) = 1
      end do
    else
      weights_x_parallel(w2h_idx) = 1.0_r_tran
      weights_x_perp(w2h_idx) = 0.0_r_tran
      indices_x(w2h_idx) = 1
      do i = 2, ndata
        weights_x_parallel(w2h_idx+i-1) = 0.0_r_tran
        weights_x_perp(w2h_idx+i-1) = 0.0_r_tran
        indices_x(w2h_idx+i-1) = 1
      end do
    end if

  end do  ! faces

contains

  !-----------------------------------------------------------------------------
  !> @brief Compute the Lagrange donor-stencil interpolation weights/indices
  !!        tangential to the panel edge, for the target W2H DoF "target_df".
  !!        This mirrors panel_edge_weights_kernel_mod's panel_edge_weights_1d,
  !!        but evaluated at a W2H face node instead of a W3/Wtheta cell node.
  subroutine compute_lagrange_weights(target_df, is_x_direction,               &
                                       swapped_panel, wgt, id)

    implicit none

    integer(kind=i_def), intent(in)  :: target_df
    logical(kind=l_def), intent(in)  :: is_x_direction
    integer(kind=i_def), intent(in)  :: swapped_panel
    real(kind=r_def),    intent(out) :: wgt(4)
    integer(kind=i_def), intent(out) :: id(4)

    integer(kind=i_def) :: branch_l, branch_r
    integer(kind=i_def) :: i
    integer(kind=i_def) :: wx_size_l, wx_size_r
    integer(kind=i_def) :: pid_size_l, pid_size_r
    integer(kind=i_def) :: ncells_in_stencil
    integer(kind=i_def) :: wx_stencil_1d(ndf_wx, 2*wx_max_length-1)
    integer(kind=i_def) :: pid_stencil_1d(2*pid_max_length-1)
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

    ! Direction of interpolation (alpha or beta), based on owned/swapped panels
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

    ! Interpolation point x0: the target W2H DoF's own extended-mesh coords,
    ! converted back into the owned panel's native coordinates
    abh = 0.0_r_def
    if (is_x_direction) then
      do dfc = 1, ndf_wx_2d
        abh(1) = abh(1) + alpha_x(map_wx_2d(dfc))*basis_wx_2d(1,dfc,target_df)
        abh(2) = abh(2) + beta_x(map_wx_2d(dfc))*basis_wx_2d(1,dfc,target_df)
      end do
    else
      do dfc = 1, ndf_wx_2d
        abh(1) = abh(1) + alpha_y(map_wx_2d(dfc))*basis_wx_2d(1,dfc,target_df)
        abh(2) = abh(2) + beta_y(map_wx_2d(dfc))*basis_wx_2d(1,dfc,target_df)
      end do
    end if

    call alphabetar2xyz(                                                       &
            abh(1), abh(2), unit_radius, swapped_panel, xyz(1), xyz(2), xyz(3) &
    )
    call xyz2alphabetar(                                                       &
            xyz(1), xyz(2), xyz(3), owned_panel, abh(1), abh(2), h_dummy       &
    )
    x0 = abh(interp_dir)

    ! Locations of donor points, in owned-panel native coordinates
    do i = 1, ncells_in_stencil

      panel_of_donor = int(panel_id(pid_stencil_1d(i)), i_def)

      if ( owned_panel /= panel_of_donor ) then
        x1(i) = LARGE_REAL_POSITIVE
      else
        abh = 0.0_r_def
        do dfc = 1, ndf_wx
          call chi2xyz(                                                        &
                  chi_1(wx_stencil_1d(dfc,i)),                                 &
                  chi_2(wx_stencil_1d(dfc,i)),                                 &
                  unit_radius, owned_panel, geometry, topology, coord_system,  &
                  scaled_radius, xyz(1), xyz(2), xyz(3)                       &
          )
          if (to_rotate) then
            xyz = matmul(inverse_rot_matrix, xyz)
          end if
          if (to_stretch) then
            xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
          end if
          call xyz2alphabetar(                                                 &
                  xyz(1), xyz(2), xyz(3), owned_panel,                        &
                  alpha_owned, beta_owned, h_dummy                            &
          )
          abh(1) = abh(1) + alpha_owned * basis_wx(1,dfc,target_df)
          abh(2) = abh(2) + beta_owned * basis_wx(1,dfc,target_df)
        end do
        x1(i) = abh(interp_dir)
      end if
    end do

    ! Lagrange interpolation weights and donor indices
    dx = ABS(x1 - x0)
    id1 = MINLOC(dx(1:ncells_in_stencil), 1)
    dx(id1) = LARGE_REAL_POSITIVE
    id2 = MINLOC(dx(1:ncells_in_stencil), 1)

    wgt = 0.0_r_def
    id = 1

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

      id(1) = id1
      id(2) = id2
      id(3) = id3
      id(4) = id4
    else
      wgt(1) = (x0 - x1(id2))/(x1(id1) - x1(id2))
      wgt(2) = (x0 - x1(id1))/(x1(id2) - x1(id1))

      id(1) = id1
      id(2) = id2
    end if

  end subroutine compute_lagrange_weights

  !-----------------------------------------------------------------------------
  !> @brief Compute the cos/sin rotation factors between the native and eave
  !!        coordinate frames at the target W2H DoF "target_df". This mirrors
  !!        the rotation-factor calculation of peregrin_weights_kernel_mod.
  subroutine compute_rotation_factors(target_df, is_x_direction,               &
                                       neighbour_panel_id,                     &
                                       parallel_factor, perp_factor)

    implicit none

    integer(kind=i_def), intent(in)  :: target_df
    logical(kind=l_def), intent(in)  :: is_x_direction
    integer(kind=i_def), intent(in)  :: neighbour_panel_id
    real(kind=r_def),    intent(out) :: parallel_factor, perp_factor

    integer(kind=i_def) :: dfc, idx
    real(kind=r_def)    :: chi_1_df, chi_2_df, chi_3_df
    real(kind=r_def)    :: alpha_native, beta_native
    real(kind=r_def)    :: alpha_eave, beta_eave
    real(kind=r_def)    :: native_coord_parallel, native_coord_perp
    real(kind=r_def)    :: eave_coord, direction_sign
    real(kind=r_def)    :: xyz(3), r

    chi_1_df = 0.0_r_def
    chi_2_df = 0.0_r_def
    chi_3_df = 0.0_r_def
    do dfc = 1, ndf_wx
      chi_1_df = chi_1_df + chi_1(map_wx(dfc)) * basis_wx(1, dfc, target_df)
      chi_2_df = chi_2_df + chi_2(map_wx(dfc)) * basis_wx(1, dfc, target_df)
      chi_3_df = chi_3_df + chi_3(map_wx(dfc)) * basis_wx(1, dfc, target_df)
    end do

    call chi2xyz(                                                              &
        chi_1_df, chi_2_df, chi_3_df, owned_panel,                             &
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
        xyz(1), xyz(2), xyz(3), owned_panel, alpha_native, beta_native, r      &
    )

    alpha_eave = 0.0_r_def
    beta_eave = 0.0_r_def
    if (is_x_direction) then
      do dfc = 1, ndf_wx_2d
        idx = map_wx_2d(dfc)
        alpha_eave = alpha_eave + alpha_x(idx) * basis_wx_2d(1, dfc, target_df)
        beta_eave = beta_eave + beta_x(idx) * basis_wx_2d(1, dfc, target_df)
      end do
    else
      do dfc = 1, ndf_wx_2d
        idx = map_wx_2d(dfc)
        alpha_eave = alpha_eave + alpha_y(idx) * basis_wx_2d(1, dfc, target_df)
        beta_eave = beta_eave + beta_y(idx) * basis_wx_2d(1, dfc, target_df)
      end do
    end if

    call alphabetar2xyz(                                                       &
        alpha_eave, beta_eave, 1.0_r_def, neighbour_panel_id,                  &
        xyz(1), xyz(2), xyz(3)                                                 &
    )
    call xyz2alphabetar(                                                       &
        xyz(1), xyz(2), xyz(3), owned_panel, alpha_eave, beta_eave, r          &
    )

    if (is_x_direction) then
      native_coord_parallel = beta_native
      native_coord_perp = alpha_native
      eave_coord = beta_eave
      direction_sign = -1.0_r_def
    else
      native_coord_parallel = alpha_native
      native_coord_perp = beta_native
      eave_coord = alpha_eave
      direction_sign = 1.0_r_def
    end if

    parallel_factor = COS(eave_coord) / COS(native_coord_parallel)
    perp_factor = - direction_sign * SIN(eave_coord) / COS(native_coord_perp)

    write(log_scratch_space, *) &
      'W2H ROT target_df=', target_df, ' xdir=', is_x_direction,           &
      ' neighbour=', neighbour_panel_id,                                    &
      ' alpha_native=', alpha_native,                                      &
      ' beta_native=', beta_native,                                        &
      ' alpha_eave=', alpha_eave,                                          &
      ' beta_eave=', beta_eave,                                            &
      ' parallel_factor=', parallel_factor,                                &
      ' perp_factor=', perp_factor
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  end subroutine compute_rotation_factors

end subroutine panel_edge_remap_w2h_weights_code

end module panel_edge_remap_w2h_weights_kernel_mod
