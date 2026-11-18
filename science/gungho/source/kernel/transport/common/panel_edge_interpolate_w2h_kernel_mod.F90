!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Apply the pure Lagrange donor-stencil interpolation computed by
!!        panel_edge_remap_w2h_weights_kernel_mod to a W2H flux/wind
!!        field near the edges of cubed-sphere panels, with no rotation
!!        or area scaling applied (that is done later by
!!        panel_edge_remap_w2h_kernel_mod).
!!
!!          For each selector-owned W2H face, this kernel interpolates every
!!          active panel-labelled value using signed donor offsets and the
!!          computing cell's existing CROSS2D flux stencil. The target face
!!          orientation is determined separately for each panel label, which
!!          handles rotated panel edges without changing face ownership.
!!
!!          The perpendicular value is formed from two vertex interpolations
!!          using the matching panel-labelled weights and offsets. The W3
!!          staging fields remain x/y placeholders for the later W3 path.
!!
!!          parl_p_at_w3_{x,y} is the staging output for the "perp_w2 ->
!!          W3" donor set (set 3 in the instructions doc). In this phase it
!!          is still written as a zero placeholder only, so that the field
!!          exists with stable meta_args ahead of a later phase filling it in.
module panel_edge_interpolate_w2h_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER, GH_LOGICAL,              &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_9,                    &
                                 ANY_DISCONTINUOUS_SPACE_10,                   &
                                 CELL_COLUMN,                                  &
                                 STENCIL, CROSS2D
use fs_continuity_mod,     only: W2h, W3
use constants_mod,         only: r_def, r_tran, i_def, l_def
use reference_element_mod, only: W, S, N, E
use panel_edge_support_mod, only: panel_neighbour, panel_axis_label
use sci_face_selector_support_mod, only: face_from_face_selector

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> @brief Kernel metadata for the W2H panel-edge interpolation stage.
!> @details The type declaration for the kernel. Contains the metadata
!!          needed by the Psy layer.
type, public, extends(kernel_type) :: panel_edge_interpolate_w2h_kernel_type
  private
    type(arg_type) :: meta_args(21) = (/                                         &
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, W2h),                         & ! parl_w2_at_w2
      arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, W2h),                         & ! perp_w2_at_w2
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                          & ! parl_p_at_w3_x (placeholder)
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                          & ! parl_p_at_w3_y (placeholder)
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h, STENCIL(CROSS2D)),       & ! field_for_x (area-normalised flux)
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h, STENCIL(CROSS2D)),       & ! field_for_y
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),   & ! interp_w2h_weights
      arg_type(GH_FIELD*3, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_5),   & ! interp_w2h_indices
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_10),  & ! interp_perp_weights
      arg_type(GH_FIELD*3, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_10),  & ! interp_perp_indices
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_weights_x
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_indices_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_weights_y
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_indices_y
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   & ! panel_edge_dist_W/E/S/N
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   & ! panel_id
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   & ! face_selector_ew
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   & ! face_selector_ns
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                               & ! depth
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                               & ! ndata
      arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                                & ! include_perp_p
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: panel_edge_interpolate_w2h_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_interpolate_w2h_code

contains

!> @brief Apply the "parl_w2" Lagrange interpolation to a W2H field near
!!        cubed-sphere panel edges.
!> @param[in] nlayers Number of layers
!> @param[in,out] parl_w2_at_w2_1 Interpolated parl values for panel-axis label 1
!> @param[in,out] parl_w2_at_w2_2 Interpolated parl values for panel-axis label 2
!> @param[in,out] parl_w2_at_w2_3 Interpolated parl values for panel-axis label 3
!> @param[in,out] perp_w2_at_w2_1 Interpolated perp values for panel-axis label 1
!> @param[in,out] perp_w2_at_w2_2 Interpolated perp values for panel-axis label 2
!> @param[in,out] perp_w2_at_w2_3 Interpolated perp values for panel-axis label 3
!> @param[in,out] parl_p_at_w3_x Placeholder (zero) W3 staging output (x)
!> @param[in,out] parl_p_at_w3_y Placeholder (zero) W3 staging output (y)
!> @param[in] field_for_x Field to be interpolated in the x direction
!> @param[in] wsx_stencil_size Size of field stencil for x (num cells)
!> @param[in] wsx_max_length Maximum stencil branch length for x
!> @param[in] wsx_stencil_map DoF map for the field stencil (x)
!> @param[in] field_for_y Field to be interpolated in the y direction
!> @param[in] wsy_stencil_size Size of field stencil for y (num cells)
!> @param[in] wsy_max_length Maximum stencil branch length for y
!> @param[in] wsy_stencil_map DoF map for the field stencil (y)
!> @param[in] interp_w2h_weights_1 PARL donor weights for label 1
!> @param[in] interp_w2h_weights_2 PARL donor weights for label 2
!> @param[in] interp_w2h_weights_3 PARL donor weights for label 3
!> @param[in] interp_w2h_indices_1 Signed PARL donor offsets for label 1
!> @param[in] interp_w2h_indices_2 Signed PARL donor offsets for label 2
!> @param[in] interp_w2h_indices_3 Signed PARL donor offsets for label 3
!> @param[in] interp_perp_weights_1 Perpendicular donor weights for label 1
!> @param[in] interp_perp_weights_2 Perpendicular donor weights for label 2
!> @param[in] interp_perp_weights_3 Perpendicular donor weights for label 3
!> @param[in] interp_perp_indices_1 Signed perpendicular donor offsets for label 1
!> @param[in] interp_perp_indices_2 Signed perpendicular donor offsets for label 2
!> @param[in] interp_perp_indices_3 Signed perpendicular donor offsets for label 3
!> @param[in] interp_w3_weights_x Lagrange donor weights for the W3
!!                               x remap (set 3)
!> @param[in] interp_w3_indices_x Donor cell indices for the W3 x remap
!> @param[in] interp_w3_weights_y Lagrange donor weights for the W3
!!                               y remap (set 3)
!> @param[in] interp_w3_indices_y Donor cell indices for the W3 y remap
!> @param[in] panel_edge_dist_W Distance of this column to the W panel edge
!> @param[in] panel_edge_dist_E Distance of this column to the E panel edge
!> @param[in] panel_edge_dist_S Distance of this column to the S panel edge
!> @param[in] panel_edge_dist_N Distance of this column to the N panel edge
!> @param[in] panel_id Panel ID for this column
!> @param[in] face_selector_ew East-West face selector, used to ensure
!!                            each W/E face is only ever written by one of
!!                            its two neighbouring columns
!> @param[in] face_selector_ns North-South face selector, as
!!                            face_selector_ew for S/N faces
!> @param[in] depth Maximum halo depth to consider
!> @param[in] ndata Number of donor points in the interpolation
!> @param[in] ndf_wr Num DoFs per cell for the W2H fields (the outputs and
!!                  field_for_x/y all share this space, so PSyclone generates
!!                  a single dofmap group for all of them)
!> @param[in] undf_wr Num DoFs for this partition for the W2H fields
!> @param[in] map_wr DoF map for the W2H fields
!> @param[in] ndf_w3 Num DoFs per cell for the W3 placeholder
!> @param[in] undf_w3 Num DoFs for this partition for the W3 placeholder
!> @param[in] map_w3 DoF map for the W3 placeholder
!> @param[in] ndf_w3w Num DoFs per cell for the W3 centre weights/indices
!!                  (set 3), distinct from ndf_w3 since this is a multidata
!!                  space
!> @param[in] undf_w3w Num DoFs for this partition for the W3 centre
!!                   weights/indices
!> @param[in] map_w3w DoF map for the W3 centre weights/indices
!> @param[in] ndf_ww Num DoFs per cell for interpolation weights
!> @param[in] undf_ww Num DoFs for this partition for weights
!> @param[in] map_ww DoF map for interpolation weights
!> @param[in] ndf_wp Num DoFs per cell for the perp_w2 weights/indices
!> @param[in] undf_wp Num DoFs for this partition for the perp_w2
!!                    weights/indices
!> @param[in] map_wp DoF map for the perp_w2 weights/indices
!> @param[in] ndf_pid Num DoFs per cell for panel ID field
!> @param[in] undf_pid Num DoFs for this partition for panel ID
!> @param[in] map_pid DoF map for panel ID field
subroutine panel_edge_interpolate_w2h_code( nlayers,                           &
                                            parl_w2_at_w2_1,                  &
                                            parl_w2_at_w2_2,                  &
                                            parl_w2_at_w2_3,                  &
                                            perp_w2_at_w2_1,                  &
                                            perp_w2_at_w2_2,                  &
                                            perp_w2_at_w2_3,                  &
                                            parl_p_at_w3_x,                   &
                                            parl_p_at_w3_y,                   &
                                            field_for_x,                      &
                                            wsx_stencil_size,                 &
                                            wsx_max_length,                   &
                                            wsx_stencil_map,                  &
                                            field_for_y,                      &
                                            wsy_stencil_size,                 &
                                            wsy_max_length,                   &
                                            wsy_stencil_map,                  &
                                            interp_w2h_weights_1,             &
                                            interp_w2h_weights_2,             &
                                            interp_w2h_weights_3,             &
                                            interp_w2h_indices_1,             &
                                            interp_w2h_indices_2,             &
                                            interp_w2h_indices_3,             &
                                            interp_perp_weights_1,            &
                                            interp_perp_weights_2,            &
                                            interp_perp_weights_3,            &
                                            interp_perp_indices_1,            &
                                            interp_perp_indices_2,            &
                                            interp_perp_indices_3,            &
                                            interp_w3_weights_x,              &
                                            interp_w3_indices_x,              &
                                            interp_w3_weights_y,              &
                                            interp_w3_indices_y,              &
                                            panel_edge_dist_W,                &
                                            panel_edge_dist_E,                &
                                            panel_edge_dist_S,                &
                                            panel_edge_dist_N,                &
                                            panel_id,                         &
                                            face_selector_ew,                 &
                                            face_selector_ns,                 &
                                            remap_depth,                      &
                                            ndata,                            &
                                            include_perp_p,                   &
                                            ndf_w2h, undf_w2h, map_w2h,          &
                                            ndf_w3, undf_w3, map_w3,          &
                                            ndf_ww, undf_ww, map_ww,          &
                                            ndf_wp, undf_wp, map_wp,          &
                                            ndf_w3w, undf_w3w, map_w3w,       &
                                            ndf_pid, undf_pid, map_pid        &
  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: wsx_stencil_size(4), wsy_stencil_size(4)
  integer(kind=i_def), intent(in)    :: wsx_max_length, wsy_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3w, undf_w3w
  integer(kind=i_def), intent(in)    :: ndf_wp, undf_wp
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: remap_depth
  integer(kind=i_def), intent(in)    :: ndata
  logical(kind=l_def), intent(in)    :: include_perp_p

  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w3w(ndf_w3w)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wp(ndf_wp)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wsx_stencil_map(ndf_w2h, wsx_max_length, 4)
  integer(kind=i_def), intent(in)    :: wsy_stencil_map(ndf_w2h, wsy_max_length, 4)

  real(kind=r_tran),   intent(inout) :: parl_w2_at_w2_1(undf_w2h)
  real(kind=r_tran),   intent(inout) :: parl_w2_at_w2_2(undf_w2h)
  real(kind=r_tran),   intent(inout) :: parl_w2_at_w2_3(undf_w2h)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_1(undf_w2h)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_2(undf_w2h)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_3(undf_w2h)
  real(kind=r_tran),   intent(inout) :: parl_p_at_w3_x(undf_w3)
  real(kind=r_tran),   intent(inout) :: parl_p_at_w3_y(undf_w3)
  real(kind=r_tran),   intent(in)    :: field_for_x(undf_w2h)
  real(kind=r_tran),   intent(in)    :: field_for_y(undf_w2h)
  real(kind=r_tran),   intent(in)    :: interp_w2h_weights_1(undf_ww)
  real(kind=r_tran),   intent(in)    :: interp_w2h_weights_2(undf_ww)
  real(kind=r_tran),   intent(in)    :: interp_w2h_weights_3(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_1(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_2(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_3(undf_ww)
  real(kind=r_tran),   intent(in)    :: interp_perp_weights_1(undf_wp)
  real(kind=r_tran),   intent(in)    :: interp_perp_weights_2(undf_wp)
  real(kind=r_tran),   intent(in)    :: interp_perp_weights_3(undf_wp)
  integer(kind=i_def), intent(in)    :: interp_perp_indices_1(undf_wp)
  integer(kind=i_def), intent(in)    :: interp_perp_indices_2(undf_wp)
  integer(kind=i_def), intent(in)    :: interp_perp_indices_3(undf_wp)
  real(kind=r_tran),   intent(in)    :: interp_w3_weights_x(undf_w3w)
  integer(kind=i_def), intent(in)    :: interp_w3_indices_x(undf_w3w)
  real(kind=r_tran),   intent(in)    :: interp_w3_weights_y(undf_w3w)
  integer(kind=i_def), intent(in)    :: interp_w3_indices_y(undf_w3w)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)

  ! Internal variables
  integer(kind=i_def) :: face, pid_idx, nvert
  integer(kind=i_def) :: ew_sel, ns_sel, j, label, target_face
  integer(kind=i_def) :: w3_idx, w3w_idx, w2h_idx, ww_idx
  integer(kind=i_def) :: owned_panel, other_panel, other_face, own_label
  integer(kind=i_def) :: other_label, host_sign, edge_face_dist
  integer(kind=i_def) :: edge_direction_x, edge_direction_y
  logical(kind=l_def) :: label_active, skip_perp
  logical(kind=l_def) :: near_panel_edge_x, near_panel_edge_y, near_panel_edge

  ! Number of vertical levels to loop over (W2H is layer-based, like W3)
  nvert = nlayers - 1

  pid_idx = map_pid(1)
  owned_panel = INT(panel_id(pid_idx), i_def)

  if (ABS(panel_edge_dist_W(pid_idx)) < ABS(panel_edge_dist_E(pid_idx))) then
    edge_direction_x = W
  else
    edge_direction_x = E
  end if
  if (ABS(panel_edge_dist_S(pid_idx)) < ABS(panel_edge_dist_N(pid_idx))) then
    edge_direction_y = S
  else
    edge_direction_y = N
  end if

  near_panel_edge_x = (                                                        &
      ABS(panel_edge_dist_W(pid_idx)) <= remap_depth+1                         &
      .or. ABS(panel_edge_dist_E(pid_idx)) <= remap_depth+1                    &
  )
  near_panel_edge_y = (                                                        &
      ABS(panel_edge_dist_S(pid_idx)) <= remap_depth+1                         &
      .or. ABS(panel_edge_dist_N(pid_idx)) <= remap_depth+1                    &
  )
  near_panel_edge = (near_panel_edge_x .or. near_panel_edge_y)

  ! Set 3: "parl_p donors -> W3 centre targets" interpolation (Route B of
  ! the instructions doc). There is a single W3 DoF per column; the donor
  ! sampling DoF is fixed (S for the x remap, W for the y remap - see
  ! panel_edge_remap_w2h_weights_kernel_mod). The weights default to the
  ! identity away from a panel edge, so this is unconditional
  w3_idx = map_w3(1)
  w3w_idx = map_w3w(1)

  if (near_panel_edge) then

    if (include_perp_p .and. near_panel_edge_x) then
      call panel_edge_interpolate_w3_1d(                                           &
              nlayers, parl_p_at_w3_x, w3_idx, w3w_idx, field_for_x, S,           &
              wsx_max_length,                                                     &
              wsx_stencil_size(S), wsx_stencil_map(:,:,S),                        &
              wsx_stencil_size(N), wsx_stencil_map(:,:,N),                        &
              interp_w3_weights_x, interp_w3_indices_x, ndata,                    &
              ndf_w2h, undf_w2h, undf_w3w                                           &
      )
    else
      parl_p_at_w3_x(w3_idx : w3_idx+nlayers-1) = -17000_r_tran
    end if
    if (include_perp_p .and. near_panel_edge_y) then
      call panel_edge_interpolate_w3_1d(                                           &
              nlayers, parl_p_at_w3_y, w3_idx, w3w_idx, field_for_y, W,           &
              wsy_max_length,                                                     &
              wsy_stencil_size(W), wsy_stencil_map(:,:,W),                        &
              wsy_stencil_size(E), wsy_stencil_map(:,:,E),                        &
              interp_w3_weights_y, interp_w3_indices_y, ndata,                    &
              ndf_w2h, undf_w2h, undf_w3w                                           &
      )
    else
      parl_p_at_w3_y(w3_idx : w3_idx+nlayers-1) = -17000_r_tran
    end if
  end if  ! near_panel_edge

  ew_sel = face_selector_ew(pid_idx)
  ns_sel = face_selector_ns(pid_idx)

  do j = 1, ABS(ew_sel) + ABS(ns_sel)
    face = face_from_face_selector(j, ew_sel, ns_sel)
    ww_idx = map_ww(face)
    w2h_idx = map_w2h(face)

    select case (face)
    case (W)
      edge_face_dist = panel_edge_dist_W(pid_idx)
    case (E)
      edge_face_dist = panel_edge_dist_E(pid_idx)
    case (S)
      edge_face_dist = panel_edge_dist_S(pid_idx)
    case (N)
      edge_face_dist = panel_edge_dist_N(pid_idx)
    end select

    other_panel = owned_panel
    other_face = face
    if (edge_face_dist == 1) then
      other_panel = panel_neighbour(owned_panel, face)
      other_face = face_direction_to_panel(other_panel, owned_panel)
      if (owned_panel < other_panel) then
        host_sign = 1_i_def
      else
        host_sign = -1_i_def
      end if
    else if (face == S .or. face == W) then
      host_sign = 1_i_def
    else
      host_sign = -1_i_def
    end if

    own_label = panel_axis_label(owned_panel)
    other_label = panel_axis_label(other_panel)
    parl_w2_at_w2_1(w2h_idx:w2h_idx+nlayers-1) = 0.0_r_tran
    parl_w2_at_w2_2(w2h_idx:w2h_idx+nlayers-1) = 0.0_r_tran
    parl_w2_at_w2_3(w2h_idx:w2h_idx+nlayers-1) = 0.0_r_tran
    perp_w2_at_w2_1(w2h_idx:w2h_idx+nlayers-1) = 0.0_r_tran
    perp_w2_at_w2_2(w2h_idx:w2h_idx+nlayers-1) = 0.0_r_tran
    perp_w2_at_w2_3(w2h_idx:w2h_idx+nlayers-1) = 0.0_r_tran

    do label = 1, 3
      target_face = face
      if (label /= own_label) then
        if (edge_face_dist /= 1 .or. label /= other_label) cycle
        target_face = other_face
      end if

      skip_perp = .false.
      if (include_perp_p) then
        if (edge_face_dist == 1) then
          skip_perp = .true.
        else if ((face == W .or. face == E) .and. near_panel_edge_x) then
          skip_perp = (label == panel_axis_label(                            &
              panel_neighbour(owned_panel, edge_direction_x)))
        else if ((face == S .or. face == N) .and. near_panel_edge_y) then
          skip_perp = (label == panel_axis_label(                            &
              panel_neighbour(owned_panel, edge_direction_y)))
        end if
      end if

      select case (label)
      case (1)
        label_active = ANY(interp_w2h_indices_1(                              &
          ww_idx:ww_idx+ndata-1) /= 0_i_def)
        if (label_active) then
            call interpolate_labelled_face(parl_w2_at_w2_1, perp_w2_at_w2_1,   &
              target_face, face, host_sign, skip_perp,                         &
              field_for_x, field_for_y,                                      &
              wsx_max_length, wsx_stencil_size, wsx_stencil_map,             &
              wsy_max_length, wsy_stencil_size, wsy_stencil_map,             &
              interp_w2h_weights_1, interp_w2h_indices_1,                     &
              interp_perp_weights_1, interp_perp_indices_1,                   &
              ndata, nlayers, ndf_w2h, undf_w2h, map_w2h,                    &
              ndf_ww, undf_ww, map_ww, ndf_wp, undf_wp, map_wp)
        end if
      case (2)
        label_active = ANY(interp_w2h_indices_2(                              &
          ww_idx:ww_idx+ndata-1) /= 0_i_def)
        if (label_active) then
            call interpolate_labelled_face(parl_w2_at_w2_2, perp_w2_at_w2_2,   &
              target_face, face, host_sign, skip_perp,                         &
              field_for_x, field_for_y,                                      &
              wsx_max_length, wsx_stencil_size, wsx_stencil_map,             &
              wsy_max_length, wsy_stencil_size, wsy_stencil_map,             &
              interp_w2h_weights_2, interp_w2h_indices_2,                     &
              interp_perp_weights_2, interp_perp_indices_2,                   &
              ndata, nlayers, ndf_w2h, undf_w2h, map_w2h,                    &
              ndf_ww, undf_ww, map_ww, ndf_wp, undf_wp, map_wp)
        end if
      case (3)
        label_active = ANY(interp_w2h_indices_3(                              &
          ww_idx:ww_idx+ndata-1) /= 0_i_def)
        if (label_active) then
            call interpolate_labelled_face(parl_w2_at_w2_3, perp_w2_at_w2_3,   &
              target_face, face, host_sign, skip_perp,                         &
              field_for_x, field_for_y,                                      &
              wsx_max_length, wsx_stencil_size, wsx_stencil_map,             &
              wsy_max_length, wsy_stencil_size, wsy_stencil_map,             &
              interp_w2h_weights_3, interp_w2h_indices_3,                     &
              interp_perp_weights_3, interp_perp_indices_3,                   &
              ndata, nlayers, ndf_w2h, undf_w2h, map_w2h,                    &
              ndf_ww, undf_ww, map_ww, ndf_wp, undf_wp, map_wp)
        end if
      end select

    end do
  end do

end subroutine panel_edge_interpolate_w2h_code

integer(kind=i_def) function face_direction_to_panel(panel, target_panel)      &
    result(direction)

  implicit none

  integer(kind=i_def), intent(in) :: panel, target_panel
  integer(kind=i_def) :: directions(4), i

  directions = (/ W, S, E, N /)
  direction = 0_i_def
  do i = 1, 4
    if (panel_neighbour(panel, directions(i)) == target_panel) then
      direction = directions(i)
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

subroutine interpolate_labelled_face(parl_output, perp_output, target_face,   &
                                     host_face, host_sign, skip_perp,         &
                                     field_for_x, field_for_y,               &
                                     wsx_max_length, wsx_stencil_size,       &
                                     wsx_stencil_map, wsy_max_length,        &
                                     wsy_stencil_size, wsy_stencil_map,      &
                                     weights, offsets, perp_weights,         &
                                     perp_offsets, ndata, nlayers,            &
                                     ndf_w2h, undf_w2h, map_w2h,             &
                                     ndf_ww, undf_ww, map_ww,                &
                                     ndf_wp, undf_wp, map_wp)

  implicit none

  integer(kind=i_def), intent(in) :: target_face, host_face, host_sign
  logical(kind=l_def), intent(in) :: skip_perp
  integer(kind=i_def), intent(in) :: wsx_max_length, wsy_max_length
  integer(kind=i_def), intent(in) :: wsx_stencil_size(4), wsy_stencil_size(4)
  integer(kind=i_def), intent(in) :: ndata, nlayers
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_ww, undf_ww, ndf_wp, undf_wp
  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in) :: map_ww(ndf_ww), map_wp(ndf_wp)
  integer(kind=i_def), intent(in) :: wsx_stencil_map(                        &
      ndf_w2h, wsx_max_length, 4)
  integer(kind=i_def), intent(in) :: wsy_stencil_map(                        &
      ndf_w2h, wsy_max_length, 4)
  real(kind=r_tran), intent(inout) :: parl_output(undf_w2h)
  real(kind=r_tran), intent(inout) :: perp_output(undf_w2h)
  real(kind=r_tran), intent(in) :: field_for_x(undf_w2h)
  real(kind=r_tran), intent(in) :: field_for_y(undf_w2h)
  real(kind=r_tran), intent(in) :: weights(undf_ww)
  integer(kind=i_def), intent(in) :: offsets(undf_ww)
  real(kind=r_tran), intent(in) :: perp_weights(undf_wp)
  integer(kind=i_def), intent(in) :: perp_offsets(undf_wp)

  if (target_face == S .or. target_face == N) then
    call interpolate_labelled_face_1d(                                      &
        parl_output, perp_output, field_for_x, target_face, host_face,        &
      host_sign, skip_perp,                                                &
        wsx_max_length, wsx_stencil_size,                                     &
        wsx_stencil_map, weights, offsets, perp_weights, perp_offsets,        &
        ndata, nlayers, ndf_w2h, undf_w2h, map_w2h,                          &
        ndf_ww, undf_ww, map_ww, ndf_wp, undf_wp, map_wp)
  else
    call interpolate_labelled_face_1d(                                      &
        parl_output, perp_output, field_for_y, target_face, host_face,        &
      host_sign, skip_perp,                                                &
        wsy_max_length, wsy_stencil_size,                                     &
        wsy_stencil_map, weights, offsets, perp_weights, perp_offsets,        &
        ndata, nlayers, ndf_w2h, undf_w2h, map_w2h,                          &
        ndf_ww, undf_ww, map_ww, ndf_wp, undf_wp, map_wp)
  end if

end subroutine interpolate_labelled_face

subroutine interpolate_labelled_face_1d(parl_output, perp_output, field,     &
                                         target_face, host_face, host_sign,   &
                                         skip_perp,                           &
                                         ws_max_length,                      &
                                         ws_stencil_size, ws_stencil_map,     &
                                         weights, offsets, perp_weights,     &
                                         perp_offsets, ndata, nlayers,        &
                                         ndf_w2h, undf_w2h, map_w2h,          &
                                         ndf_ww, undf_ww, map_ww,             &
                                         ndf_wp, undf_wp, map_wp)

  implicit none

  integer(kind=i_def), intent(in) :: target_face, host_face, host_sign
  logical(kind=l_def), intent(in) :: skip_perp
  integer(kind=i_def), intent(in) :: ws_max_length, ws_stencil_size(4)
  integer(kind=i_def), intent(in) :: ndata, nlayers
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_ww, undf_ww, ndf_wp, undf_wp
  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in) :: map_ww(ndf_ww), map_wp(ndf_wp)
  integer(kind=i_def), intent(in) :: ws_stencil_map(                        &
      ndf_w2h, ws_max_length, 4)
  real(kind=r_tran), intent(inout) :: parl_output(undf_w2h)
  real(kind=r_tran), intent(inout) :: perp_output(undf_w2h)
  real(kind=r_tran), intent(in) :: field(undf_w2h)
  real(kind=r_tran), intent(in) :: weights(undf_ww)
  integer(kind=i_def), intent(in) :: offsets(undf_ww)
  real(kind=r_tran), intent(in) :: perp_weights(undf_wp)
  integer(kind=i_def), intent(in) :: perp_offsets(undf_wp)

  integer(kind=i_def) :: i, k, source_face, source_face_a, source_face_b
  integer(kind=i_def) :: w_idx, wp_idx, output_idx, donor_idx, offset
  integer(kind=i_def) :: nvert

  nvert = nlayers - 1
  w_idx = map_ww(host_face)
  wp_idx = map_wp(host_face)
  output_idx = map_w2h(host_face)
  parl_output(output_idx:output_idx+nvert) = 0.0_r_tran
  perp_output(output_idx:output_idx+nvert) = 0.0_r_tran

  do i = 0, ndata - 1
    offset = offsets(w_idx+i)
    if (offset == 0_i_def) cycle
    donor_idx = donor_dof(target_face, offset)
    if (donor_idx == 0_i_def) cycle
    parl_output(output_idx:output_idx+nvert) =                                &
        parl_output(output_idx:output_idx+nvert)                              &
        + weights(w_idx+i) * field(donor_idx:donor_idx+nvert)
  end do

        if (skip_perp) return

  if (target_face == S .or. target_face == N) then
    source_face_a = W
    source_face_b = E
  else
    source_face_a = S
    source_face_b = N
  end if

  do k = 0, 1
    if (k == 0) then
      source_face = source_face_a
    else
      source_face = source_face_b
    end if
    do i = 0, ndata - 1
      offset = perp_offsets(wp_idx+k*ndata+i)
      if (offset == 0_i_def) cycle
      donor_idx = donor_dof(source_face, offset)
      if (donor_idx == 0_i_def) cycle
      perp_output(output_idx:output_idx+nvert) =                              &
          perp_output(output_idx:output_idx+nvert)                            &
          + 0.5_r_tran * perp_weights(wp_idx+k*ndata+i)                       &
          * field(donor_idx:donor_idx+nvert)
    end do
  end do

contains

  integer(kind=i_def) function donor_dof(source_df, donor_offset) result(dof)

    implicit none

    integer(kind=i_def), intent(in) :: source_df, donor_offset
    integer(kind=i_def) :: branch, position

    if (donor_offset == host_sign) then
      dof = map_w2h(source_df)
      return
    else if (donor_offset*host_sign < 0_i_def) then
      branch = host_face
      position = ABS(donor_offset) + 1
    else
      branch = opposite_face(host_face)
      position = ABS(donor_offset)
    end if

    if (position > ws_stencil_size(branch)) then
      dof = 0_i_def
    else
      dof = ws_stencil_map(source_df, position, branch)
    end if

  end function donor_dof

end subroutine interpolate_labelled_face_1d

!> @brief Private routine to interpolate the single W3 cell-centre target
!!        from a donor-cell stencil tangential to the relevant panel edge,
!!        using the "parl_p donors -> W3 centre targets" Lagrange weights/
!!        indices (set 3, Route B of the instructions doc). The donor
!!        sampling DoF ("sample_df") is fixed: S for the x remap, W for the
!!        y remap (see panel_edge_remap_w2h_weights_kernel_mod).
subroutine panel_edge_interpolate_w3_1d( nlayers,                              &
                                         interpolated_field, w3_idx, w3w_idx, &
                                         field, sample_df,                    &
                                         ws_max_length,                       &
                                         ws_stencil_size_l, ws_stencil_l,      &
                                         ws_stencil_size_r, ws_stencil_r,      &
                                         weights, indices, ndata,             &
                                         ndf_w2h, undf_w2h, undf_w3w            &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: w3_idx, w3w_idx
  integer(kind=i_def), intent(in)    :: sample_df
  integer(kind=i_def), intent(in)    :: ws_stencil_size_l
  integer(kind=i_def), intent(in)    :: ws_stencil_size_r
  integer(kind=i_def), intent(in)    :: ws_max_length
  integer(kind=i_def), intent(in)    :: ndata
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h, undf_w3w

  integer(kind=i_def), intent(in)    :: ws_stencil_l(ndf_w2h, ws_stencil_size_l)
  integer(kind=i_def), intent(in)    :: ws_stencil_r(ndf_w2h, ws_stencil_size_r)

  real(kind=r_tran),   intent(inout) :: interpolated_field(:)
  real(kind=r_tran),   intent(in)    :: field(undf_w2h)
  real(kind=r_tran),   intent(in)    :: weights(undf_w3w)
  integer(kind=i_def), intent(in)    :: indices(undf_w3w)

  ! Internal variables
  integer(kind=i_def) :: i
  integer(kind=i_def) :: f_idx
  integer(kind=i_def) :: nvert
  integer(kind=i_def) :: ncells_in_stencil
  integer(kind=i_def) :: stencil_1d(2*ws_max_length-1)

  nvert = nlayers - 1

  ! -------------------------------------------------------------------------- !
  ! Create a 1D donor stencil (unify the two branches of the CROSS2D stencil),
  ! sampling the donor's own face orientation ("sample_df"), from the donor
  ! cells
  ! -------------------------------------------------------------------------- !
  ncells_in_stencil = ws_stencil_size_l + ws_stencil_size_r - 1

  stencil_1d(:) = 0
  do i = 1, ws_stencil_size_l
    stencil_1d(i) = ws_stencil_l(sample_df, i)
  end do
  ! Omit 1 from the second part, as the central cell is already included
  do i = 1, ws_stencil_size_r - 1
    stencil_1d(i+ws_stencil_size_l) = ws_stencil_r(sample_df, i+1)
  end do

  ! -------------------------------------------------------------------------- !
  ! Interpolate
  ! -------------------------------------------------------------------------- !
  interpolated_field(w3_idx : w3_idx+nvert) = 0.0_r_tran

  ! Loop through (multidata) interpolation points
  do i = 0, ndata - 1
    f_idx = stencil_1d(indices(w3w_idx+i))

    interpolated_field(w3_idx : w3_idx+nvert) = interpolated_field(w3_idx : w3_idx+nvert) &
      + weights(w3w_idx+i) * field(f_idx : f_idx+nvert)
  end do

end subroutine panel_edge_interpolate_w3_1d

end module panel_edge_interpolate_w2h_kernel_mod
