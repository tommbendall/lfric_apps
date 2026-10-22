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
!> @details This is the W2H analogue of panel_edge_remap_kernel_mod (the
!!          scalar W3/Wtheta remap), restricted to the "parl_w2" donor set
!!          (set 1 in the instructions doc): donor values are taken from the
!!          same face orientation as the target DoF, and interpolated
!!          tangential to the panel edge using
!!          interp_w2h_weights_{x,y} / interp_w2h_indices_{x,y}.
!!
!!          Two independent interpolations are produced, mirroring the x/y
!!          split used throughout this kernel family:
!!            - parl_w2_at_w2_x: interpolated S/N (v-type) face values, valid
!!              near a W or E panel edge. The W/E (u-type) face values of
!!              this output are simply copied from the input field.
!!            - parl_w2_at_w2_y: interpolated W/E (u-type) face values, valid
!!              near a S or N panel edge. The S/N face values of this output
!!              are simply copied from the input field.
!!          Away from the relevant edge, the weights already default to the
!!          identity (see panel_edge_remap_w2h_weights_kernel_mod), so both
!!          outputs equal the input field everywhere they are not actively
!!          interpolating a face value.
!!
!!          perp_w2_at_w2_{x,y} are the "perp_w2" outputs (set 2 in the
!!          instructions doc): at each parl_p face (S/N for x, W/E for y) the
!!          native perpendicular donor DoFs (W and E faces for x, from
!!          field_for_x; S and N faces for y, from field_for_y) are
!!          interpolated tangential to the panel edge to each of the face's
!!          two vertices, using interp_perp_weights_{x,y} /
!!          interp_perp_indices_{x,y}, and the two results are averaged.
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
  type(arg_type) :: meta_args(26) = (/                                         &
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2h),                         & ! parl_w2_at_w2_x
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2h),                         & ! parl_w2_at_w2_y
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2h),                         & ! perp_w2_at_w2_x
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2h),                         & ! perp_w2_at_w2_y
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                          & ! parl_p_at_w3_x (placeholder)
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                          & ! parl_p_at_w3_y (placeholder)
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h, STENCIL(CROSS2D)),       & ! field_for_x (area-normalised flux)
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h, STENCIL(CROSS2D)),       & ! field_for_y
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),   & ! interp_w2h_weights_x
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_5),   & ! interp_w2h_indices_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),   & ! interp_w2h_weights_y
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_5),   & ! interp_w2h_indices_y
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_10),  & ! interp_perp_weights_x
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_10),  & ! interp_perp_indices_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_10),  & ! interp_perp_weights_y
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_10),  & ! interp_perp_indices_y
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_weights_x
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_indices_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_weights_y
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_9),   & ! interp_w3_indices_y
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   & ! panel_edge_dist_W/E/S/N
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
!> @param[in,out] parl_w2_at_w2_x Interpolated S/N values for the x
!!                                remap
!> @param[in,out] parl_w2_at_w2_y Interpolated W/E values for the y
!!                                remap
!> @param[in,out] perp_w2_at_w2_x Interpolated perp_w2 values at the S/N
!!                                faces for the x remap
!> @param[in,out] perp_w2_at_w2_y Interpolated perp_w2 values at the W/E
!!                                faces for the y remap
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
!> @param[in] interp_w2h_weights_x Lagrange donor weights for the x remap
!> @param[in] interp_w2h_indices_x Donor cell indices for the x remap
!> @param[in] interp_w2h_weights_y Lagrange donor weights for the y remap
!> @param[in] interp_w2h_indices_y Donor cell indices for the y remap
!> @param[in] interp_perp_weights_x Lagrange donor weights for the perp_w2
!!                                x remap (two sets of ndata per face)
!> @param[in] interp_perp_indices_x Donor cell indices for the perp_w2 x remap
!> @param[in] interp_perp_weights_y Lagrange donor weights for the perp_w2
!!                                y remap (two sets of ndata per face)
!> @param[in] interp_perp_indices_y Donor cell indices for the perp_w2 y remap
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
                                            parl_w2_at_w2_x,                  &
                                            parl_w2_at_w2_y,                  &
                                            perp_w2_at_w2_x,                  &
                                            perp_w2_at_w2_y,                  &
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
                                            interp_w2h_weights_x,             &
                                            interp_w2h_indices_x,             &
                                            interp_w2h_weights_y,             &
                                            interp_w2h_indices_y,             &
                                            interp_perp_weights_x,            &
                                            interp_perp_indices_x,            &
                                            interp_perp_weights_y,            &
                                            interp_perp_indices_y,            &
                                            interp_w3_weights_x,              &
                                            interp_w3_indices_x,              &
                                            interp_w3_weights_y,              &
                                            interp_w3_indices_y,              &
                                            panel_edge_dist_W,                &
                                            panel_edge_dist_E,                &
                                            panel_edge_dist_S,                &
                                            panel_edge_dist_N,                &
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

  real(kind=r_tran),   intent(inout) :: parl_w2_at_w2_x(undf_w2h)
  real(kind=r_tran),   intent(inout) :: parl_w2_at_w2_y(undf_w2h)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_x(undf_w2h)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_y(undf_w2h)
  real(kind=r_tran),   intent(inout) :: parl_p_at_w3_x(undf_w3)
  real(kind=r_tran),   intent(inout) :: parl_p_at_w3_y(undf_w3)
  real(kind=r_tran),   intent(in)    :: field_for_x(undf_w2h)
  real(kind=r_tran),   intent(in)    :: field_for_y(undf_w2h)
  real(kind=r_tran),   intent(in)    :: interp_w2h_weights_x(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_x(undf_ww)
  real(kind=r_tran),   intent(in)    :: interp_w2h_weights_y(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_y(undf_ww)
  real(kind=r_tran),   intent(in)    :: interp_perp_weights_x(undf_wp)
  integer(kind=i_def), intent(in)    :: interp_perp_indices_x(undf_wp)
  real(kind=r_tran),   intent(in)    :: interp_perp_weights_y(undf_wp)
  integer(kind=i_def), intent(in)    :: interp_perp_indices_y(undf_wp)
  real(kind=r_tran),   intent(in)    :: interp_w3_weights_x(undf_w3w)
  integer(kind=i_def), intent(in)    :: interp_w3_indices_x(undf_w3w)
  real(kind=r_tran),   intent(in)    :: interp_w3_weights_y(undf_w3w)
  integer(kind=i_def), intent(in)    :: interp_w3_indices_y(undf_w3w)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)

  ! Internal variables
  integer(kind=i_def) :: face, pid_idx, nvert
  integer(kind=i_def) :: ew_sel, ns_sel, j
  integer(kind=i_def) :: w3_idx, w3w_idx, w2h_idx
  logical(kind=l_def) :: near_panel_edge_x, near_panel_edge_y, near_panel_edge

  ! Number of vertical levels to loop over (W2H is layer-based, like W3)
  nvert = nlayers - 1

  pid_idx = map_pid(1)

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

    ew_sel = face_selector_ew(pid_idx)
    ns_sel = face_selector_ns(pid_idx)

    ! Loop only over the faces this column owns (per face_selector_ew/ns), so
    ! that each shared W2H face is interpolated by exactly one of its two
    ! neighbouring columns (the stored donor indices are only meaningful to
    ! the column that computed them)
    do j = 1, ABS(ew_sel) + ABS(ns_sel)
      face = face_from_face_selector(j, ew_sel, ns_sel)
      w2h_idx = map_w2h(face)
      select case (face)
      case (S, N)
        if (near_panel_edge_x .or. include_perp_p) then
          call panel_edge_interpolate_w2h_1d(                                    &
                  nlayers, parl_w2_at_w2_x, field_for_x, face,                  &
                  wsx_max_length,                                                &
                  wsx_stencil_size(S), wsx_stencil_map(:,:,S),                   &
                  wsx_stencil_size(N), wsx_stencil_map(:,:,N),                   &
                  interp_w2h_weights_x, interp_w2h_indices_x, ndata,            &
                  ndf_w2h, undf_w2h, map_w2h,                                       &
                  ndf_ww, undf_ww, map_ww                                        &
          )
        else
          parl_w2_at_w2_x(w2h_idx : w2h_idx+nlayers-1) = -8100_r_tran
        end if

        if (near_panel_edge_x) then
          call panel_edge_interpolate_perp_1d(                                 &
                  nlayers, perp_w2_at_w2_x, field_for_x, face, W, E,           &
                  wsx_max_length,                                              &
                  wsx_stencil_size(S), wsx_stencil_map(:,:,S),                 &
                  wsx_stencil_size(N), wsx_stencil_map(:,:,N),                 &
                  interp_perp_weights_x, interp_perp_indices_x, ndata,         &
                  ndf_w2h, undf_w2h, map_w2h,                                  &
                  ndf_wp, undf_wp, map_wp                                      &
          )
        else
          perp_w2_at_w2_x(w2h_idx : w2h_idx+nlayers-1) = 0.0_r_tran
        end if

      case (W, E)
        if (near_panel_edge_y .or. include_perp_p) then
          call panel_edge_interpolate_w2h_1d(                                    &
                  nlayers, parl_w2_at_w2_y, field_for_y, face,                  &
                  wsy_max_length,                                                &
                  wsy_stencil_size(W), wsy_stencil_map(:,:,W),                   &
                  wsy_stencil_size(E), wsy_stencil_map(:,:,E),                   &
                  interp_w2h_weights_y, interp_w2h_indices_y, ndata,            &
                  ndf_w2h, undf_w2h, map_w2h,                                       &
                  ndf_ww, undf_ww, map_ww                                        &
          )
        else
          parl_w2_at_w2_y(w2h_idx : w2h_idx+nlayers-1) = -8100_r_tran
        end if

        if (near_panel_edge_y) then
          call panel_edge_interpolate_perp_1d(                                 &
                  nlayers, perp_w2_at_w2_y, field_for_y, face, S, N,           &
                  wsy_max_length,                                              &
                  wsy_stencil_size(W), wsy_stencil_map(:,:,W),                 &
                  wsy_stencil_size(E), wsy_stencil_map(:,:,E),                 &
                  interp_perp_weights_y, interp_perp_indices_y, ndata,         &
                  ndf_w2h, undf_w2h, map_w2h,                                  &
                  ndf_wp, undf_wp, map_wp                                      &
          )
        else
          perp_w2_at_w2_y(w2h_idx : w2h_idx+nlayers-1) = 0.0_r_tran
        end if
      end select
    end do
  end if  ! near panel edge

end subroutine panel_edge_interpolate_w2h_code

!> @brief Private routine to interpolate a single W2H target DoF ("target_df")
!!        from a donor-cell stencil tangential to the relevant panel edge,
!!        using the "parl_w2" Lagrange weights/indices (same-orientation
!!        donors only, no rotation).
subroutine panel_edge_interpolate_w2h_1d( nlayers,                             &
                                          interpolated_field, field,          &
                                          target_df,                          &
                                          ws_max_length,                      &
                                          ws_stencil_size_l, ws_stencil_l,     &
                                          ws_stencil_size_r, ws_stencil_r,     &
                                          weights, indices, ndata,            &
                                          ndf_w2h, undf_w2h, map_w2h,            &
                                          ndf_ww, undf_ww, map_ww             &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: target_df
  integer(kind=i_def), intent(in)    :: ws_stencil_size_l
  integer(kind=i_def), intent(in)    :: ws_stencil_size_r
  integer(kind=i_def), intent(in)    :: ws_max_length
  integer(kind=i_def), intent(in)    :: ndata
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h

  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: ws_stencil_l(ndf_w2h, ws_stencil_size_l)
  integer(kind=i_def), intent(in)    :: ws_stencil_r(ndf_w2h, ws_stencil_size_r)

  real(kind=r_tran),   intent(inout) :: interpolated_field(undf_w2h)
  real(kind=r_tran),   intent(in)    :: field(undf_w2h)
  real(kind=r_tran),   intent(in)    :: weights(undf_ww)
  integer(kind=i_def), intent(in)    :: indices(undf_ww)

  ! Internal variables
  integer(kind=i_def) :: i
  integer(kind=i_def) :: w_idx, r_idx, f_idx
  integer(kind=i_def) :: nvert
  integer(kind=i_def) :: ncells_in_stencil
  integer(kind=i_def) :: stencil_1d(2*ws_max_length-1)

  nvert = nlayers - 1

  ! -------------------------------------------------------------------------- !
  ! Create a 1D donor stencil (unify the two branches of the CROSS2D stencil),
  ! sampling the target's own face orientation, from the donor cells
  ! -------------------------------------------------------------------------- !
  ncells_in_stencil = ws_stencil_size_l + ws_stencil_size_r - 1

  stencil_1d(:) = 0
  do i = 1, ws_stencil_size_l
    stencil_1d(i) = ws_stencil_l(target_df, i)
  end do
  ! Omit 1 from the second part, as the central cell is already included
  do i = 1, ws_stencil_size_r - 1
    stencil_1d(i+ws_stencil_size_l) = ws_stencil_r(target_df, i+1)
  end do

  ! -------------------------------------------------------------------------- !
  ! Interpolate
  ! -------------------------------------------------------------------------- !
  w_idx = map_ww(target_df)
  r_idx = map_w2h(target_df)

  interpolated_field(r_idx : r_idx+nvert) = 0.0_r_tran

  ! Loop through (multidata) interpolation points
  do i = 0, ndata - 1
    f_idx = stencil_1d(indices(w_idx+i))

    interpolated_field(r_idx : r_idx+nvert) = interpolated_field(r_idx : r_idx+nvert) &
      + weights(w_idx+i) * field(f_idx : f_idx+nvert)
  end do

end subroutine panel_edge_interpolate_w2h_1d

!> @brief Private routine to interpolate the perp_w2 value at a single W2H
!!        target DoF (target_df), tangential to the relevant panel edge. Two
!!        interpolations are made, from the donor DoFs source_df_a and
!!        source_df_b to the two vertices of the target face, using the two
!!        consecutive sets of ndata weights/indices stored for the target,
!!        and the results are averaged.
subroutine panel_edge_interpolate_perp_1d( nlayers,                            &
                                           interpolated_field, field,         &
                                           target_df,                         &
                                           source_df_a, source_df_b,          &
                                           ws_max_length,                     &
                                           ws_stencil_size_l, ws_stencil_l,    &
                                           ws_stencil_size_r, ws_stencil_r,    &
                                           weights, indices, ndata,           &
                                           ndf_w2h, undf_w2h, map_w2h,        &
                                           ndf_wp, undf_wp, map_wp            &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: target_df
  integer(kind=i_def), intent(in)    :: source_df_a, source_df_b
  integer(kind=i_def), intent(in)    :: ws_stencil_size_l
  integer(kind=i_def), intent(in)    :: ws_stencil_size_r
  integer(kind=i_def), intent(in)    :: ws_max_length
  integer(kind=i_def), intent(in)    :: ndata
  integer(kind=i_def), intent(in)    :: ndf_wp, undf_wp
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h

  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_wp(ndf_wp)
  integer(kind=i_def), intent(in)    :: ws_stencil_l(ndf_w2h, ws_stencil_size_l)
  integer(kind=i_def), intent(in)    :: ws_stencil_r(ndf_w2h, ws_stencil_size_r)

  real(kind=r_tran),   intent(inout) :: interpolated_field(undf_w2h)
  real(kind=r_tran),   intent(in)    :: field(undf_w2h)
  real(kind=r_tran),   intent(in)    :: weights(undf_wp)
  integer(kind=i_def), intent(in)    :: indices(undf_wp)

  ! Internal variables
  integer(kind=i_def) :: i, k
  integer(kind=i_def) :: w_idx, r_idx, f_idx
  integer(kind=i_def) :: nvert
  integer(kind=i_def) :: source_df
  integer(kind=i_def) :: stencil_1d(2*ws_max_length-1)

  nvert = nlayers - 1

  w_idx = map_wp(target_df)
  r_idx = map_w2h(target_df)

  interpolated_field(r_idx : r_idx+nvert) = 0.0_r_tran

  ! k = 0: vertex a (donor DoFs source_df_a); k = 1: vertex b (source_df_b)
  do k = 0, 1
    if (k == 0) then
      source_df = source_df_a
    else
      source_df = source_df_b
    end if

    ! Create a 1D donor stencil (unify the two branches of the CROSS2D
    ! stencil), sampling the donor's own face orientation
    stencil_1d(:) = 0
    do i = 1, ws_stencil_size_l
      stencil_1d(i) = ws_stencil_l(source_df, i)
    end do
    ! Omit 1 from the second part, as the central cell is already included
    do i = 1, ws_stencil_size_r - 1
      stencil_1d(i+ws_stencil_size_l) = ws_stencil_r(source_df, i+1)
    end do

    ! Loop through (multidata) interpolation points, averaging the two vertices
    do i = 0, ndata - 1
      f_idx = stencil_1d(indices(w_idx+k*ndata+i))

      interpolated_field(r_idx : r_idx+nvert) =                                &
        interpolated_field(r_idx : r_idx+nvert)                                &
        + 0.5_r_tran * weights(w_idx+k*ndata+i) * field(f_idx : f_idx+nvert)
    end do
  end do

end subroutine panel_edge_interpolate_perp_1d

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
