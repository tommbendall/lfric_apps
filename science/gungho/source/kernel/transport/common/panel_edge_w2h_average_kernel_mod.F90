!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Average the W1V vertex staging values onto the "parl_p" W2H faces
!!        near cubed-sphere panel edges, giving the "perp_w2" staged value
!!        consumed by panel_edge_remap_w2h_kernel_mod.
!!
!> @details This implements both routes of the instructions doc (section 2.4)
!!          for staging the "perp_w2" value onto a cell's final W2H faces:
!!          - Route A: perp_p_at_w1v_x/y already hold the "perp_w2" donor
!!            values interpolated onto the W1V vertices by
!!            panel_edge_interpolate_w2h_kernel_mod. For the S face of this
!!            column, average the SW and SE vertices of perp_p_at_w1v_x; for
!!            the N face, average NW and NE. Likewise for the W face, average
!!            the SW and NW vertices of perp_p_at_w1v_y; for the E face,
!!            average SE and NE. These are the "parl_p" targets for each
!!            direction's remap.
!!          - Route B: parl_p_at_w3_x/y hold the "parl_p" value interpolated
!!            onto the W3 cell centres. For the W/E faces (the "perp_p"
!!            targets of the x remap), each column adds half of its own
!!            parl_p_at_w3_x centre value to both its W and E faces; the
!!            neighbouring column across each face independently adds its own
!!            half when it is visited (no stencil read of the neighbour is
!!            needed). Likewise for the S/N faces (the "perp_p" targets of
!!            the y remap) using parl_p_at_w3_y.
!!
!!          Only near a panel edge does either contribution differ from the
!!          identity, so each direction is gated behind its own
!!          `panel_edge_dist`/`remap_depth` test (mirroring the weights
!!          kernel's `near_panel_edge_x/y`): perp_w2_at_w2_x (S/N faces via
!!          Route A, W/E faces via Route B) is gated on proximity to a W/E
!!          panel edge, perp_w2_at_w2_y (W/E faces via Route A, S/N faces via
!!          Route B) on proximity to a S/N panel edge.
!!
!!          perp_w2_at_w2_{x,y} use GH_INC rather than GH_WRITE, both so that
!!          Route B's contribution can be added after Route A's, and so that
!!          each face's two neighbouring columns can each add their own Route
!!          B half independently.
module panel_edge_w2h_average_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                   &
                                 GH_FIELD, GH_SCALAR,                        &
                                 GH_REAL, GH_INTEGER, GH_LOGICAL,             &
                                 GH_READ, GH_INC,                            &
                                 ANY_SPACE_1,                                &
                                 ANY_DISCONTINUOUS_SPACE_3,                  &
                                 CELL_COLUMN
use fs_continuity_mod,     only: W2h, W3
use constants_mod,         only: r_tran, i_def, l_def
use reference_element_mod, only: W, S, E, N, SW, SE, NW, NE
use sci_face_selector_support_mod, only: face_from_face_selector

! TODO: to remove
use log_mod, only: log_scratch_space, log_event, LOG_LEVEL_INFO

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!! Psy layer.
type, public, extends(kernel_type) :: panel_edge_w2h_average_kernel_type
  private
  type(arg_type) :: meta_args(11) = (/                                         &
      arg_type(GH_FIELD,   GH_REAL,    GH_INC,   W2h),                         & ! perp_w2_at_w2_x
      arg_type(GH_FIELD,   GH_REAL,    GH_INC,   W2h),                         & ! perp_w2_at_w2_y
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_SPACE_1),                 & ! perp_p_at_w1v_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_SPACE_1),                 & ! perp_p_at_w1v_y
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                          & ! parl_p_at_w3_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                          & ! parl_p_at_w3_y
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,                                &
               ANY_DISCONTINUOUS_SPACE_3),                                     & ! panel_edge_dist_W/E/S/N
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,                                &
               ANY_DISCONTINUOUS_SPACE_3),                                     & ! face_selector_ew
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,                                &
               ANY_DISCONTINUOUS_SPACE_3),                                     & ! face_selector_ns
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                               & ! remap_depth
      arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                                & ! include_perp_p
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: panel_edge_w2h_average_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_w2h_average_code

contains

!> @brief Average the W1V vertex values onto the "parl_p" W2H faces.
!> @param[in]     nlayers           Number of layers
!> @param[in,out] perp_w2_at_w2_x   "perp_w2" staged value for the x remap,
!!                                  incremented at S/N faces
!> @param[in,out] perp_w2_at_w2_y   "perp_w2" staged value for the y remap,
!!                                  incremented at W/E faces
!> @param[in]     perp_p_at_w1v_x   "perp_w2" value interpolated onto the W1V
!!                                  vertices for the x remap
!> @param[in]     perp_p_at_w1v_y   "perp_w2" value interpolated onto the W1V
!!                                  vertices for the y remap
!> @param[in]     parl_p_at_w3_x    "parl_p" value interpolated onto the W3
!!                                  cell centre for the x remap
!> @param[in]     parl_p_at_w3_y    "parl_p" value interpolated onto the W3
!!                                  cell centre for the y remap
!> @param[in]     panel_edge_dist_W Distance of this column to the W panel edge
!> @param[in]     panel_edge_dist_E Distance of this column to the E panel edge
!> @param[in]     panel_edge_dist_S Distance of this column to the S panel edge
!> @param[in]     panel_edge_dist_N Distance of this column to the N panel edge
!> @param[in]     face_selector_ew  East-West face selector, used to ensure
!!                                  each W/E face is only ever written by one
!!                                  of its two neighbouring columns
!> @param[in]     face_selector_ns  North-South face selector, as
!!                                  face_selector_ew for S/N faces
!> @param[in]     remap_depth       Depth (in columns) to which remapping
!!                                  applies
!> @param[in]     ndf_wr            Num DoFs per cell for the W2H field
!> @param[in]     undf_wr           Num DoFs for this partition for the W2H
!!                                  field
!> @param[in]     map_wr            DoF map for the W2H field
!> @param[in]     ndf_wv            Num DoFs per cell for the W1V field
!> @param[in]     undf_wv           Num DoFs for this partition for the W1V
!!                                  field
!> @param[in]     map_wv            DoF map for the W1V field
!> @param[in]     ndf_pid           Num DoFs per cell for panel ID field
!> @param[in]     undf_pid          Num DoFs for this partition for panel ID
!> @param[in]     map_pid           DoF map for panel ID field
!> @param[in]     ndf_w3            Num DoFs per cell for the W3 field
!> @param[in]     undf_w3           Num DoFs for this partition for the W3
!!                                  field
!> @param[in]     map_w3            DoF map for the W3 field
subroutine panel_edge_w2h_average_code( nlayers,                              &
                                        perp_w2_at_w2_x,                      &
                                        perp_w2_at_w2_y,                      &
                                        perp_p_at_w1v_x,                      &
                                        perp_p_at_w1v_y,                      &
                                        parl_p_at_w3_x,                       &
                                        parl_p_at_w3_y,                       &
                                        panel_edge_dist_W,                    &
                                        panel_edge_dist_E,                    &
                                        panel_edge_dist_S,                    &
                                        panel_edge_dist_N,                    &
                                        face_selector_ew,                     &
                                        face_selector_ns,                     &
                                        remap_depth,                          &
                                        include_perp_p,                       &
                                        ndf_wr, undf_wr, map_wr,              &
                                        ndf_wv, undf_wv, map_wv,              &
                                        ndf_pid, undf_pid, map_pid,           &
                                        ndf_w3, undf_w3, map_w3              &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_wv, undf_wv
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_wv(ndf_wv)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_x(undf_wr)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_y(undf_wr)
  real(kind=r_tran),   intent(in)    :: perp_p_at_w1v_x(undf_wv)
  real(kind=r_tran),   intent(in)    :: perp_p_at_w1v_y(undf_wv)
  real(kind=r_tran),   intent(in)    :: parl_p_at_w3_x(undf_w3)
  real(kind=r_tran),   intent(in)    :: parl_p_at_w3_y(undf_w3)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)
  integer(kind=i_def), intent(in)    :: remap_depth
  logical(kind=l_def), intent(in)    :: include_perp_p

  ! Internal variables
  integer(kind=i_def) :: face, pid_idx, nvert
  integer(kind=i_def) :: ew_sel, ns_sel, j
  integer(kind=i_def) :: wr_idx, sw_idx, se_idx, nw_idx, ne_idx
  integer(kind=i_def) :: own_w3_idx
  logical(kind=l_def) :: near_panel_edge_x, near_panel_edge_y

  nvert = nlayers - 1

  pid_idx = map_pid(1)

  sw_idx = map_wv(SW - SW + 1)
  se_idx = map_wv(SE - SW + 1)
  nw_idx = map_wv(NW - SW + 1)
  ne_idx = map_wv(NE - SW + 1)

  ew_sel = face_selector_ew(pid_idx)
  ns_sel = face_selector_ns(pid_idx)

  ! Away from a panel edge the W1V staging reduces to the identity value and
  ! this average would be a no-op, so gate each direction behind its own
  ! panel_edge_dist/remap_depth test (mirroring the weights kernel's
  ! near_panel_edge_x/y): perp_w2_at_w2_x (S/N faces) only needs correcting
  ! near a W/E panel edge, perp_w2_at_w2_y (W/E faces) only near a S/N one
  near_panel_edge_x = (                                                     &
      ABS(panel_edge_dist_W(pid_idx)) <= remap_depth+1                        &
      .or. ABS(panel_edge_dist_E(pid_idx)) <= remap_depth+1                   &
  )
  near_panel_edge_y = (                                                     &
      ABS(panel_edge_dist_S(pid_idx)) <= remap_depth+1                        &
      .or. ABS(panel_edge_dist_N(pid_idx)) <= remap_depth+1                   &
  )

  own_w3_idx = map_w3(1)

  if (include_perp_p) then
    if (near_panel_edge_x) then
      wr_idx = map_wr(W)
      perp_w2_at_w2_x(wr_idx : wr_idx+nvert) =                                  &
          perp_w2_at_w2_x(wr_idx : wr_idx+nvert)                               &
          + 0.5_r_tran * parl_p_at_w3_x(own_w3_idx : own_w3_idx+nvert)
      wr_idx = map_wr(E)
      perp_w2_at_w2_x(wr_idx : wr_idx+nvert) =                                  &
          perp_w2_at_w2_x(wr_idx : wr_idx+nvert)                               &
          + 0.5_r_tran * parl_p_at_w3_x(own_w3_idx : own_w3_idx+nvert)
    end if
    if (near_panel_edge_y) then
      wr_idx = map_wr(S)
      perp_w2_at_w2_y(wr_idx : wr_idx+nvert) =                                  &
          perp_w2_at_w2_y(wr_idx : wr_idx+nvert)                               &
          + 0.5_r_tran * parl_p_at_w3_y(own_w3_idx : own_w3_idx+nvert)
      wr_idx = map_wr(N)
      perp_w2_at_w2_y(wr_idx : wr_idx+nvert) =                                  &
          perp_w2_at_w2_y(wr_idx : wr_idx+nvert)                               &
          + 0.5_r_tran * parl_p_at_w3_y(own_w3_idx : own_w3_idx+nvert)
    end if
  end if

  ! Loop only over the faces this column owns (per face_selector_ew/ns), so
  ! that each shared W2H face is written by exactly one of its two
  ! neighbouring columns. perp_w2_at_w2_{x,y} use GH_INC, so the average is
  ! added to (rather than overwriting) the field's current value
  do j = 1, ABS(ew_sel) + ABS(ns_sel)
    face = face_from_face_selector(j, ew_sel, ns_sel)

    ! TODO: to remove
    write(log_scratch_space, *) 'DEBUG AVERAGE', pid_idx, face
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    select case (face)
    case (S)
      if (near_panel_edge_x) then
        wr_idx = map_wr(face)

        write(log_scratch_space, *) 'DEBUG AVERAGE S SE', pid_idx, face, se_idx
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        write(log_scratch_space, *) 'DEBUG AVERAGE S SE', perp_p_at_w1v_x(se_idx : se_idx+nvert)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        write(log_scratch_space, *) 'DEBUG AVERAGE S SW', pid_idx, face, sw_idx
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        write(log_scratch_space, *) 'DEBUG AVERAGE S SW', perp_p_at_w1v_x(sw_idx : sw_idx+nvert)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        ! TODO: to remove
        write(log_scratch_space, *) 'DEBUG AVERAGE', pid_idx, face, &
          perp_p_at_w1v_x(sw_idx : sw_idx+nvert), perp_p_at_w1v_x(se_idx : se_idx+nvert)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        perp_w2_at_w2_x(wr_idx : wr_idx+nvert) =                             &
            perp_w2_at_w2_x(wr_idx : wr_idx+nvert)                          &
            + 0.5_r_tran * ( perp_p_at_w1v_x(sw_idx : sw_idx+nvert)         &
                           + perp_p_at_w1v_x(se_idx : se_idx+nvert) )
      end if
    case (N)
      if (near_panel_edge_x) then
        wr_idx = map_wr(face)

        ! TODO: to remove
        write(log_scratch_space, *) 'DEBUG AVERAGE', pid_idx, face, &
          perp_p_at_w1v_x(nw_idx : nw_idx+nvert), perp_p_at_w1v_x(ne_idx : ne_idx+nvert)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        perp_w2_at_w2_x(wr_idx : wr_idx+nvert) =                             &
            perp_w2_at_w2_x(wr_idx : wr_idx+nvert)                          &
            + 0.5_r_tran * ( perp_p_at_w1v_x(nw_idx : nw_idx+nvert)         &
                           + perp_p_at_w1v_x(ne_idx : ne_idx+nvert) )
      end if

    case (W)
      if (near_panel_edge_y) then
        wr_idx = map_wr(face)

        ! TODO: to remove
        write(log_scratch_space, *) 'DEBUG AVERAGE', pid_idx, face, &
          perp_p_at_w1v_y(sw_idx : sw_idx+nvert), perp_p_at_w1v_y(nw_idx : nw_idx+nvert)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        perp_w2_at_w2_y(wr_idx : wr_idx+nvert) =                             &
            perp_w2_at_w2_y(wr_idx : wr_idx+nvert)                          &
            + 0.5_r_tran * ( perp_p_at_w1v_y(sw_idx : sw_idx+nvert)         &
                           + perp_p_at_w1v_y(nw_idx : nw_idx+nvert) )
      end if
    case (E)
      if (near_panel_edge_y) then
        wr_idx = map_wr(face)

        ! TODO: to remove
        write(log_scratch_space, *) 'DEBUG AVERAGE', pid_idx, face, &
          perp_p_at_w1v_y(se_idx : se_idx+nvert), perp_p_at_w1v_y(ne_idx : ne_idx+nvert)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        perp_w2_at_w2_y(wr_idx : wr_idx+nvert) =                             &
            perp_w2_at_w2_y(wr_idx : wr_idx+nvert)                          &
            + 0.5_r_tran * ( perp_p_at_w1v_y(se_idx : se_idx+nvert)         &
                           + perp_p_at_w1v_y(ne_idx : ne_idx+nvert) )
      end if
    end select
  end do

end subroutine panel_edge_w2h_average_code

end module panel_edge_w2h_average_kernel_mod
