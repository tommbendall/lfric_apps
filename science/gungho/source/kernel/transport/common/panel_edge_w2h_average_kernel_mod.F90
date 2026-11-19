!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Add the W3 centre staging values onto the "perp_p" W2H faces near
!!        cubed-sphere panel edges, giving part of the "perp_w2" staged value
!!        consumed by panel_edge_remap_w2h_kernel_mod.
!!
!> @details This implements Route B of the instructions doc (section 2.4):
!!          parl_p_at_w3_x/y hold the "parl_p" value interpolated onto each
!!          cell's W3 centres. On interior faces, each adjacent cell adds
!!          half of its value to its own panel label, preserving the existing
!!          average. At a panel edge the adjacent cells have different labels,
!!          so each adds its full value to its own labelled face.
!!
!!          The "parl_p" faces (S/N for x, W/E for y) are written directly by
!!          panel_edge_interpolate_w2h_kernel_mod.
!!
!!          Each direction remains gated by its own panel-edge-distance test.
!!          Labelled outputs use GH_INC; no halo exchange is added here.
module panel_edge_w2h_average_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                   &
                                 GH_FIELD, GH_SCALAR,                        &
                                 GH_REAL, GH_INTEGER, GH_LOGICAL,             &
                                 GH_READ, GH_INC,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                  &
                                 CELL_COLUMN
use fs_continuity_mod,     only: W2h, W3
use constants_mod,         only: r_def, r_tran, i_def, l_def
use reference_element_mod, only: W, S, E, N
use panel_edge_support_mod, only: panel_neighbour, panel_axis_label

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!! Psy layer.
type, public, extends(kernel_type) :: panel_edge_w2h_average_kernel_type
  private
    type(arg_type) :: meta_args(7) = (/                                          &
      arg_type(GH_FIELD*3, GH_REAL,    GH_INC,   W2h),                         & ! perp_w2_at_w2
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                          & ! parl_p_at_w3_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                          & ! parl_p_at_w3_y
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,                                &
               ANY_DISCONTINUOUS_SPACE_3),                                     & ! panel_edge_dist_W/E/S/N
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   & ! panel_id
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

!> @brief Add the W3 centre staging values onto the "perp_p" W2H faces.
!> @param[in]     nlayers           Number of layers
!> @param[in,out] perp_w2_at_w2_1   "perp_w2" staged value for label 1
!> @param[in,out] perp_w2_at_w2_2   "perp_w2" staged value for label 2
!> @param[in,out] perp_w2_at_w2_3   "perp_w2" staged value for label 3
!> @param[in]     parl_p_at_w3_x    "parl_p" value interpolated onto the W3
!!                                  cell centre for the x remap
!> @param[in]     parl_p_at_w3_y    "parl_p" value interpolated onto the W3
!!                                  cell centre for the y remap
!> @param[in]     panel_edge_dist_W Distance of this column to the W panel edge
!> @param[in]     panel_edge_dist_E Distance of this column to the E panel edge
!> @param[in]     panel_edge_dist_S Distance of this column to the S panel edge
!> @param[in]     panel_edge_dist_N Distance of this column to the N panel edge
!> @param[in]     panel_id          Panel ID for this column
!> @param[in]     remap_depth       Depth (in columns) to which remapping
!!                                  applies
!> @param[in]     include_perp_p    Whether to include the perp_p contribution
!> @param[in]     ndf_wr            Num DoFs per cell for the W2H field
!> @param[in]     undf_wr           Num DoFs for this partition for the W2H
!!                                  field
!> @param[in]     map_wr            DoF map for the W2H field
!> @param[in]     ndf_w3            Num DoFs per cell for the W3 field
!> @param[in]     undf_w3           Num DoFs for this partition for the W3
!!                                  field
!> @param[in]     map_w3            DoF map for the W3 field
!> @param[in]     ndf_pid           Num DoFs per cell for panel ID field
!> @param[in]     undf_pid          Num DoFs for this partition for panel ID
!> @param[in]     map_pid           DoF map for panel ID field
subroutine panel_edge_w2h_average_code( nlayers,                              &
                                        perp_w2_at_w2_1,                      &
                                        perp_w2_at_w2_2,                      &
                                        perp_w2_at_w2_3,                      &
                                        parl_p_at_w3_x,                       &
                                        parl_p_at_w3_y,                       &
                                        panel_edge_dist_W,                    &
                                        panel_edge_dist_E,                    &
                                        panel_edge_dist_S,                    &
                                        panel_edge_dist_N,                    &
                                        panel_id,                             &
                                        remap_depth,                          &
                                        include_perp_p,                       &
                                        ndf_wr, undf_wr, map_wr,              &
                                        ndf_w3, undf_w3, map_w3,              &
                                        ndf_pid, undf_pid, map_pid            &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_1(undf_wr)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_2(undf_wr)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_3(undf_wr)
  real(kind=r_tran),   intent(in)    :: parl_p_at_w3_x(undf_w3)
  real(kind=r_tran),   intent(in)    :: parl_p_at_w3_y(undf_w3)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: remap_depth
  logical(kind=l_def), intent(in)    :: include_perp_p

  ! Internal variables
  integer(kind=i_def) :: pid_idx, nvert, owned_panel, face
  integer(kind=i_def) :: wr_idx, face_edge_dist
  integer(kind=i_def) :: edge_direction_x, edge_direction_y, label_x, label_y
  real(kind=r_tran) :: contribution(nlayers)
  real(kind=r_tran) :: contribution_scale
  integer(kind=i_def) :: own_w3_idx
  logical(kind=l_def) :: near_panel_edge_x, near_panel_edge_y

  nvert = nlayers - 1

  pid_idx = map_pid(1)
  owned_panel = INT(panel_id(pid_idx), i_def)

  ! Away from a panel edge the contribution is the identity, so gate each
  ! direction behind its own panel_edge_dist/remap_depth test (mirroring the
  ! weights kernel's near_panel_edge_x/y)
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
      if (ABS(panel_edge_dist_W(pid_idx)) < ABS(panel_edge_dist_E(pid_idx))) then
        edge_direction_x = W
      else
        edge_direction_x = E
      end if
      label_x = panel_axis_label(panel_neighbour(owned_panel, edge_direction_x))
      do face = W, E, E-W
        if (face == W) then
          face_edge_dist = panel_edge_dist_W(pid_idx)
        else
          face_edge_dist = panel_edge_dist_E(pid_idx)
        end if
        contribution_scale = 0.5_r_tran
        if (face_edge_dist == 1) contribution_scale = 1.0_r_tran
        contribution = contribution_scale *                                  &
            parl_p_at_w3_x(own_w3_idx : own_w3_idx+nvert)
        wr_idx = map_wr(face)
        call add_labelled_contribution(label_x, wr_idx, contribution)
      end do
    end if
    if (near_panel_edge_y) then
      if (ABS(panel_edge_dist_S(pid_idx)) < ABS(panel_edge_dist_N(pid_idx))) then
        edge_direction_y = S
      else
        edge_direction_y = N
      end if
      label_y = panel_axis_label(panel_neighbour(owned_panel, edge_direction_y))
      do face = S, N, N-S
        if (face == S) then
          face_edge_dist = panel_edge_dist_S(pid_idx)
        else
          face_edge_dist = panel_edge_dist_N(pid_idx)
        end if
        contribution_scale = 0.5_r_tran
        if (face_edge_dist == 1) contribution_scale = 1.0_r_tran
        contribution = contribution_scale *                                  &
            parl_p_at_w3_y(own_w3_idx : own_w3_idx+nvert)
        wr_idx = map_wr(face)
        call add_labelled_contribution(label_y, wr_idx, contribution)
      end do
    end if
  end if

contains

  subroutine add_labelled_contribution(label, idx, values)

    implicit none

    integer(kind=i_def), intent(in) :: label, idx
    real(kind=r_tran), intent(in) :: values(nvert+1)

    select case (label)
    case (1)
      perp_w2_at_w2_1(idx:idx+nvert) =                                       &
          perp_w2_at_w2_1(idx:idx+nvert) + values
    case (2)
      perp_w2_at_w2_2(idx:idx+nvert) =                                       &
          perp_w2_at_w2_2(idx:idx+nvert) + values
    case (3)
      perp_w2_at_w2_3(idx:idx+nvert) =                                       &
          perp_w2_at_w2_3(idx:idx+nvert) + values
    end select

  end subroutine add_labelled_contribution

end subroutine panel_edge_w2h_average_code

end module panel_edge_w2h_average_kernel_mod
