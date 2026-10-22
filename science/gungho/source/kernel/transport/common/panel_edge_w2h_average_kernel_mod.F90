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
!!          parl_p_at_w3_x/y hold the "parl_p" value interpolated
!!          onto the W3 cell centres. For the W/E faces (the "perp_p"
!!          targets of the x remap), each column adds half of its own
!!          parl_p_at_w3_x centre value to both its W and E faces; the
!!          neighbouring column across each face independently adds its own
!!          half when it is visited (no stencil read of the neighbour is
!!          needed). Likewise for the S/N faces (the "perp_p" targets of
!!          the y remap) using parl_p_at_w3_y.
!!
!!          The "parl_p" faces (S/N for x, W/E for y) are written directly by
!!          panel_edge_interpolate_w2h_kernel_mod.
!!
!!          Only near a panel edge does the contribution differ from the
!!          identity, so each direction is gated behind its own
!!          `panel_edge_dist`/`remap_depth` test (mirroring the weights
!!          kernel's `near_panel_edge_x/y`).
!!
!!          perp_w2_at_w2_{x,y} use GH_INC so that each face's two
!!          neighbouring columns can each add their own half independently.
module panel_edge_w2h_average_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                   &
                                 GH_FIELD, GH_SCALAR,                        &
                                 GH_REAL, GH_INTEGER, GH_LOGICAL,             &
                                 GH_READ, GH_INC,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                  &
                                 CELL_COLUMN
use fs_continuity_mod,     only: W2h, W3
use constants_mod,         only: r_tran, i_def, l_def
use reference_element_mod, only: W, S, E, N

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
      arg_type(GH_FIELD,   GH_REAL,    GH_INC,   W2h),                         & ! perp_w2_at_w2_x
      arg_type(GH_FIELD,   GH_REAL,    GH_INC,   W2h),                         & ! perp_w2_at_w2_y
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                          & ! parl_p_at_w3_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                          & ! parl_p_at_w3_y
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,                                &
               ANY_DISCONTINUOUS_SPACE_3),                                     & ! panel_edge_dist_W/E/S/N
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
!> @param[in,out] perp_w2_at_w2_x   "perp_w2" staged value for the x remap,
!!                                  incremented at W/E faces
!> @param[in,out] perp_w2_at_w2_y   "perp_w2" staged value for the y remap,
!!                                  incremented at S/N faces
!> @param[in]     parl_p_at_w3_x    "parl_p" value interpolated onto the W3
!!                                  cell centre for the x remap
!> @param[in]     parl_p_at_w3_y    "parl_p" value interpolated onto the W3
!!                                  cell centre for the y remap
!> @param[in]     panel_edge_dist_W Distance of this column to the W panel edge
!> @param[in]     panel_edge_dist_E Distance of this column to the E panel edge
!> @param[in]     panel_edge_dist_S Distance of this column to the S panel edge
!> @param[in]     panel_edge_dist_N Distance of this column to the N panel edge
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
                                        perp_w2_at_w2_x,                      &
                                        perp_w2_at_w2_y,                      &
                                        parl_p_at_w3_x,                       &
                                        parl_p_at_w3_y,                       &
                                        panel_edge_dist_W,                    &
                                        panel_edge_dist_E,                    &
                                        panel_edge_dist_S,                    &
                                        panel_edge_dist_N,                    &
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

  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_x(undf_wr)
  real(kind=r_tran),   intent(inout) :: perp_w2_at_w2_y(undf_wr)
  real(kind=r_tran),   intent(in)    :: parl_p_at_w3_x(undf_w3)
  real(kind=r_tran),   intent(in)    :: parl_p_at_w3_y(undf_w3)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: remap_depth
  logical(kind=l_def), intent(in)    :: include_perp_p

  ! Internal variables
  integer(kind=i_def) :: pid_idx, nvert
  integer(kind=i_def) :: wr_idx
  integer(kind=i_def) :: own_w3_idx
  logical(kind=l_def) :: near_panel_edge_x, near_panel_edge_y

  nvert = nlayers - 1

  pid_idx = map_pid(1)

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

end subroutine panel_edge_w2h_average_code

end module panel_edge_w2h_average_kernel_mod
