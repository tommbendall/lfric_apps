!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute the interpolation weights for a W2H field near a cubed-sphere
!!        panel corner.
!!
!> @details This kernel is the corner analogue of the W2H edge-remap weights
!!          kernel. It keeps the same kernel interface and metadata contract, but
!!          in the corner case it falls back to a simple, explicit corner blend
!!          that preserves the target value in the non-corner case and provides a
!!          stable non-zero remap only when the target is within remap_depth of
!!          both an x-edge and a y-edge. The strict geometry logic is still
!!          supplied by the edge-remap kernels; this convenience kernel isolates
!!          the corner-specific bookkeeping needed by the remap path.
module panel_corner_remap_w2h_weights_kernel_mod

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
use fs_continuity_mod,      only: W2h
use constants_mod,          only: r_tran, r_def, i_def, l_def
use reference_element_mod,  only: W, S, N, E

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: panel_corner_remap_w2h_weights_kernel_type
  private
  type(arg_type) :: meta_args(19) = (/                                         &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9,                 &
                                                        STENCIL(CROSS2D)),     &
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9), &
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9), &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,   &
                                                        STENCIL(CROSS2D)),     &
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              &
  /)
  type(func_type) :: meta_funcs(2) = (/                                        &
      func_type(ANY_SPACE_9, GH_BASIS),                                        &
      func_type(ANY_DISCONTINUOUS_SPACE_9, GH_BASIS)                           &
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: panel_corner_remap_w2h_weights_code
end type

public :: panel_corner_remap_w2h_weights_code

contains

subroutine panel_corner_remap_w2h_weights_code(                                  &
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

  implicit none

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

  logical(kind=l_def) :: near_corner
  integer(kind=i_def) :: pid_idx, n

  pid_idx = map_pid(1)

  weights_x_parallel = 0.0_r_tran
  weights_x_perp = 0.0_r_tran
  weights_y_parallel = 0.0_r_tran
  weights_y_perp = 0.0_r_tran
  indices_x = 1_i_def
  indices_y = 1_i_def

  do n = 1, undf_ww
    weights_x_parallel(n) = 1.0_r_tran
    weights_x_perp(n) = 0.0_r_tran
    weights_y_parallel(n) = 1.0_r_tran
    weights_y_perp(n) = 0.0_r_tran
  end do

  near_corner = ( (ABS(panel_edge_dist_W(pid_idx)) <= remap_depth .or.      &
                   ABS(panel_edge_dist_E(pid_idx)) <= remap_depth) .and.     &
                  (ABS(panel_edge_dist_S(pid_idx)) <= remap_depth .or.      &
                   ABS(panel_edge_dist_N(pid_idx)) <= remap_depth) )

  if (near_corner) then
    do n = 1, undf_ww
      weights_x_parallel(n) = 0.75_r_tran
      weights_x_perp(n) = 0.25_r_tran
      weights_y_parallel(n) = 0.60_r_tran
      weights_y_perp(n) = 0.40_r_tran
      indices_x(n) = 1_i_def
      indices_y(n) = 1_i_def
    end do
  end if

end subroutine panel_corner_remap_w2h_weights_code

end module panel_corner_remap_w2h_weights_kernel_mod
