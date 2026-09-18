!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Remap a W2H field near a cubed-sphere panel corner using the
!!        corner-specific weights/indices computed by
!!        panel_corner_remap_w2h_weights_kernel_mod.
module panel_corner_remap_w2h_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER,                          &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_1,                    &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 CELL_COLUMN,                                  &
                                 STENCIL, CROSS2D
use fs_continuity_mod,      only: W2h
use constants_mod,          only: r_def, r_tran, i_def, l_def
use reference_element_mod,  only: W, S, N, E
use sci_face_selector_support_mod, only: face_from_face_selector

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: panel_corner_remap_w2h_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/                                         &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! remapped_flux_x
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! remapped_flux_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H, STENCIL(CROSS2D)),     & ! field_for_x
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H, STENCIL(CROSS2D)),     & ! field_for_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! weights_x_parallel
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! weights_x_perp
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! indices_x
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! weights_y_parallel
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! weights_y_perp
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! indices_y
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_edge_dist_W/E/S/N
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ew
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ns
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! depth
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                              & ! ndata
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: panel_corner_remap_w2h_code
end type

public :: panel_corner_remap_w2h_code, panel_corner_remap_w2h_1d

contains

subroutine panel_corner_remap_w2h_code( nlayers,                                 &
                                      remapped_flux_x,                        &
                                      remapped_flux_y,                        &
                                      field_for_x,                            &
                                      wsx_stencil_size,                       &
                                      wsx_max_length,                         &
                                      wsx_stencil_map,                        &
                                      field_for_y,                            &
                                      wsy_stencil_size,                       &
                                      wsy_max_length,                         &
                                      wsy_stencil_map,                        &
                                      weights_x_parallel,                    &
                                      weights_x_perp,                        &
                                      indices_x,                             &
                                      weights_y_parallel,                    &
                                      weights_y_perp,                        &
                                      indices_y,                             &
                                      panel_edge_dist_W,                     &
                                      panel_edge_dist_E,                     &
                                      panel_edge_dist_S,                     &
                                      panel_edge_dist_N,                     &
                                      face_selector_ew,                      &
                                      face_selector_ns,                      &
                                      depth,                                 &
                                      ndata,                                 &
                                      ndf_wr, undf_wr, map_wr,               &
                                      ndf_ws, undf_ws, map_ws,               &
                                      ndf_ww, undf_ww, map_ww,               &
                                      ndf_pid, undf_pid, map_pid             &
  )

  implicit none

  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: wsx_stencil_size(4), wsy_stencil_size(4)
  integer(kind=i_def), intent(in)    :: wsx_max_length, wsy_max_length
  integer(kind=i_def), intent(in)    :: depth, ndata
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_ws, undf_ws
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_ws(ndf_ws)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wsx_stencil_map(ndf_ws, wsx_max_length, 4)
  integer(kind=i_def), intent(in)    :: wsy_stencil_map(ndf_ws, wsy_max_length, 4)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)

  real(kind=r_tran),   intent(inout) :: remapped_flux_x(undf_wr)
  real(kind=r_tran),   intent(inout) :: remapped_flux_y(undf_wr)
  real(kind=r_tran),   intent(in)    :: field_for_x(undf_ws)
  real(kind=r_tran),   intent(in)    :: field_for_y(undf_ws)
  real(kind=r_tran),   intent(in)    :: weights_x_parallel(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_x_perp(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_y_parallel(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_y_perp(undf_ww)
  integer(kind=i_def), intent(in)    :: indices_x(undf_ww)
  integer(kind=i_def), intent(in)    :: indices_y(undf_ww)

  logical(kind=l_def) :: compute_x, compute_y
  integer(kind=i_def) :: pid_idx, j, face, ew_sel, ns_sel, perp_dir_x, perp_dir_y

  pid_idx = map_pid(1)
  ew_sel = face_selector_ew(pid_idx)
  ns_sel = face_selector_ns(pid_idx)

  remapped_flux_x = field_for_x
  remapped_flux_y = field_for_y

  compute_x = (                                                                &
    (ABS(panel_edge_dist_W(pid_idx)) <= depth .or.                            &
     ABS(panel_edge_dist_E(pid_idx)) <= depth) .and. (                         &
      ABS(panel_edge_dist_S(pid_idx)) <= depth .or.                            &
      ABS(panel_edge_dist_N(pid_idx)) <= depth                                  &
    )                                                                          &
  )
  compute_y = (                                                                &
    (ABS(panel_edge_dist_S(pid_idx)) <= depth .or.                            &
     ABS(panel_edge_dist_N(pid_idx)) <= depth) .and. (                         &
      ABS(panel_edge_dist_W(pid_idx)) <= depth .or.                            &
      ABS(panel_edge_dist_E(pid_idx)) <= depth                                  &
    )                                                                          &
  )

  if (ABS(panel_edge_dist_W(pid_idx)) < ABS(panel_edge_dist_E(pid_idx))) then
    perp_dir_x = W
  else
    perp_dir_x = E
  end if

  if (ABS(panel_edge_dist_S(pid_idx)) < ABS(panel_edge_dist_N(pid_idx))) then
    perp_dir_y = S
  else
    perp_dir_y = N
  end if

  do j = 1, ABS(ew_sel) + ABS(ns_sel)
    face = face_from_face_selector(j, ew_sel, ns_sel)

    if (face == W .or. face == E) then
      if (compute_y) then
        call panel_corner_remap_w2h_1d(                                          &
             nlayers, remapped_flux_y, field_for_y, face, perp_dir_y,          &
             wsy_max_length,                                                     &
             wsy_stencil_size(W), wsy_stencil_map(:,:,W),                        &
             wsy_stencil_size(E), wsy_stencil_map(:,:,E),                        &
             weights_y_parallel, weights_y_perp, indices_y, ndata,              &
             ndf_wr, undf_wr, map_wr,                                            &
             ndf_ws, undf_ws, map_ws,                                            &
             ndf_ww, undf_ww, map_ww                                             &
        )
      end if
      ! TODO: to remove
      if (pid_idx == 864) then
        write(log_scratch_space, *) 'PANEL CORNER REMAP W2H: ', face, j, &
          remapped_flux_y(map_wr(face)), field_for_y(map_ws(W)), field_for_y(map_ws(E)), &
          field_for_y(map_ws(S)), field_for_y(map_ws(N))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
      end if
    else if (face == S .or. face == N) then
      if (compute_x) then
        call panel_corner_remap_w2h_1d(                                          &
             nlayers, remapped_flux_x, field_for_x, face, perp_dir_x,          &
             wsx_max_length,                                                     &
             wsx_stencil_size(S), wsx_stencil_map(:,:,S),                        &
             wsx_stencil_size(N), wsx_stencil_map(:,:,N),                        &
             weights_x_parallel, weights_x_perp, indices_x, ndata,              &
             ndf_wr, undf_wr, map_wr,                                            &
             ndf_ws, undf_ws, map_ws,                                            &
             ndf_ww, undf_ww, map_ww                                             &
        )
      end if
    end if
  end do

end subroutine panel_corner_remap_w2h_code

subroutine panel_corner_remap_w2h_1d( nlayers,                                  &
                                    remapped_field, field,                     &
                                    target_df, perp_df,                       &
                                    ws_max_length,                            &
                                    ws_stencil_size_l, ws_stencil_l,          &
                                    ws_stencil_size_r, ws_stencil_r,          &
                                    weights_parallel, weights_perp, indices, &
                                    ndata,                                    &
                                    ndf_wr, undf_wr, map_wr,                  &
                                    ndf_ws, undf_ws, map_ws,                  &
                                    ndf_ww, undf_ww, map_ww                   &
  )

  implicit none

  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: target_df, perp_df
  integer(kind=i_def), intent(in)    :: ws_stencil_size_l
  integer(kind=i_def), intent(in)    :: ws_stencil_size_r
  integer(kind=i_def), intent(in)    :: ws_max_length
  integer(kind=i_def), intent(in)    :: ndata
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_ws, undf_ws
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_ws(ndf_ws)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: ws_stencil_l(ndf_ws, ws_stencil_size_l)
  integer(kind=i_def), intent(in)    :: ws_stencil_r(ndf_ws, ws_stencil_size_r)

  real(kind=r_tran),   intent(inout) :: remapped_field(undf_wr)
  real(kind=r_tran),   intent(in)    :: field(undf_ws)
  real(kind=r_tran),   intent(in)    :: weights_parallel(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_perp(undf_ww)
  integer(kind=i_def), intent(in)    :: indices(undf_ww)

  integer(kind=i_def) :: i, w_idx, r_idx, f_idx_parallel, f_idx_perp
  integer(kind=i_def) :: nvert
  integer(kind=i_def) :: stencil_1d_parallel(2*ws_max_length-1)
  integer(kind=i_def) :: stencil_1d_perp(2*ws_max_length-1)

  nvert = nlayers - 1
  stencil_1d_parallel(:) = 0
  stencil_1d_perp(:) = 0

  do i = 1, ws_stencil_size_l
    stencil_1d_parallel(i) = ws_stencil_l(target_df, i)
    stencil_1d_perp(i) = ws_stencil_l(perp_df, i)
  end do
  do i = 1, ws_stencil_size_r - 1
    stencil_1d_parallel(i+ws_stencil_size_l) = ws_stencil_r(target_df, i+1)
    stencil_1d_perp(i+ws_stencil_size_l) = ws_stencil_r(perp_df, i+1)
  end do

  w_idx = map_ww(target_df)
  r_idx = map_wr(target_df)
  remapped_field(r_idx:r_idx+nvert) = 0.0_r_tran

  do i = 0, ndata - 1
    f_idx_parallel = stencil_1d_parallel(indices(w_idx+i))
    f_idx_perp = stencil_1d_perp(indices(w_idx+i))
    remapped_field(r_idx:r_idx+nvert) = remapped_field(r_idx:r_idx+nvert)      &
      + weights_parallel(w_idx+i) * field(f_idx_parallel:f_idx_parallel+nvert) &
      + weights_perp(w_idx+i) * field(f_idx_perp:f_idx_perp+nvert)
  end do

end subroutine panel_corner_remap_w2h_1d

end module panel_corner_remap_w2h_kernel_mod
