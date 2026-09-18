!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! GitHub Copilot.
!-------------------------------------------------------------------------------
!> @brief Remap a W2H flux/wind field at the edges of cubed-sphere panels to
!!        the corresponding neighbouring panels, using the weights/indices
!!        computed by panel_edge_remap_w2h_weights_kernel_mod.
!!
!> @details Two independent remaps are produced, mirroring
!!          panel_edge_remap_kernel_mod's scalar x/y split:
!!            - remapped_flux_x: corrected S/N (v-type) face values, valid
!!              near a W or E panel edge. The W/E (u-type) face values of
!!              this output are simply copied from the input field, since the
!!              x remap does not apply to them.
!!            - remapped_flux_y: corrected W/E (u-type) face values, valid
!!              near a S or N panel edge. The S/N face values of this output
!!              are simply copied from the input field.
!!          Away from the relevant edge, the weights already default to the
!!          identity (see panel_edge_remap_w2h_weights_kernel_mod), so both
!!          outputs equal the input field everywhere they are not actively
!!          correcting a face value. This means the two outputs can later be
!!          combined without any additional merge/select logic, via:
!!             corrected_flux = remapped_flux_x + remapped_flux_y - ref_flux
!!
!!          At each target S/N (or W/E) DoF, the remapped value is built from
!!          a donor-cell stencil tangential to the relevant edge:
!!             remapped = sum_n [ weights_parallel(n) * donor_same_face(n)    &
!!                               + weights_perp(n)     * donor_cross_face(n) ]
!!          where donor_same_face is the donor cell's own value at the same
!!          face orientation as the target DoF, and donor_cross_face is the
!!          donor cell's value at the reference face of the opposite
!!          orientation matching the panel edge actually being crossed (W or
!!          E for the x remap, S or N for the y remap), i.e. the same
!!          edge_direction used by panel_edge_remap_w2h_weights_kernel_mod to
!!          select the neighbouring panel for the rotation factors.
module panel_edge_remap_w2h_kernel_mod

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
! TODO: to remove
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_DEBUG

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_edge_remap_w2h_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/                                         &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! remapped_flux_x
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! remapped_flux_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H, STENCIL(CROSS2D)),     & ! field_for_x (ref_flux)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H, STENCIL(CROSS2D)),     & ! field_for_y (ref_flux)
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
  procedure, nopass :: panel_edge_remap_w2h_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_remap_w2h_code

contains

!> @brief Remap a W2H flux/wind field at the edges of cubed-sphere panels.
!> @param[in]     nlayers               Number of layers
!> @param[in,out] remapped_flux_x       Field remapped for the x (S/N) correction
!> @param[in,out] remapped_flux_y       Field remapped for the y (W/E) correction
!> @param[in]     field_for_x           Field to be remapped in the x direction
!> @param[in]     wsx_stencil_size      Size of field stencil for x (num cells)
!> @param[in]     wsx_max_length        Maximum stencil branch length for x
!> @param[in]     wsx_stencil_map       DoF map for the field stencil (x)
!> @param[in]     field_for_y           Field to be remapped in the y direction
!> @param[in]     wsy_stencil_size      Size of field stencil for y (num cells)
!> @param[in]     wsy_max_length        Maximum stencil branch length for y
!> @param[in]     wsy_stencil_map       DoF map for the field stencil (y)
!> @param[in]     weights_x_parallel    Same-orientation donor weights for x remap
!> @param[in]     weights_x_perp        Cross-orientation donor weights for x remap
!> @param[in]     indices_x             Donor cell indices for the x remap
!> @param[in]     weights_y_parallel    Same-orientation donor weights for y remap
!> @param[in]     weights_y_perp        Cross-orientation donor weights for y remap
!> @param[in]     indices_y             Donor cell indices for the y remap
!> @param[in]     panel_edge_dist_W     Distance of this column to the W panel edge
!> @param[in]     panel_edge_dist_E     Distance of this column to the E panel edge
!> @param[in]     panel_edge_dist_S     Distance of this column to the S panel edge
!> @param[in]     panel_edge_dist_N     Distance of this column to the N panel edge
!> @param[in]     face_selector_ew      East-West face selector, used to ensure
!!                                      each W/E face is only ever written by
!!                                      one of its two neighbouring columns
!> @param[in]     face_selector_ns      North-South face selector, as
!!                                      face_selector_ew for S/N faces
!> @param[in]     depth                 Maximum halo depth to consider
!> @param[in]     ndata                 Number of donor points in the remapping
!> @param[in]     ndf_wr                Num DoFs per cell for remapped fields
!> @param[in]     undf_wr               Num DoFs for this partition for remapped
!!                                      fields
!> @param[in]     map_wr                DoF map for remapped fields
!> @param[in]     ndf_ws                Num DoFs per cell for input flux field
!> @param[in]     undf_ws               Num DoFs for this partition for input
!!                                      flux field
!> @param[in]     map_ws                DoF map for input flux field
!> @param[in]     ndf_ww                Num DoFs per cell for remapping weights
!> @param[in]     undf_ww               Num DoFs for this partition for weights
!> @param[in]     map_ww                DoF map for remapping weights
!> @param[in]     ndf_pid               Num DoFs per cell for panel ID field
!> @param[in]     undf_pid              Num DoFs for this partition for panel ID
!> @param[in]     map_pid               DoF map for panel ID field
subroutine panel_edge_remap_w2h_code( nlayers,                                 &
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
                                      face_selector_ew,                     &
                                      face_selector_ns,                     &
                                      depth,                                 &
                                      ndata,                                 &
                                      ndf_wr, undf_wr, map_wr,               &
                                      ndf_ws, undf_ws, map_ws,               &
                                      ndf_ww, undf_ww, map_ww,               &
                                      ndf_pid, undf_pid, map_pid             &
  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: wsx_stencil_size(4), wsy_stencil_size(4)
  integer(kind=i_def), intent(in)    :: wsx_max_length, wsy_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_ws, undf_ws
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: depth
  integer(kind=i_def), intent(in)    :: ndata

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_ws(ndf_ws)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wsx_stencil_map(ndf_ws, wsx_max_length, 4)
  integer(kind=i_def), intent(in)    :: wsy_stencil_map(ndf_ws, wsy_max_length, 4)

  real(kind=r_tran),   intent(inout) :: remapped_flux_x(undf_wr)
  real(kind=r_tran),   intent(inout) :: remapped_flux_y(undf_wr)
  real(kind=r_tran),   intent(in)    :: field_for_x(undf_ws)
  real(kind=r_tran),   intent(in)    :: field_for_y(undf_ws)
  real(kind=r_tran),   intent(in)    :: weights_x_parallel(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_x_perp(undf_ww)
  integer(kind=i_def), intent(in)    :: indices_x(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_y_parallel(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_y_perp(undf_ww)
  integer(kind=i_def), intent(in)    :: indices_y(undf_ww)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)

  ! Internal variables
  integer(kind=i_def) :: face, pid_idx, nvert
  integer(kind=i_def) :: ew_sel, ns_sel, j
  logical(kind=l_def) :: compute_x, compute_y

  ! Number of vertical levels to loop over (W2H is layer-based, like W3)
  nvert = nlayers - 1

  pid_idx = map_pid(1)

  ! Determine whether to perform remapping --------------------------------------
  compute_x = ABS(panel_edge_dist_W(pid_idx)) <= depth .or.                    &
              ABS(panel_edge_dist_E(pid_idx)) <= depth

  compute_y = ABS(panel_edge_dist_N(pid_idx)) <= depth .or.                    &
              ABS(panel_edge_dist_S(pid_idx)) <= depth

  ! Remap at x edge of panel (corrects S/N face values) ------------------------
  ew_sel = face_selector_ew(pid_idx)
  ns_sel = face_selector_ns(pid_idx)

  ! Loop only over the faces this column owns (per face_selector_ew/ns), so
  ! that each shared W2H face is written by exactly one of its two
  ! neighbouring columns
  do j = 1, ABS(ew_sel) + ABS(ns_sel)
    face = face_from_face_selector(j, ew_sel, ns_sel)
    select case (face)
    case (S, N)
      if (compute_x) then
        write(log_scratch_space, '(A,I0,A,I0,A,E20.12)')                       &
          'PEREGRIN W2H PRE  col=', pid_idx, ' face=', face,                   &
          ' field_for_x=', field_for_x(map_ws(face))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

        call panel_edge_remap_w2h_1d(                                          &
                nlayers, remapped_flux_x, field_for_x, face, W, E,             &
                wsx_max_length,                                                &
                wsx_stencil_size(S), wsx_stencil_map(:,:,S),                   &
                wsx_stencil_size(N), wsx_stencil_map(:,:,N),                   &
                weights_x_parallel, weights_x_perp, indices_x, ndata,          &
                ndf_wr, undf_wr, map_wr,                                       &
                ndf_ws, undf_ws, map_ws,                                       &
                ndf_ww, undf_ww, map_ww                                        &
        )

        write(log_scratch_space, '(A,I0,A,I0,A,E20.12)')                       &
          'PEREGRIN W2H POST col=', pid_idx, ' face=', face,                   &
          ' remapped_flux_x=', remapped_flux_x(map_wr(face))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
      else
        ! Not near a W/E edge: the x remap is the identity
        remapped_flux_x(map_wr(face):map_wr(face)+nvert) =                    &
            field_for_x(map_ws(face):map_ws(face)+nvert)
      end if
      ! The y remap does not apply to S/N faces: copy the field through
      remapped_flux_y(map_wr(face):map_wr(face)+nvert) =                      &
          field_for_y(map_ws(face):map_ws(face)+nvert)

    case (W, E)
      if (compute_y) then
        write(log_scratch_space, '(A,I0,A,I0,A,E20.12)')                       &
          'PEREGRIN W2H PRE  col=', pid_idx, ' face=', face,                   &
          ' field_for_y=', field_for_y(map_ws(face))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

        call panel_edge_remap_w2h_1d(                                          &
                nlayers, remapped_flux_y, field_for_y, face, S, N,             &
                wsy_max_length,                                                &
                wsy_stencil_size(W), wsy_stencil_map(:,:,W),                    &
                wsy_stencil_size(E), wsy_stencil_map(:,:,E),                    &
                weights_y_parallel, weights_y_perp, indices_y, ndata,          &
                ndf_wr, undf_wr, map_wr,                                       &
                ndf_ws, undf_ws, map_ws,                                       &
                ndf_ww, undf_ww, map_ww                                        &
        )

        ! TODO: to remove
        if (pid_idx == 864 .or. pid_idx == 848) then
          write(log_scratch_space, *) 'PANEL EDGE REMAP W2H: ', pid_idx, face, j, &
            remapped_flux_y(map_wr(face)), field_for_y(map_ws(W)), field_for_y(map_ws(E)), &
            field_for_y(map_ws(S)), field_for_y(map_ws(N))
          call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
        end if

        write(log_scratch_space, '(A,I0,A,I0,A,E20.12)')                       &
          'PEREGRIN W2H POST col=', pid_idx, ' face=', face,                   &
          ' remapped_flux_y=', remapped_flux_y(map_wr(face))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
      else
        ! Not near a S/N edge: the y remap is the identity
        remapped_flux_y(map_wr(face):map_wr(face)+nvert) =                    &
            field_for_y(map_ws(face):map_ws(face)+nvert)
      end if
      ! The y remap does not touch S/N faces: copy the x-remapped field here
      remapped_flux_x(map_wr(face):map_wr(face)+nvert) =                      &
          field_for_x(map_ws(face):map_ws(face)+nvert)
    end select
  end do

end subroutine panel_edge_remap_w2h_code

!> @brief Private routine to remap a single W2H target DoF ("target_df") from
!!        a donor-cell stencil tangential to the relevant panel edge.
!!        remapped_field/field are both full W2H fields (all 4 faces); this
!!        routine only writes/reads the "target_df" and "perp_df" faces.
subroutine panel_edge_remap_w2h_1d( nlayers,                                   &
                                    remapped_field, field,                     &
                                    target_df, perp_df1, perp_df2,             &
                                    ws_max_length,                             &
                                    ws_stencil_size_l, ws_stencil_l,           &
                                    ws_stencil_size_r, ws_stencil_r,           &
                                    weights_parallel, weights_perp, indices,   &
                                    ndata,                                     &
                                    ndf_wr, undf_wr, map_wr,                   &
                                    ndf_ws, undf_ws, map_ws,                   &
                                    ndf_ww, undf_ww, map_ww                    &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: target_df, perp_df1, perp_df2
  integer(kind=i_def), intent(in)    :: ws_stencil_size_l
  integer(kind=i_def), intent(in)    :: ws_stencil_size_r
  integer(kind=i_def), intent(in)    :: ws_max_length
  integer(kind=i_def), intent(in)    :: ndata
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_ws, undf_ws

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_ws(ndf_ws)
  integer(kind=i_def), intent(in)    :: ws_stencil_l(ndf_ws, ws_stencil_size_l)
  integer(kind=i_def), intent(in)    :: ws_stencil_r(ndf_ws, ws_stencil_size_r)

  real(kind=r_tran),   intent(inout) :: remapped_field(undf_wr)
  real(kind=r_tran),   intent(in)    :: field(undf_ws)
  real(kind=r_tran),   intent(in)    :: weights_parallel(undf_ww)
  real(kind=r_tran),   intent(in)    :: weights_perp(undf_ww)
  integer(kind=i_def), intent(in)    :: indices(undf_ww)

  ! Internal variables
  integer(kind=i_def) :: i
  integer(kind=i_def) :: w_idx, r_idx
  integer(kind=i_def) :: f_idx_parallel, f_idx_perp1, f_idx_perp2
  integer(kind=i_def) :: nvert
  integer(kind=i_def) :: ncells_in_stencil
  integer(kind=i_def) :: stencil_1d_parallel(2*ws_max_length-1)
  integer(kind=i_def) :: stencil_1d_perp1(2*ws_max_length-1)
  integer(kind=i_def) :: stencil_1d_perp2(2*ws_max_length-1)

  nvert = nlayers - 1

  ! -------------------------------------------------------------------------- !
  ! Create 1D donor stencils (unify the two branches of the CROSS2D stencil),
  ! one sampling the target's own face orientation ("parallel") and one
  ! sampling the fixed cross-orientation reference face ("perp"), from the
  ! same donor cells
  ! -------------------------------------------------------------------------- !
  ncells_in_stencil = ws_stencil_size_l + ws_stencil_size_r - 1

  stencil_1d_parallel(:) = 0
  stencil_1d_perp1(:) = 0
  stencil_1d_perp2(:) = 0
  do i = 1, ws_stencil_size_l
    stencil_1d_parallel(i) = ws_stencil_l(target_df, i)
    stencil_1d_perp1(i) = ws_stencil_l(perp_df1, i)
    stencil_1d_perp2(i) = ws_stencil_l(perp_df2, i)
  end do
  ! Omit 1 from the second part, as the central cell is already included
  do i = 1, ws_stencil_size_r - 1
    stencil_1d_parallel(i+ws_stencil_size_l) = ws_stencil_r(target_df, i+1)
    stencil_1d_perp1(i+ws_stencil_size_l) = ws_stencil_r(perp_df1, i+1)
    stencil_1d_perp2(i+ws_stencil_size_l) = ws_stencil_r(perp_df2, i+1)
  end do

  ! -------------------------------------------------------------------------- !
  ! Remap
  ! -------------------------------------------------------------------------- !
  w_idx = map_ww(target_df)
  r_idx = map_wr(target_df)

  remapped_field(r_idx : r_idx+nvert) = 0.0_r_tran

  ! Loop through (multidata) interpolation points
  do i = 0, ndata - 1

    f_idx_parallel = stencil_1d_parallel(indices(w_idx+i))
    f_idx_perp1 = stencil_1d_perp1(indices(w_idx+i))
    f_idx_perp2 = stencil_1d_perp2(indices(w_idx+i))

    write(log_scratch_space,                                                   &
      '(A,I4,A,I4,A,E16.8,A,E16.8,A,E16.8,A,E16.8,A,E16.8)')                  &
      '  W2H 1D i=', i, ' idx=', indices(w_idx+i),                            &
      ' wp=', weights_parallel(w_idx+i), ' wq=', weights_perp(w_idx+i),       &
      ' fp=', field(f_idx_parallel), &
      ' fq1=', field(f_idx_perp1), ' fq2=', field(f_idx_perp2)
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

    remapped_field(r_idx : r_idx+nvert) = remapped_field(r_idx : r_idx+nvert)  &
      + weights_parallel(w_idx+i) * field(f_idx_parallel : f_idx_parallel+nvert) &
      + 0.5_r_tran * weights_perp(w_idx+i) * field(f_idx_perp1 : f_idx_perp1+nvert) &
      + 0.5_r_tran * weights_perp(w_idx+i) * field(f_idx_perp2 : f_idx_perp2+nvert)
  end do

end subroutine panel_edge_remap_w2h_1d

end module panel_edge_remap_w2h_kernel_mod
