!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Remap a W2H flux/wind field at the edges of cubed-sphere panels to
!!        the corresponding neighbouring panels, using the staged
!!        interpolation values and rotation factors computed by
!!        panel_edge_remap_w2h_weights_kernel_mod and
!!        panel_edge_interpolate_w2h_kernel_mod.
!!
!> @details This kernel applies the final weighted sum (see the instructions
!!          doc, Sec. 5) for every active panel label at each selected face:
!!             remapped(label) = (parl_rot(label)*parl(label)
!!                                + perp_rot(label)*perp(label))*dA_eave
!!
!!          The labelled staged values come from
!!          panel_edge_interpolate_w2h_kernel_mod; the labelled rotation
!!          factors come from panel_edge_remap_w2h_weights_kernel_mod.
!!
!!          Away from the relevant panel edge, parl_w2_rot defaults to 1 and
!!          perp_w2_rot defaults to 0, and parl_w2_at_w2/perp_w2_at_w2 default
!!          to the input field's own value there, so the formula above
!!          reduces to the identity remap with no special-case branch
!!          required. Do not reintroduce copy branches here.
module panel_edge_remap_w2h_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER, GH_LOGICAL,              &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_1,                    &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_6,                    &
                                 CELL_COLUMN
use fs_continuity_mod,      only: W2h, W2
use constants_mod,          only: r_def, r_tran, i_def, l_def
use reference_element_mod,  only: W, S, N, E
use sci_face_selector_support_mod, only: face_from_face_selector

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_edge_remap_w2h_kernel_type
  private
    type(arg_type) :: meta_args(15) = (/                                        &
       arg_type(GH_FIELD*3, GH_REAL, GH_WRITE,                                 &
                ANY_DISCONTINUOUS_SPACE_1),                                    &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  W2h),                        &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  W2h),                        &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_6),  &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_6),  &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
      arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
      arg_type(GH_FIELD*3, GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2),                         &
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                              &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                               &
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
!> @param[in,out] remapped_flux_1      Remapped field for eaves
!!                                      whose neighbouring panels are 1 or 3
!> @param[in,out] remapped_flux_2      As remapped_flux_1, for panels
!!                                      2 or 4
!> @param[in,out] remapped_flux_3      As remapped_flux_1, for panels
!!                                      5 or 6
!> @param[in]     parl_w2_at_w2_1       "parl_w2" staged value for label 1
!> @param[in]     parl_w2_at_w2_2       "parl_w2" staged value for label 2
!> @param[in]     parl_w2_at_w2_3       "parl_w2" staged value for label 3
!> @param[in]     perp_w2_at_w2_1       "perp_w2" staged value for label 1
!> @param[in]     perp_w2_at_w2_2       "perp_w2" staged value for label 2
!> @param[in]     perp_w2_at_w2_3       "perp_w2" staged value for label 3
!> @param[in]     parl_w2_rot_1         Parallel rotation factor for label 1
!> @param[in]     parl_w2_rot_2         Parallel rotation factor for label 2
!> @param[in]     parl_w2_rot_3         Parallel rotation factor for label 3
!> @param[in]     perp_w2_rot_1         Perpendicular rotation factor for label 1
!> @param[in]     perp_w2_rot_2         Perpendicular rotation factor for label 2
!> @param[in]     perp_w2_rot_3         Perpendicular rotation factor for label 3
!> @param[in]     eave_lon_1            Diagnostic eave longitude for label 1
!> @param[in]     eave_lon_2            Diagnostic eave longitude for label 2
!> @param[in]     eave_lon_3            Diagnostic eave longitude for label 3
!> @param[in]     eave_lat_1            Diagnostic eave latitude for label 1
!> @param[in]     eave_lat_2            Diagnostic eave latitude for label 2
!> @param[in]     eave_lat_3            Diagnostic eave latitude for label 3
!> @param[in]     correct_flux_1        Diagnostic analytic flux for label 1
!> @param[in]     correct_flux_2        Diagnostic analytic flux for label 2
!> @param[in]     correct_flux_3        Diagnostic analytic flux for label 3
!> @param[in]     interp_w2h_indices_1  Signed donor offsets; zero means inactive
!> @param[in]     interp_w2h_indices_2  Signed donor offsets; zero means inactive
!> @param[in]     interp_w2h_indices_3  Signed donor offsets; zero means inactive
!> @param[in]     dA_at_w2_eave         W2 face-area field, used to scale the
!!                                      remapped output back into the target
!!                                      eave's area (area varies per layer,
!!                                      unlike the rotation factors which are 2D)
!> @param[in]     face_selector_ew      East-West face selector, used to ensure
!!                                      each W/E face is only ever written by
!!                                      one of its two neighbouring columns
!> @param[in]     face_selector_ns      North-South face selector, as
!!                                      face_selector_ew for S/N faces
!> @param[in]     ndf_wr                Num DoFs per cell for remapped fields
!> @param[in]     undf_wr               Num DoFs for this partition for remapped
!!                                      fields
!> @param[in]     map_wr                DoF map for remapped fields
!> @param[in]     ndf_w2h               Num DoFs per cell for the staged values
!> @param[in]     undf_w2h              Num DoFs for this partition for the
!!                                      staged values
!> @param[in]     map_w2h               DoF map for the staged values
!> @param[in]     ndf_wrot              Num DoFs per cell for the (2D) rotation
!!                                      factors
!> @param[in]     undf_wrot             Num DoFs for this partition for the
!!                                      rotation factors
!> @param[in]     map_wrot              DoF map for the rotation factors
!> @param[in]     ndf_ww                Num DoFs per cell for diagnostic and mask fields
!> @param[in]     undf_ww               Num DoFs for this partition for the
!!                                      diagnostic fields
!> @param[in]     map_ww                DoF map for the diagnostic fields
!> @param[in]     ndf_w2                Num DoFs per cell for dA_at_w2_eave
!> @param[in]     undf_w2               Num DoFs for this partition for dA_at_w2_eave
!> @param[in]     map_w2                DoF map for dA_at_w2_eave
!> @param[in]     ndf_pid               Num DoFs per cell for panel ID field
!> @param[in]     undf_pid              Num DoFs for this partition for panel ID
!> @param[in]     map_pid               DoF map for panel ID field
subroutine panel_edge_remap_w2h_code( nlayers,                                 &
                                      remapped_flux_1,                         &
                                      remapped_flux_2,                         &
                                      remapped_flux_3,                         &
                                      parl_w2_at_w2_1,                        &
                                      parl_w2_at_w2_2,                        &
                                      parl_w2_at_w2_3,                        &
                                      perp_w2_at_w2_1,                        &
                                      perp_w2_at_w2_2,                        &
                                      perp_w2_at_w2_3,                        &
                                      parl_w2_rot_1,                          &
                                      parl_w2_rot_2,                          &
                                      parl_w2_rot_3,                          &
                                      perp_w2_rot_1,                          &
                                      perp_w2_rot_2,                          &
                                      perp_w2_rot_3,                          &
                                      eave_lon_1, eave_lon_2, eave_lon_3,     &
                                      eave_lat_1, eave_lat_2, eave_lat_3,     &
                                      correct_flux_1, correct_flux_2,         &
                                      correct_flux_3,                         &
                                      interp_w2h_indices_1,                    &
                                      interp_w2h_indices_2,                    &
                                      interp_w2h_indices_3,                    &
                                      dA_at_w2_eave,                         &
                                      panel_edge_dist_W, panel_edge_dist_E,  &
                                      panel_edge_dist_S, panel_edge_dist_N,  &
                                      face_selector_ew,                      &
                                      face_selector_ns,                      &
                                      remap_depth, include_perp_p,           &
                                      ndf_wr, undf_wr, map_wr,               &
                                      ndf_w2h, undf_w2h, map_w2h,            &
                                      ndf_wrot, undf_wrot, map_wrot,         &
                                      ndf_ww, undf_ww, map_ww,               &
                                      ndf_w2, undf_w2, map_w2,               &
                                      ndf_pid, undf_pid, map_pid             &
  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_wrot, undf_wrot
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_wrot(ndf_wrot)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)

  integer(kind=i_def), intent(in)    :: remap_depth
  logical(kind=l_def), intent(in)    :: include_perp_p

  real(kind=r_tran),   intent(inout) :: remapped_flux_1(undf_wr)
  real(kind=r_tran),   intent(inout) :: remapped_flux_2(undf_wr)
  real(kind=r_tran),   intent(inout) :: remapped_flux_3(undf_wr)
  real(kind=r_tran),   intent(in)    :: parl_w2_at_w2_1(undf_w2h)
  real(kind=r_tran),   intent(in)    :: parl_w2_at_w2_2(undf_w2h)
  real(kind=r_tran),   intent(in)    :: parl_w2_at_w2_3(undf_w2h)
  real(kind=r_tran),   intent(in)    :: perp_w2_at_w2_1(undf_w2h)
  real(kind=r_tran),   intent(in)    :: perp_w2_at_w2_2(undf_w2h)
  real(kind=r_tran),   intent(in)    :: perp_w2_at_w2_3(undf_w2h)
  real(kind=r_tran),   intent(in)    :: parl_w2_rot_1(undf_wrot)
  real(kind=r_tran),   intent(in)    :: parl_w2_rot_2(undf_wrot)
  real(kind=r_tran),   intent(in)    :: parl_w2_rot_3(undf_wrot)
  real(kind=r_tran),   intent(in)    :: perp_w2_rot_1(undf_wrot)
  real(kind=r_tran),   intent(in)    :: perp_w2_rot_2(undf_wrot)
  real(kind=r_tran),   intent(in)    :: perp_w2_rot_3(undf_wrot)
  real(kind=r_tran),   intent(in)    :: eave_lon_1(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lon_2(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lon_3(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lat_1(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lat_2(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lat_3(undf_ww)
  real(kind=r_tran),   intent(in)    :: correct_flux_1(undf_ww)
  real(kind=r_tran),   intent(in)    :: correct_flux_2(undf_ww)
  real(kind=r_tran),   intent(in)    :: correct_flux_3(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_1(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_2(undf_ww)
  integer(kind=i_def), intent(in)    :: interp_w2h_indices_3(undf_ww)
  real(kind=r_def),    intent(in)    :: dA_at_w2_eave(undf_w2)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)

  ! Internal variables
  logical(kind=l_def) :: near_panel_edge_x, near_panel_edge_y, near_panel_edge
  integer(kind=i_def) :: face, pid_idx
  integer(kind=i_def) :: ew_sel, ns_sel, j, label
  logical(kind=l_def) :: face_active

  pid_idx = map_pid(1)

  ew_sel = face_selector_ew(pid_idx)
  ns_sel = face_selector_ns(pid_idx)

  near_panel_edge_x = (                                                        &
      ABS(panel_edge_dist_W(pid_idx)) <= remap_depth                           &
      .or. ABS(panel_edge_dist_E(pid_idx)) <= remap_depth                      &
  )
  near_panel_edge_y = (                                                        &
      ABS(panel_edge_dist_S(pid_idx)) <= remap_depth                           &
      .or. ABS(panel_edge_dist_N(pid_idx)) <= remap_depth                      &
  )
  near_panel_edge = (near_panel_edge_x .or. near_panel_edge_y)

  ! Loop only over the faces this column owns (per face_selector_ew/ns), so
  ! that each shared W2H face is written by exactly one of its two
  ! neighbouring columns
  do j = 1, ABS(ew_sel) + ABS(ns_sel)
    face = face_from_face_selector(j, ew_sel, ns_sel)

    if (.not. near_panel_edge) then
      ! Missing data value
      remapped_flux_1(map_wr(face)) = -1900.0_r_tran
      remapped_flux_2(map_wr(face)) = -1900.0_r_tran
      remapped_flux_3(map_wr(face)) = -1900.0_r_tran
      cycle
    end if

    select case (face)
    case (W, E)
      face_active = near_panel_edge_y .or.                                  &
                    (near_panel_edge_x .and. include_perp_p)
    case (S, N)
      face_active = near_panel_edge_x .or.                                  &
                    (near_panel_edge_y .and. include_perp_p)
    end select
    if (.not. face_active) cycle

    do label = 1, 3
      select case (label)
      case (1)
        if (interp_w2h_indices_1(map_ww(face)) == 0_i_def) cycle
        call remap_label(label, map_wr(face), map_w2h(face), map_wrot(face),  &
            map_w2(face), parl_w2_at_w2_1, perp_w2_at_w2_1,                  &
            parl_w2_rot_1, perp_w2_rot_1)
      case (2)
        if (interp_w2h_indices_2(map_ww(face)) == 0_i_def) cycle
        call remap_label(label, map_wr(face), map_w2h(face), map_wrot(face),  &
            map_w2(face), parl_w2_at_w2_2, perp_w2_at_w2_2,                  &
            parl_w2_rot_2, perp_w2_rot_2)
      case (3)
        if (interp_w2h_indices_3(map_ww(face)) == 0_i_def) cycle
        call remap_label(label, map_wr(face), map_w2h(face), map_wrot(face),  &
            map_w2(face), parl_w2_at_w2_3, perp_w2_at_w2_3,                  &
            parl_w2_rot_3, perp_w2_rot_3)
      end select

    end do
  end do

contains

  !> @brief Store the first nlayers values of the remapped field with axis
  !!        label "label", starting at index "r_idx"
  subroutine store_remapped(label, r_idx, values)
    implicit none
    integer(kind=i_def), intent(in) :: label, r_idx
    real(kind=r_tran),   intent(in) :: values(nlayers)

    select case (label)
    case (1)
      remapped_flux_1(r_idx : r_idx + nlayers - 1) = values
    case (2)
      remapped_flux_2(r_idx : r_idx + nlayers - 1) = values
    case (3)
      remapped_flux_3(r_idx : r_idx + nlayers - 1) = values
    end select

  end subroutine store_remapped

  subroutine remap_label(label, r_idx, w2h_idx, rot_idx, w2_idx, parl, perp, &
                         parl_rot, perp_rot)

    implicit none

    integer(kind=i_def), intent(in) :: label, r_idx, w2h_idx, rot_idx, w2_idx
    real(kind=r_tran), intent(in) :: parl(undf_w2h), perp(undf_w2h)
    real(kind=r_tran), intent(in) :: parl_rot(undf_wrot), perp_rot(undf_wrot)

    real(kind=r_tran) :: values(nlayers)

    values = (parl_rot(rot_idx) * parl(w2h_idx:w2h_idx+nlayers-1)             &
        + perp_rot(rot_idx) * perp(w2h_idx:w2h_idx+nlayers-1))                 &
        * dA_at_w2_eave(w2_idx:w2_idx+nlayers-1)
    call store_remapped(label, r_idx, values)

  end subroutine remap_label

end subroutine panel_edge_remap_w2h_code

end module panel_edge_remap_w2h_kernel_mod
