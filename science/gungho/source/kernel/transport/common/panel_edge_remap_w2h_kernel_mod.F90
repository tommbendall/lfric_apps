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
!!          doc, Sec. 5) at every DoF this column owns:
!!             remapped = parl_w2_rot * parl_w2_at_w2 + perp_w2_rot * perp_w2_at_w2
!!             remapped = remapped * dA_at_w2_eave
!!          for both the x (S/N target) and y (W/E target) remaps.
!!
!!          parl_w2_at_w2_{x,y} is the pure Lagrange interpolation already
!!          computed by panel_edge_interpolate_w2h_kernel_mod; perp_w2_at_w2_{x,y}
!!          is the staged value from the (not yet implemented) W1V/W3 average
!!          kernel. parl_w2_rot_{x,y} / perp_w2_rot_{x,y} are the rotation
!!          factors computed by panel_edge_remap_w2h_weights_kernel_mod.
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
! Diagnostic logging retained deliberately for this phase of the rewrite (see
! the instructions doc, Step 4.2) - gate or remove once the rewrite is fully
! validated.
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_DEBUG

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_edge_remap_w2h_kernel_type
  private
  type(arg_type) :: meta_args(22) = (/                                        &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! remapped_flux_x
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! remapped_flux_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h),                       & ! parl_w2_at_w2_x
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h),                       & ! parl_w2_at_w2_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h),                       & ! perp_w2_at_w2_x
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2h),                       & ! perp_w2_at_w2_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_6), & ! parl_w2_rot_x
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_6), & ! perp_w2_rot_x
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_6), & ! parl_w2_rot_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_6), & ! perp_w2_rot_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! eave_lon_x (diagnostic)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! eave_lat_x (diagnostic)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! correct_flux_x (diagnostic)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! eave_lon_y (diagnostic)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! eave_lat_y (diagnostic)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! correct_flux_y (diagnostic)
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2),                       & ! dA_at_w2_eave
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_edge_dist_W/E/S/N
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ew
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ns
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! remap_depth
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                              & ! include_perp_p
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
!> @param[in]     parl_w2_at_w2_x       "parl_w2" staged value for the x remap
!> @param[in]     parl_w2_at_w2_y       "parl_w2" staged value for the y remap
!> @param[in]     perp_w2_at_w2_x       "perp_w2" staged value for the x remap
!> @param[in]     perp_w2_at_w2_y       "perp_w2" staged value for the y remap
!> @param[in]     parl_w2_rot_x         "parl_w2" rotation factor for the x remap
!> @param[in]     perp_w2_rot_x         "perp_w2" rotation factor for the x remap
!> @param[in]     parl_w2_rot_y         "parl_w2" rotation factor for the y remap
!> @param[in]     perp_w2_rot_y         "perp_w2" rotation factor for the y remap
!> @param[in]     eave_lon_x            Diagnostic: true eave longitude (x remap)
!> @param[in]     eave_lat_x            Diagnostic: true eave latitude (x remap)
!> @param[in]     correct_flux_x        Diagnostic: analytic "correct" flux (x remap)
!> @param[in]     eave_lon_y            Diagnostic: true eave longitude (y remap)
!> @param[in]     eave_lat_y            Diagnostic: true eave latitude (y remap)
!> @param[in]     correct_flux_y        Diagnostic: analytic "correct" flux (y remap)
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
!> @param[in]     ndf_ww                Num DoFs per cell for the diagnostic fields
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
                                      remapped_flux_x,                        &
                                      remapped_flux_y,                        &
                                      parl_w2_at_w2_x,                        &
                                      parl_w2_at_w2_y,                        &
                                      perp_w2_at_w2_x,                        &
                                      perp_w2_at_w2_y,                        &
                                      parl_w2_rot_x,                          &
                                      perp_w2_rot_x,                         &
                                      parl_w2_rot_y,                          &
                                      perp_w2_rot_y,                         &
                                      eave_lon_x, eave_lat_x, correct_flux_x, &
                                      eave_lon_y, eave_lat_y, correct_flux_y, &
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

  real(kind=r_tran),   intent(inout) :: remapped_flux_x(undf_wr)
  real(kind=r_tran),   intent(inout) :: remapped_flux_y(undf_wr)
  real(kind=r_tran),   intent(in)    :: parl_w2_at_w2_x(undf_w2h)
  real(kind=r_tran),   intent(in)    :: parl_w2_at_w2_y(undf_w2h)
  real(kind=r_tran),   intent(in)    :: perp_w2_at_w2_x(undf_w2h)
  real(kind=r_tran),   intent(in)    :: perp_w2_at_w2_y(undf_w2h)
  real(kind=r_tran),   intent(in)    :: parl_w2_rot_x(undf_wrot)
  real(kind=r_tran),   intent(in)    :: perp_w2_rot_x(undf_wrot)
  real(kind=r_tran),   intent(in)    :: parl_w2_rot_y(undf_wrot)
  real(kind=r_tran),   intent(in)    :: perp_w2_rot_y(undf_wrot)
  real(kind=r_tran),   intent(in)    :: eave_lon_x(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lat_x(undf_ww)
  real(kind=r_tran),   intent(in)    :: correct_flux_x(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lon_y(undf_ww)
  real(kind=r_tran),   intent(in)    :: eave_lat_y(undf_ww)
  real(kind=r_tran),   intent(in)    :: correct_flux_y(undf_ww)
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
  integer(kind=i_def) :: ew_sel, ns_sel, j
  integer(kind=i_def) :: w_idx, wr_idx
  real(kind=r_tran)   :: diff, pct_diff

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

    if (near_panel_edge) then

      select case (face)
      case (W, E)
        if (near_panel_edge_y .or. (near_panel_edge_x .and. include_perp_p)) then
          call panel_edge_remap_w2h_1d(                                        &
                  nlayers, remapped_flux_y,                                    &
                  parl_w2_at_w2_y, perp_w2_at_w2_y,                            &
                  parl_w2_rot_y, perp_w2_rot_y,                                &
                  dA_at_w2_eave, face,                                         &
                  ndf_wr, undf_wr, map_wr,                                     &
                  ndf_w2h, undf_w2h, map_w2h,                                  &
                  ndf_wrot, undf_wrot, map_wrot,                               &
                  ndf_w2, undf_w2, map_w2                                      &
          )


        end if

      case (S, N)
        if (near_panel_edge_x .or. (near_panel_edge_y .and. include_perp_p)) then
          call panel_edge_remap_w2h_1d(                                        &
                  nlayers, remapped_flux_x,                                    &
                  parl_w2_at_w2_x, perp_w2_at_w2_x,                            &
                  parl_w2_rot_x, perp_w2_rot_x,                                &
                  dA_at_w2_eave, face,                                         &
                  ndf_wr, undf_wr, map_wr,                                     &
                  ndf_w2h, undf_w2h, map_w2h,                                  &
                  ndf_wrot, undf_wrot, map_wrot,                               &
                  ndf_w2, undf_w2, map_w2                                      &
          )
        end if

      end select

      ! ====================================================================== !
      ! Diagnostic, to remove
      ! ====================================================================== !
      if ((face == W .or. face == E) .and.                                     &
           (near_panel_edge_y .or. (near_panel_edge_x .and. include_perp_p))   &
           .or.                                                                &
           ((face == N .or. face == S) .and.                                   &
            (near_panel_edge_x .or. (near_panel_edge_y .and. include_perp_p)))) then

        write(log_scratch_space, '(A,I0,A,I0,A,E20.12,A,E20.12)')                  &
          'PEREGRIN W2H PRE  col=', pid_idx, ' face=', face,                      &
          ' parl_w2_at_w2_x=', parl_w2_at_w2_x(map_w2h(face)),                    &
          ' parl_w2_at_w2_y=', parl_w2_at_w2_y(map_w2h(face))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

        wr_idx = map_wr(face)
        ! TEMP DEBUG T4: to remove
        write(log_scratch_space, *) 'T4 KR col=', pid_idx, ' face=', face,    &
          ' dist_WESN=', panel_edge_dist_W(pid_idx), panel_edge_dist_E(pid_idx),&
          panel_edge_dist_S(pid_idx), panel_edge_dist_N(pid_idx),              &
          ' remap_depth=', remap_depth,                                        &
          ' parl_rot_x=', parl_w2_rot_x(map_wrot(face)),                       &
          ' perp_rot_x=', perp_w2_rot_x(map_wrot(face)),                       &
          ' perp_at_x=', perp_w2_at_w2_x(map_w2h(face)),                       &
          ' parl_rot_y=', parl_w2_rot_y(map_wrot(face)),                       &
          ' perp_rot_y=', perp_w2_rot_y(map_wrot(face)),                       &
          ' perp_at_y=', perp_w2_at_w2_y(map_w2h(face)),                       &
          ' dA_eave=', dA_at_w2_eave(map_w2(face))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
        write(log_scratch_space, '(A,I0,A,I0,A,E20.12,A,E20.12)')                  &
          'PEREGRIN W2H POST col=', pid_idx, ' face=', face,                      &
          ' remapped_flux_x=', remapped_flux_x(wr_idx),                          &
          ' remapped_flux_y=', remapped_flux_y(wr_idx)
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

        ! Diagnostic: compare against the analytic solid-body-rotation "correct"
        ! flux at this eave, for both directions independently - each is only
        ! ever non-identity near its own relevant panel edge (near_panel_edge_x
        ! gates the x remap, near_panel_edge_y the y remap)
        w_idx = map_ww(face)

        diff = remapped_flux_x(wr_idx) - correct_flux_x(w_idx)
        if (ABS(correct_flux_x(w_idx)) > EPSILON(0.0_r_tran)) then
          pct_diff = 100.0_r_tran * diff / correct_flux_x(w_idx)
        else
          pct_diff = 0.0_r_tran
        end if
        write(log_scratch_space,                                                   &
          '(A,I0,A,I0,A,E16.8,A,E16.8,A,E16.8,A,E16.8,A,E16.8)')                  &
          'PEREGRIN W2H DIAG col=', pid_idx, ' face=', face,                       &
          ' eave_lon_x=', eave_lon_x(w_idx), ' eave_lat_x=', eave_lat_x(w_idx),   &
          ' correct_x=', correct_flux_x(w_idx), ' diff_x=', diff,                 &
          ' pct_diff_x=', pct_diff
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

        diff = remapped_flux_y(wr_idx) - correct_flux_y(w_idx)
        if (ABS(correct_flux_y(w_idx)) > EPSILON(0.0_r_tran)) then
          pct_diff = 100.0_r_tran * diff / correct_flux_y(w_idx)
        else
          pct_diff = 0.0_r_tran
        end if
        write(log_scratch_space,                                                   &
          '(A,I0,A,I0,A,E16.8,A,E16.8,A,E16.8,A,E16.8,A,E16.8)')                  &
          'PEREGRIN W2H DIAG col=', pid_idx, ' face=', face,                       &
          ' eave_lon_y=', eave_lon_y(w_idx), ' eave_lat_y=', eave_lat_y(w_idx),   &
          ' correct_y=', correct_flux_y(w_idx), ' diff_y=', diff,                 &
          ' pct_diff_y=', pct_diff
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
      end if

    else
      remapped_flux_x(map_w2h(face)) = -1900.0_r_tran  ! Missing data value
    end if
  end do

end subroutine panel_edge_remap_w2h_code

!> @brief Private routine to remap a single W2H target DoF ("face") by
!!        combining its "parl_w2" and "perp_w2" staged values with the
!!        rotation factors, then scaling by dA_at_w2_eave. No stencil
!!        handling here - that is done by panel_edge_interpolate_w2h_kernel_mod
!!        and (in later phases) panel_edge_w2h_average_kernel_mod.
subroutine panel_edge_remap_w2h_1d( nlayers,                                   &
                                    remapped_field,                           &
                                    parl_w2_at_w2, perp_w2_at_w2,              &
                                    parl_w2_rot, perp_w2_rot,                  &
                                    dA_at_w2_eave, face,                       &
                                    ndf_wr, undf_wr, map_wr,                   &
                                    ndf_w2h, undf_w2h, map_w2h,                &
                                    ndf_wrot, undf_wrot, map_wrot,             &
                                    ndf_w2, undf_w2, map_w2                    &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: face
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_wrot, undf_wrot
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_wrot(ndf_wrot)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)

  real(kind=r_tran),   intent(inout) :: remapped_field(undf_wr)
  real(kind=r_tran),   intent(in)    :: parl_w2_at_w2(undf_w2h)
  real(kind=r_tran),   intent(in)    :: perp_w2_at_w2(undf_w2h)
  real(kind=r_tran),   intent(in)    :: parl_w2_rot(undf_wrot)
  real(kind=r_tran),   intent(in)    :: perp_w2_rot(undf_wrot)
  real(kind=r_def),    intent(in)    :: dA_at_w2_eave(undf_w2)

  ! Internal variables
  integer(kind=i_def) :: nvert, r_idx, w2h_idx, rot_idx, w2_idx

  nvert = nlayers - 1
  r_idx = map_wr(face)
  w2h_idx = map_w2h(face)
  rot_idx = map_wrot(face)
  w2_idx = map_w2(face)

  ! The rotation factors are 2D (one value per column), so are not indexed
  ! over the vertical
  remapped_field(r_idx : r_idx+nvert) =                                        &
      ( parl_w2_rot(rot_idx) * parl_w2_at_w2(w2h_idx : w2h_idx+nvert)          &
      + perp_w2_rot(rot_idx) * perp_w2_at_w2(w2h_idx : w2h_idx+nvert) )        &
      * dA_at_w2_eave(w2_idx : w2_idx+nvert)

end subroutine panel_edge_remap_w2h_1d

end module panel_edge_remap_w2h_kernel_mod
