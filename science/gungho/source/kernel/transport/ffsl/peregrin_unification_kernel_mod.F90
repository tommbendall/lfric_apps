!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------


module peregrin_unification_kernel_mod

  use argument_mod,                  only: arg_type,                           &
                                           GH_FIELD,                           &
                                           GH_READWRITE, GH_READ,              &
                                           GH_REAL, GH_INTEGER,                &
                                           ANY_DISCONTINUOUS_SPACE_1,          &
                                           ANY_DISCONTINUOUS_SPACE_2,          &
                                           STENCIL, CROSS2D,                   &
                                           CELL_COLUMN
  use constants_mod,                 only: i_def, r_def, r_tran
  use kernel_mod,                    only: kernel_type
  use reference_element_mod,         only: W, E, S, N
  use sci_face_selector_support_mod, only: face_from_face_selector

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! PSy layer
  type, public, extends(kernel_type) :: peregrin_unification_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_READWRITE,                         &
                                                   ANY_DISCONTINUOUS_SPACE_2), & ! flux_unified
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                            STENCIL(CROSS2D)), & ! flux_read
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                            STENCIL(CROSS2D)), & ! flux_edge_downwind
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), &
        arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! panel_edge_dist
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! face_selector_ew
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1)  & ! face_selector_ns
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: peregrin_unification_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: peregrin_unification_code

contains

  !> @brief Spread the flux/flux_edge_downwind discrepancy near a panel edge
  !!        across a downwind stencil
  !> @param[in]     nlayers               Number of layers
  !> @param[in,out] flux_unified          The unified flux to be updated
  !> @param[in]     flux_read             Large halo version of mass flux
  !> @param[in]     stencil_sizes_f       Sizes of branches of the cross
  !!                                      stencil for `flux`
  !> @param[in]     stencil_max_f         Maximum size of a cross stencil
  !!                                      branch for `flux`
  !> @param[in]     stencil_map_f         Dofmap for the `flux` stencil
  !> @param[in]     flux_edge_downwind    Downwind mass flux immediately next
  !!                                      to a panel edge
  !> @param[in]     stencil_sizes_d       Sizes of branches of the cross
  !!                                      stencil for `flux_edge_downwind`
  !> @param[in]     stencil_max_d         Maximum size of a cross stencil
  !!                                      branch for `flux_edge_downwind`
  !> @param[in]     stencil_map_d         Dofmap for the `flux_edge_downwind`
  !!                                      stencil
  !> @param[in]     dep_dist              Horizontal departure distances, read
  !!                                      at the column's own face to
  !!                                      determine whether that face is
  !!                                      downwind of the edge (no stencil
  !!                                      needed)
  !> @param[in]     panel_id              Panel IDs in the CROSS2D stencil
  !> @param[in]     stencil_sizes_p       Sizes of branches of the panel
  !!                                     ID stencil
  !> @param[in]     stencil_max_p         Maximum size of a panel ID
  !!                                     stencil branch
  !> @param[in]     stencil_map_p         Dofmap for the panel ID stencil
  !> @param[in]     panel_edge_dist_W     2D field: distance to the panel edge
  !!                                      to the West
  !> @param[in]     panel_edge_dist_E     2D field: distance to the panel edge
  !!                                      to the East
  !> @param[in]     panel_edge_dist_S     2D field: distance to the panel edge
  !!                                      to the South
  !> @param[in]     panel_edge_dist_N     2D field: distance to the panel edge
  !!                                      to the North
  !> @param[in]     face_selector_ew      East-West face selector, used to
  !!                                      ensure each W/E face is only ever
  !!                                      updated from one of its two
  !!                                      neighbouring columns
  !> @param[in]     face_selector_ns      North-South face selector, as
  !!                                      `face_selector_ew` for S/N faces
  !> @param[in]     ndf_w2h               Num of DoFs for W2h per cell
  !> @param[in]     undf_w2h              Num of DoFs for W2h in this partition
  !> @param[in]     map_w2h               Map for W2h
  !> @param[in]     ndf_pid               Num of DoFs for panel ID field per cell
  !> @param[in]     undf_pid              Num DoFs for this partition for panel_id
  !> @param[in]     map_pid               Map for panel ID field
  subroutine peregrin_unification_code(                                        &
      nlayers,                                                                 &
      flux_unified,                                                            &
      flux_read, stencil_sizes_f, stencil_max_f, stencil_map_f,                &
      flux_edge_downwind, stencil_sizes_d, stencil_max_d, stencil_map_d,       &
      dep_dist,                                                                &
      panel_id, stencil_sizes_p, stencil_max_p, stencil_map_p,                 &
      panel_edge_dist_W, panel_edge_dist_E,                                    &
      panel_edge_dist_S, panel_edge_dist_N,                                    &
      face_selector_ew, face_selector_ns,                                      &
      ndf_w2h, undf_w2h, map_w2h,                                              &
      ndf_pid, undf_pid, map_pid                                               &
  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_pid
    integer(kind=i_def), intent(in) :: ndf_pid
    integer(kind=i_def), intent(in) :: stencil_max_f
    integer(kind=i_def), intent(in) :: stencil_max_d
    integer(kind=i_def), intent(in) :: stencil_max_p
    integer(kind=i_def), intent(in) :: stencil_sizes_f(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_d(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_p(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in) :: stencil_map_f(ndf_w2h, stencil_max_f, 4)
    integer(kind=i_def), intent(in) :: stencil_map_d(ndf_w2h, stencil_max_d, 4)
    integer(kind=i_def), intent(in) :: stencil_map_p(ndf_pid, stencil_max_p, 4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux_unified(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux_read(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux_edge_downwind(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)

    ! Internal arguments
    integer(kind=i_def) :: edge_dist_W, edge_dist_E
    integer(kind=i_def) :: edge_dist_S, edge_dist_N
    integer(kind=i_def) :: edge_dist
    integer(kind=i_def) :: ew_sel, ns_sel
    integer(kind=i_def) :: j, face, branch_direction
    integer(kind=i_def) :: central_panel, sample_panel, panel_edge
    integer(kind=i_def) :: seam_face, seam_sign

    integer(kind=i_def), parameter :: cw_face_map(4)        = (/ S, E, N, W /)
    integer(kind=i_def), parameter :: cw_component_sign(4)  = (/ -1, 1, -1, 1 /)
    integer(kind=i_def), parameter :: acw_face_map(4)       = (/ N, W, S, E /)
    integer(kind=i_def), parameter :: acw_component_sign(4) = (/ 1, -1, 1, -1 /)

    edge_dist_W = ABS(panel_edge_dist_W(map_pid(1)))
    edge_dist_E = ABS(panel_edge_dist_E(map_pid(1)))
    edge_dist_S = ABS(panel_edge_dist_S(map_pid(1)))
    edge_dist_N = ABS(panel_edge_dist_N(map_pid(1)))

    ew_sel = face_selector_ew(map_pid(1))
    ns_sel = face_selector_ns(map_pid(1))
    central_panel = INT(panel_id(map_pid(1)), i_def)

    do j = 1, ABS(ew_sel) + ABS(ns_sel)
      face = face_from_face_selector(j, ew_sel, ns_sel)

      ! Determine distance from this face to panel edge
      ! Edge distance of zero means it lies on the panel edge itself
      select case (face)
      case (W)
        edge_dist = MIN(edge_dist_W, edge_dist_E + 1) - 1
      case (E)
        edge_dist = MIN(edge_dist_E, edge_dist_W + 1) - 1
      case (S)
        edge_dist = MIN(edge_dist_S, edge_dist_N + 1) - 1
      case (N)
        edge_dist = MIN(edge_dist_N, edge_dist_S + 1) - 1
      end select

          if (edge_dist > 0 .and. edge_dist < stencil_max_f) then
        ! Determine direction of stencil branch to use
        select case (face)
        case (W, E)
          if (edge_dist_E < edge_dist_W) then
            branch_direction = E
          else
            branch_direction = W
          end if

        case (S, N)
          if (edge_dist_S < edge_dist_N) then
            branch_direction = S
          else
            branch_direction = N
          end if
        end select

        if (edge_dist + 1 <= stencil_sizes_p(branch_direction)) then
          sample_panel = int(panel_id(stencil_map_p(                           &
            1, edge_dist + 1, branch_direction)), i_def)
          panel_edge = 10 * central_panel + sample_panel
          seam_face = face
          seam_sign = 1_i_def
          select case (panel_edge)
          case (41, 32, 16, 25, 64, 53)
            seam_face = cw_face_map(face)
            seam_sign = cw_component_sign(face)
          case (14, 23, 61, 52, 46, 35)
            seam_face = acw_face_map(face)
            seam_sign = acw_component_sign(face)
          end select

          ! If near panel edge, apply unification correction; otherwise
          ! do nothing
          call peregrin_unification_1d(                                        &
            nlayers, flux_unified, flux_read, flux_edge_downwind, dep_dist,    &
            undf_w2h, map_w2h(face), stencil_map_f, seam_face, seam_sign,      &
            edge_dist, branch_direction                                        &
          )
        end if
      end if
    end do

  end subroutine peregrin_unification_code

  subroutine peregrin_unification_1d(                                          &
      nlayers, flux_unified, flux_read, flux_edge_downwind, dep_dist,          &
      undf_w2h, idx, stencil_map, seam_face, seam_sign, edge_dist,             &
      branch_direction                                                         &
  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: stencil_map(:,:, :)
    integer(kind=i_def), intent(in)    :: edge_dist
    integer(kind=i_def), intent(in)    :: seam_face
    integer(kind=i_def), intent(in)    :: seam_sign
    integer(kind=i_def), intent(in)    :: branch_direction
    integer(kind=i_def), intent(in)    :: undf_w2h
    integer(kind=i_def), intent(in)    :: idx
    real(kind=r_tran),   intent(inout) :: flux_unified(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux_read(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux_edge_downwind(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)

    ! Internal variables
    integer(kind=i_def) :: nl, e_idx
    real(kind=r_tran)   :: branch_sign
    real(kind=r_tran)   :: seam(0:nlayers-1)
    real(kind=r_tran)   :: switch(0:nlayers-1)

    nl = nlayers - 1

    ! Set sign for displacement into the panel from the selected edge
    select case (branch_direction)
    case (W, N)
      branch_sign = 1.0_r_tran
    case (E, S)
      branch_sign = -1.0_r_tran
    end select

    ! Apply correction only where dep_dist+1 exceeds distance to panel edge
    switch(:) = 0.5_r_tran + sign(0.5_r_tran,                                  &
        branch_sign * dep_dist(idx:idx + nl)                                   &
        - real(edge_dist - 1, r_tran)                                          &
    )

    ! Calculate seam
    e_idx = stencil_map(seam_face, edge_dist + 1, branch_direction)
    seam(:) = real(seam_sign, r_tran) * (                                      &
      flux_read(e_idx:e_idx + nl) - flux_edge_downwind(e_idx:e_idx + nl))

    ! Apply correction, which decays away from the panel edge
    flux_unified(idx:idx+nl) = flux_unified(idx:idx+nl)                        &
        + seam(:) * switch(:)                                                  &
        * (1.0_r_tran - REAL(edge_dist, r_tran) / (REAL(edge_dist+1, r_tran)))

  end subroutine peregrin_unification_1d

end module peregrin_unification_kernel_mod
