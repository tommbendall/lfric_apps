!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Unifies `flux` and `flux_edge_downwind` near a cubed-sphere panel
!!          edge, by spreading the mass "gap" between them out over a
!!          downwind stencil, sized by the Courant number at the edge.
!> @details Near a panel edge, `peregrin_flux_kernel_mod.F90` may leave a
!!          genuine discrepancy between `flux` (the value actually used,
!!          which for an upwind-classified column is the panel-remap value)
!!          and `flux_edge_downwind` (the sphere-based value that a fully
!!          downwind-consistent calculation would have given at that same
!!          face). This kernel computes that discrepancy ("seam"):
!!
!!              seam = flux(edge face) - flux_edge_downwind(edge face)
!!
!!          at the immediate edge-adjacent face, and spreads it, using a
!!          linear ramp, across the `n_cells = ceiling(courant + 1)` cells
!!          downwind of the edge (`courant = ABS(dep_dist)`, read at the face
!!          being calculated, so the cap on the redistribution depth comes
!!          directly from `dep_dist` rather than from a separately-supplied
!!          scalar, and no stencil is needed on `dep_dist`) - the
!!          minimum-norm, non-monotone redistribution described in tracker
!!          step 3.4.4 (mathematically equivalent to
!!          `unify_flux_x`/`_unify_flux_row` in the separate swift_ml
!!          repository, but restructured for the LFRic per-column,
!!          stencil-read-only kernel model).
!!
!!          Because an LFRic kernel can only write to its own column's dofs
!!          (never to a stencil neighbour's), the calculation is inverted
!!          relative to swift_ml: rather than looping outward from the edge
!!          and writing to each downwind facet, every column near an edge
!!          works out its own distance `k` from `panel_edge_dist`, reads
!!          `flux`/`flux_edge_downwind`/`dep_dist` at the immediate
!!          edge-adjacent column via an X1D/Y1D stencil (built here from a
!!          `STENCIL(CROSS2D)` access, as in `peregrin_flux_kernel_mod.F90`),
!!          and - only if it is on the downwind side of that edge - adds its
!!          own share of the ramp to its own local edge-facing face of
!!          `flux_correction`.
!!
!!          Since `flux_correction` (like `flux`) lives on a continuous W2h
!!          space where each face dof is shared between two columns, the
!!          `face_selector_ew`/`face_selector_ns` fields (as used elsewhere in
!!          the codebase, e.g. `ffsl_flux_xy_kernel_mod.F90`) are used to
!!          ensure that each face is only ever updated from one of its two
!!          neighbouring columns, avoiding a double-write/race.
!!
!!          `flux_correction` must be initialised to zero before this kernel
!!          is called, and is intended to be added onto `flux` by the calling
!!          algorithm afterwards (e.g. via the `X_plus_Y` builtin), since a
!!          column may add a contribution to `flux_correction` for more than
!!          one face (if it is within reach of more than one panel edge).
!!
!> @note This is a first draft (tracker step 3.4.4). The following points in
!!       particular have not yet been validated against the PEREGRIN design
!!       and should be reviewed:
!!       - Whether `courant = ABS(dep_dist)` (in units of grid cells) is the
!!         correct quantity to size the redistribution stencil, or whether a
!!         separate wind/dry-mass-based Courant number is needed instead.
!!       - Whether reading `dep_dist` at the column's own face (rather than
!!         at the edge-adjacent reference face) to decide downwind-ness and
!!         the Courant number is the correct choice.
!!       - Whether `peregrin_flux_kernel_mod.F90` can leave a genuine race
!!         between the two columns either side of a shared edge face each
!!         writing their own `flux_edge_downwind` there (both are "immediate"
!!         to the same edge from their own side) - this has not been checked.
!> @note This kernel only works when field is a W3 field at lowest order,
!!       and flux/flux_edge_downwind/dep_dist are on the lowest order W2h
!!       dof layout (1=W, 2=S, 3=E, 4=N).

module peregrin_unification_kernel_mod

  use argument_mod,          only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   CELL_COLUMN, GH_READWRITE, &
                                   GH_READ,                    &
                                   GH_INTEGER,                &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   STENCIL, CROSS2D
  use constants_mod,                 only: i_def, r_tran
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
    type(arg_type) :: meta_args(7) = (/                                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_READWRITE,                         &
                                                   ANY_DISCONTINUOUS_SPACE_2), & ! flux_correction
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                            STENCIL(CROSS2D)), & ! flux
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                            STENCIL(CROSS2D)), & ! flux_edge_downwind
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
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
  !!        across a downwind stencil, sized by the Courant number there.
  !> @param[in]     nlayers               Number of layers
  !> @param[in,out] flux_correction       Correction to be added onto `flux`
  !!                                      by the calling algorithm; must be
  !!                                      initialised to zero beforehand
  !> @param[in]     flux                  The mass flux
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
  subroutine peregrin_unification_code( nlayers,             &
                                        flux_correction,     &
                                        flux,                &
                                        stencil_sizes_f,     &
                                        stencil_max_f,       &
                                        stencil_map_f,       &
                                        flux_edge_downwind,  &
                                        stencil_sizes_d,     &
                                        stencil_max_d,       &
                                        stencil_map_d,       &
                                        dep_dist,            &
                                        panel_edge_dist_W,   &
                                        panel_edge_dist_E,   &
                                        panel_edge_dist_S,   &
                                        panel_edge_dist_N,   &
                                        face_selector_ew,    &
                                        face_selector_ns,    &
                                        ndf_w2h,             &
                                        undf_w2h,            &
                                        map_w2h,             &
                                        ndf_pid,             &
                                        undf_pid,            &
                                        map_pid )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_pid
    integer(kind=i_def), intent(in) :: ndf_pid
    integer(kind=i_def), intent(in) :: stencil_max_f
    integer(kind=i_def), intent(in) :: stencil_max_d
    integer(kind=i_def), intent(in) :: stencil_sizes_f(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_d(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in) :: stencil_map_f(ndf_w2h, stencil_max_f, 4)
    integer(kind=i_def), intent(in) :: stencil_map_d(ndf_w2h, stencil_max_d, 4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux_correction(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux_edge_downwind(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)

    ! Internal arguments
    integer(kind=i_def) :: edge_dist_W, edge_dist_E, edge_dist_S, edge_dist_N
    integer(kind=i_def) :: ew_sel, ns_sel
    integer(kind=i_def) :: j, face

    edge_dist_W = panel_edge_dist_W(map_pid(1))
    edge_dist_E = panel_edge_dist_E(map_pid(1))
    edge_dist_S = panel_edge_dist_S(map_pid(1))
    edge_dist_N = panel_edge_dist_N(map_pid(1))

    ew_sel = face_selector_ew(map_pid(1))
    ns_sel = face_selector_ns(map_pid(1))

    ! Only iterate over the faces this column owns (per `face_selector_ew`/
    ! `face_selector_ns`), so that each shared face is updated from exactly
    ! one of its two neighbouring columns.
    ! W and S are the "direction = +1" sides (correction grows away from the
    ! edge in the direction of increasing stencil branch); E and N are the
    ! "direction = -1" sides.
    do j = 1, ABS(ew_sel) + ABS(ns_sel)
      face = face_from_face_selector(j, ew_sel, ns_sel)

      select case(face)
        case(W)
          call apply_unification_side(                                        &
              flux_correction, flux, flux_edge_downwind, dep_dist,            &
              map_w2h(W), stencil_map_f(W,:,W), edge_dist_W, 1_i_def,         &
              stencil_max_f, nlayers, undf_w2h                                &
          )
        case(E)
          call apply_unification_side(                                        &
              flux_correction, flux, flux_edge_downwind, dep_dist,            &
              map_w2h(E), stencil_map_f(E,:,E), edge_dist_E, -1_i_def,        &
              stencil_max_f, nlayers, undf_w2h                                &
          )
        case(S)
          call apply_unification_side(                                        &
              flux_correction, flux, flux_edge_downwind, dep_dist,            &
              map_w2h(S), stencil_map_f(S,:,S), edge_dist_S, 1_i_def,         &
              stencil_max_f, nlayers, undf_w2h                                &
          )
        case(N)
          call apply_unification_side(                                        &
              flux_correction, flux, flux_edge_downwind, dep_dist,            &
              map_w2h(N), stencil_map_f(N,:,N), edge_dist_N, -1_i_def,        &
              stencil_max_f, nlayers, undf_w2h                                &
          )
      end select
    end do

  end subroutine peregrin_unification_code

  !> @brief Apply the downwind ramp correction for a single side (W, E, S or
  !!        N) of a column, if it is downwind of a panel edge on that side
  !!        within reach of the stencil.
  !> @param[in,out] flux_correction     Correction to add the ramp term to
  !> @param[in]     flux                The mass flux
  !> @param[in]     flux_edge_downwind  Downwind mass flux next to a panel edge
  !> @param[in]     dep_dist            Horizontal departure distances, read
  !!                                    at this column's own face (`own_idx`)
  !> @param[in]     own_idx             Index of this column's own edge-facing
  !!                                    dof in `flux_correction`
  !> @param[in]     stencil_branch_f    Stencil dofmap for `flux`, along the
  !!                                    branch towards the edge, for this side's
  !!                                    own dof (also used to locate the same
  !!                                    physical dof in `flux_edge_downwind`)
  !> @param[in]     edge_dist           Distance from this column to the panel
  !!                                    edge on this side (0 if not near it)
  !> @param[in]     direction           +1 for the W/S sides, -1 for the E/N
  !!                                    sides
  !> @param[in]     stencil_max         Maximum size of the stencil branch;
  !!                                    this is the only cap on how far the
  !!                                    redistribution can reach, the actual
  !!                                    depth `n_cells` is derived purely from
  !!                                    the Courant number implied by
  !!                                    `dep_dist`
  !> @param[in]     nlayers             Number of layers
  !> @param[in]     undf_w2h            Num of DoFs for W2h in this partition
  subroutine apply_unification_side( flux_correction, flux,                   &
                                     flux_edge_downwind, dep_dist, own_idx,   &
                                     stencil_branch_f, edge_dist, direction,  &
                                     stencil_max, nlayers,                    &
                                     undf_w2h )

    implicit none

    integer(kind=i_def), intent(in)    :: own_idx
    integer(kind=i_def), intent(in)    :: stencil_max
    integer(kind=i_def), intent(in)    :: stencil_branch_f(stencil_max)
    integer(kind=i_def), intent(in)    :: edge_dist
    integer(kind=i_def), intent(in)    :: direction
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: undf_w2h
    real(kind=r_tran),   intent(inout) :: flux_correction(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux_edge_downwind(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)

    integer(kind=i_def) :: k, l, ref_idx, n_cells
    real(kind=r_tran)   :: seam, courant

    k = edge_dist
    if (k < 1) return
    if (k > 1 .and. k-1 > stencil_max) return

    ! Index of the reference (immediately edge-adjacent) column's own
    ! edge-facing dof: itself if k == 1, otherwise (k-1) steps back towards
    ! the edge along the stencil branch for this side
    if (k == 1) then
      ref_idx = own_idx
    else
      ref_idx = stencil_branch_f(k-1)
    end if

    do l = 0, nlayers-1
      seam = flux(ref_idx+l) - flux_edge_downwind(ref_idx+l)
      if (seam == 0.0_r_tran) cycle

      ! Downwind-ness and the Courant number are both read at this column's
      ! own face (not the edge-adjacent reference face), so no stencil is
      ! needed on dep_dist
      if (REAL(direction, r_tran)*dep_dist(own_idx+l) < 0.0_r_tran) cycle

      courant = ABS(dep_dist(own_idx+l))
      ! n_cells is capped only by dep_dist (via courant); stencil_max+1 is
      ! just the largest depth this kernel invocation can physically reach,
      ! not an independent redistribution-depth parameter
      n_cells = MAX(1_i_def, MIN(stencil_max+1, CEILING(courant + 1.0_r_tran, i_def)))

      if (k < n_cells) then
        flux_correction(own_idx+l) = flux_correction(own_idx+l) +            &
            REAL(direction, r_tran)*seam*REAL(n_cells-k, r_tran)/REAL(n_cells, r_tran)
      end if
    end do

  end subroutine apply_unification_side

end module peregrin_unification_kernel_mod
