!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the advective form update of the inner transported
!!          fields in x and y for PEREGRIN FFSL transport.
!> @details This is the PEREGRIN-specific counterpart of the combination of
!!          `fv_difference_x_kernel_mod.F90` / `fv_difference_y_kernel_mod.F90`
!!          and `swift_inner_update_kernel_mod.F90`: it takes the finite
!!          difference of `flux` at opposite cell faces and divides by the
!!          transported dry mass to give the updated field, for both the x
!!          and y directions in a single kernel.
!!
!!          Away from a panel edge this is identical to the standard (non-
!!          PEREGRIN) calculation. Immediately adjacent to a panel edge
!!          (`panel_edge_dist == 1`), the value used at the edge face is
!!          taken from `flux_edge_downwind` rather than `flux`, for whichever
!!          levels in the column are downwind of the edge (using the sign of
!!          `dep_dist` at that face) - since for those levels `flux` may
!!          hold the upwind (panel-remap) value rather than the downwind
!!          (sphere) one, if any other level in the column was upwind of the
!!          edge (see `peregrin_flux_kernel_mod.F90`).
!!
!> @note This is a first draft (tracker step 3.4.3). The choice of using
!!       `flux_edge_downwind` per-level at the edge face, based on the sign
!!       of `dep_dist` there, has not yet been validated against the
!!       PEREGRIN design and should be reviewed.
!> @note This kernel only works when field is a W3 field at lowest order.

module peregrin_divergence_kernel_mod

  use argument_mod,          only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   CELL_COLUMN, GH_WRITE,     &
                                   GH_READ, GH_SCALAR,        &
                                   GH_INTEGER,                &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,         only: i_def, r_tran, r_def
  use fs_continuity_mod,     only: W3
  use kernel_mod,            only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! PSy layer
  type, public, extends(kernel_type) :: peregrin_divergence_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                        & ! field_x
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                        & ! field_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                        & ! field_n
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! flux
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! flux_edge_downwind
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                        & ! dry_mass_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                        & ! dry_mass_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                        & ! dry_mass_n
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! panel_id
        arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! panel_edge_dist
        arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              & ! dt
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: peregrin_divergence_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: peregrin_divergence_code

contains

  !> @brief Compute the advective form update of the PEREGRIN inner
  !!        transported fields in x and y.
  !> @param[in]     nlayers             Number of layers
  !> @param[in,out] field_x             Updated field after the x sweep
  !> @param[in,out] field_y             Updated field after the y sweep
  !> @param[in]     field_n             Field at the start of the timestep
  !> @param[in]     flux                The mass flux
  !> @param[in]     flux_edge_downwind  Downwind mass flux immediately next
  !!                                    to a panel edge, from
  !!                                    `peregrin_flux_kernel_mod.F90`
  !> @param[in]     dry_mass_x          Transported dry mass after the x sweep
  !> @param[in]     dry_mass_y          Transported dry mass after the y sweep
  !> @param[in]     dry_mass_n          Dry mass at the start of the timestep
  !> @param[in]     dep_dist            Horizontal departure distances
  !> @param[in]     panel_id            Field containing IDs of mesh panels
  !> @param[in]     panel_edge_dist_W   2D field: distance to the panel edge
  !!                                    to the West
  !> @param[in]     panel_edge_dist_E   2D field: distance to the panel edge
  !!                                    to the East
  !> @param[in]     panel_edge_dist_S   2D field: distance to the panel edge
  !!                                    to the South
  !> @param[in]     panel_edge_dist_N   2D field: distance to the panel edge
  !!                                    to the North
  !> @param[in]     dt                  Time step
  !> @param[in]     ndf_w3              Num of DoFs for W3 per cell
  !> @param[in]     undf_w3             Num of DoFs for W3 in this partition
  !> @param[in]     map_w3              Map for W3
  !> @param[in]     ndf_w2h             Num of DoFs for W2h per cell
  !> @param[in]     undf_w2h            Num of DoFs for W2h in this partition
  !> @param[in]     map_w2h             Map for W2h
  !> @param[in]     ndf_pid             Num of DoFs for panel ID field per cell
  !> @param[in]     undf_pid            Num DoFs for this partition for panel_id
  !> @param[in]     map_pid             Map for panel ID field
  subroutine peregrin_divergence_code( nlayers,             &
                                       field_x,             &
                                       field_y,             &
                                       field_n,             &
                                       flux,                &
                                       flux_edge_downwind,  &
                                       dry_mass_x,          &
                                       dry_mass_y,          &
                                       dry_mass_n,          &
                                       dep_dist,            &
                                       panel_id,            &
                                       panel_edge_dist_W,   &
                                       panel_edge_dist_E,   &
                                       panel_edge_dist_S,   &
                                       panel_edge_dist_N,   &
                                       dt,                  &
                                       ndf_w3,              &
                                       undf_w3,             &
                                       map_w3,              &
                                       ndf_w2h,             &
                                       undf_w2h,            &
                                       map_w2h,             &
                                       ndf_pid,             &
                                       undf_pid,            &
                                       map_pid )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_pid
    integer(kind=i_def), intent(in) :: undf_pid
    real(kind=r_tran),   intent(in) :: dt

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: field_x(undf_w3)
    real(kind=r_tran),   intent(inout) :: field_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: field_n(undf_w3)
    real(kind=r_tran),   intent(in)    :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: flux_edge_downwind(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dry_mass_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_n(undf_w3)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)

    ! Internal arguments
    integer(kind=i_def) :: nl, w3_idx, W_idx, E_idx, S_idx, N_idx
    integer(kind=i_def) :: edge_dist_W, edge_dist_E, edge_dist_S, edge_dist_N
    real(kind=r_tran)   :: e_face(0:nlayers-1), w_face(0:nlayers-1)
    real(kind=r_tran)   :: n_face(0:nlayers-1), s_face(0:nlayers-1)

    ! This is based on the lowest order W2h dof map
    !
    !    ---4---
    !    |     |
    !    1     3  horizontal
    !    |     |
    !    ---2---

    w3_idx = map_w3(1)
    W_idx  = map_w2h(1)
    S_idx  = map_w2h(2)
    E_idx  = map_w2h(3)
    N_idx  = map_w2h(4)
    nl = nlayers - 1

    edge_dist_W = panel_edge_dist_W(map_pid(1))
    edge_dist_E = panel_edge_dist_E(map_pid(1))
    edge_dist_S = panel_edge_dist_S(map_pid(1))
    edge_dist_N = panel_edge_dist_N(map_pid(1))

    ! X direction ==============================================================
    e_face = flux(E_idx:E_idx+nl)
    w_face = flux(W_idx:W_idx+nl)

    if (edge_dist_E == 1) then
      ! Levels downwind of the eastern edge (dep_dist <= 0) use the downwind
      ! flux calculated at that face
      e_face = MERGE(                                                          &
          flux_edge_downwind(E_idx:E_idx+nl), e_face,                          &
          dep_dist(E_idx:E_idx+nl) <= 0.0_r_tran                               &
      )
    end if
    if (edge_dist_W == 1) then
      ! Levels downwind of the western edge (dep_dist >= 0) use the downwind
      ! flux calculated at that face
      w_face = MERGE(                                                          &
          flux_edge_downwind(W_idx:W_idx+nl), w_face,                          &
          dep_dist(W_idx:W_idx+nl) >= 0.0_r_tran                               &
      )
    end if

    field_x(w3_idx:w3_idx+nl) = (                                              &
        field_n(w3_idx:w3_idx+nl)*dry_mass_n(w3_idx:w3_idx+nl)                 &
        - dt*(e_face - w_face)                                                 &
    ) / dry_mass_x(w3_idx:w3_idx+nl)

    ! Y direction ==============================================================
    n_face = flux(N_idx:N_idx+nl)
    s_face = flux(S_idx:S_idx+nl)

    if (edge_dist_N == 1) then
      ! Levels downwind of the northern edge (dep_dist >= 0) use the downwind
      ! flux calculated at that face
      n_face = MERGE(                                                          &
          flux_edge_downwind(N_idx:N_idx+nl), n_face,                          &
          dep_dist(N_idx:N_idx+nl) >= 0.0_r_tran                               &
      )
    end if
    if (edge_dist_S == 1) then
      ! Levels downwind of the southern edge (dep_dist <= 0) use the downwind
      ! flux calculated at that face
      s_face = MERGE(                                                          &
          flux_edge_downwind(S_idx:S_idx+nl), s_face,                          &
          dep_dist(S_idx:S_idx+nl) <= 0.0_r_tran                               &
      )
    end if

    field_y(w3_idx:w3_idx+nl) = (                                              &
        field_n(w3_idx:w3_idx+nl)*dry_mass_n(w3_idx:w3_idx+nl)                 &
        - dt*(n_face - s_face)                                                 &
    ) / dry_mass_y(w3_idx:w3_idx+nl)

  end subroutine peregrin_divergence_code

end module peregrin_divergence_kernel_mod
