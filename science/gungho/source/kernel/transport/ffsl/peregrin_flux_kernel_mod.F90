!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the horizontal mass flux for FFSL near cubed-sphere
!!          panel edges, using the PEREGRIN treatment.
!> @details This kernel has no reconstruction logic of its own: it dispatches,
!!          per column and per direction (x using W/E faces, y using S/N
!!          faces), to the existing 1D flux routines:
!!            (a) not near a panel edge: `ffsl_flux_xy_1d`
!!                (`ffsl_flux_xy_kernel_mod.F90`), using `field_for_x`/
!!                `field_for_y`. Writes `flux`.
!!            (b) near a panel edge, upwind column (ANY level in the column
!!                has a departure/integer cell that has not crossed the
!!                edge; a column may contain winds of different directions
!!                at different levels, so this is not necessarily true of
!!                every level): `ffsl_flux_xy_panel_remap_1d`
!!                (`ffsl_flux_xy_panel_remap_kernel_mod.F90`), using the
!!                eave-remapped fields (`remapped_field_x`/`remapped_field_y`)
!!                for the fractional part only. Writes `flux` (unrestricted
!!                face selector).
!!            (c) near a panel edge, downwind column (no level in the column
!!                is upwind of the edge): `ffsl_flux_xy_sphere_1d`
!!                (`ffsl_flux_xy_sphere_kernel_mod.F90`), with the incoming
!!                face selector unchanged. Writes `flux`.
!!            (d) immediately adjacent to a panel edge (regardless of
!!                whether the column as a whole is treated as upwind or
!!                downwind above, since individual levels may still be
!!                downwind of the edge): `ffsl_flux_xy_sphere_1d`, with a
!!                local face selector that only differs from the incoming
!!                one by ensuring the panel-edge face is picked up. Writes
!!                `flux_edge_downwind` (this is later redistributed into
!!                `flux` by `peregrin_unification_kernel_mod.F90`, tracker
!!                step 3.4.4).
!!
!> @note This is a first draft (tracker step 3.4.2). The determination of
!!       "upwind"/"downwind" is done by comparing the sign of the departure
!!       distance at the near face against the direction of the panel edge
!!       (using the same +1/-1 x/y sign convention as the underlying 1D
!!       routines), applied per-level and combined with ANY() across the
!!       column; "immediately adjacent" means the near edge distance is 1.
!!       These rules have not yet been validated against the PEREGRIN design
!!       and should be reviewed.
!> @note This kernel only works when field is a W3 field at lowest order since
!!       it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing the
!!       relevant dofmaps.

module peregrin_flux_kernel_mod

  use argument_mod,                  only : arg_type,                          &
                                            GH_FIELD, GH_REAL,                 &
                                            CELL_COLUMN, GH_WRITE,             &
                                            GH_READ, GH_SCALAR,                &
                                            STENCIL, CROSS2D, GH_INTEGER,      &
                                            ANY_DISCONTINUOUS_SPACE_1,        &
                                            ANY_DISCONTINUOUS_SPACE_2,        &
                                            ANY_DISCONTINUOUS_SPACE_3,        &
                                            ANY_DISCONTINUOUS_SPACE_4
  use constants_mod,                 only : i_def, r_tran, r_def, l_def
  use fs_continuity_mod,             only : W3, W2h
  use kernel_mod,                    only : kernel_type
  use reference_element_mod,         only : W, E, N, S

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: peregrin_flux_kernel_type
    private
    type(arg_type) :: meta_args(22) = (/                                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! flux
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! flux_edge_downwind
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! field_for_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! remapped_field_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! dry_mass_for_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! field_for_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! remapped_field_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),      & ! dry_mass_for_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! frac_dry_flux
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! dep_lowest_k
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! dep_highest_k
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! panel_id
        arg_type(GH_FIELD*2, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! panel_edge_dist_x (W,E)
        arg_type(GH_FIELD*2, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! panel_edge_dist_y (S,N)
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ew
        arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ns
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! order
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! monotone
        arg_type(GH_SCALAR,  GH_REAL,    GH_READ),                             & ! min_val
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             & ! ndep
        arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              & ! dt
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: peregrin_flux_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: peregrin_flux_code

contains

  !> @brief Compute the PEREGRIN horizontal fluxes for FFSL.
  !> @param[in]     nlayers             Number of layers
  !> @param[in,out] flux                The output flux
  !> @param[in,out] flux_edge_downwind  The output flux immediately downwind of
  !!                                    a panel edge, to be redistributed by
  !!                                    `peregrin_unification_kernel_mod.F90`
  !> @param[in]     field_for_x         Field to use in evaluating x-flux
  !> @param[in]     stencil_sizes_x     Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_x       Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_x       Map of DoFs in the stencil for x-field
  !> @param[in]     remapped_field_x    Remapped field in x, used near edges
  !> @param[in]     stencil_sizes_rx    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_rx      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_rx      Map of DoFs in the stencil for x-field
  !> @param[in]     dry_mass_for_x      Volume or dry mass field at W3 points
  !!                                    for use in evaluating x-flux
  !> @param[in]     stencil_sizes_mx    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_mx      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_mx      Map of DoFs in the stencil for x-mass
  !> @param[in]     field_for_y         Field to use in evaluating y-flux
  !> @param[in]     stencil_sizes_y     Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_y       Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_y       Map of DoFs in the stencil for y-field
  !> @param[in]     remapped_field_y    Remapped field in y, used near edges
  !> @param[in]     stencil_sizes_ry    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_ry      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_ry      Map of DoFs in the stencil for y-field
  !> @param[in]     dry_mass_for_y      Volume or dry mass field at W3 points
  !!                                    for use in evaluating y-flux
  !> @param[in]     stencil_sizes_my    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_my      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_my      Map of DoFs in the stencil for y-mass
  !> @param[in]     dep_dist            Horizontal departure distances
  !> @param[in]     frac_dry_flux       Fractional part of the dry flux or wind
  !> @param[in]     dep_lowest_k        2D integer multidata W2H field, storing
  !!                                    the lowest model level to use in each
  !!                                    integer flux sum, for each column
  !> @param[in]     dep_highest_k       2D integer multidata W2H field, storing
  !!                                    the highest model level to use in each
  !!                                    integer flux sum, for each column
  !> @param[in]     panel_id            Field containing IDs of mesh panels
  !> @param[in]     panel_edge_dist_W   2D field: distance to the panel edge
  !!                                    to the West
  !> @param[in]     panel_edge_dist_E   2D field: distance to the panel edge
  !!                                    to the East
  !> @param[in]     stencil_sizes_cx    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_cx      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_cx      Map of DoFs in the stencil for x-case
  !> @param[in]     panel_edge_dist_S   2D field: distance to the panel edge
  !!                                    to the South
  !> @param[in]     panel_edge_dist_N   2D field: distance to the panel edge
  !!                                    to the North
  !> @param[in]     stencil_sizes_cy    Sizes of branches of the cross stencil
  !> @param[in]     stencil_max_cy      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_cy      Map of DoFs in the stencil for y-case
  !> @param[in]     face_selector_ew    2D field indicating which W/E faces to
  !!                                    loop over for this column
  !> @param[in]     face_selector_ns    2D field indicating which N/S faces to
  !!                                    loop over for this column
  !> @param[in]     order               Order of reconstruction
  !> @param[in]     monotone            Horizontal monotone option for FFSL
  !> @param[in]     min_val             Minimum value to enforce when using
  !!                                    quasi-monotone limiter
  !> @param[in]     ndep                Number of multidata points for departure
  !!                                    index fields
  !> @param[in]     dt                  Time step
  !> @param[in]     ndf_w2h             Num of DoFs for W2h per cell
  !> @param[in]     undf_w2h            Num of DoFs for W2h in this partition
  !> @param[in]     map_w2h             Map for W2h
  !> @param[in]     ndf_w3              Num of DoFs for W3 per cell
  !> @param[in]     undf_w3             Num of DoFs for W3 in this partition
  !> @param[in]     map_w3              Map for W3
  !> @param[in]     ndf_depk            Num of DoFs for dep idx fields per cell
  !> @param[in]     undf_depk           Num of DoFs for this partition for dep
  !!                                    idx fields
  !> @param[in]     map_depk            Map for departure index fields
  !> @param[in]     ndf_pid             Num of DoFs for panel ID field per cell
  !> @param[in]     undf_pid            Num DoFs for this partition for panel_id
  !> @param[in]     map_pid             Map for panel ID field
  !> @param[in]     ndf_w3_2d           Num of DoFs for 2D W3 per cell
  !> @param[in]     undf_w3_2d          Num of DoFs for this partition for 2D W3
  !> @param[in]     map_w3_2d           Map for 2D W3
  subroutine peregrin_flux_code( nlayers,             &
                                 flux,                &
                                 flux_edge_downwind,  &
                                 field_for_x,         &
                                 stencil_sizes_x,     &
                                 stencil_max_x,       &
                                 stencil_map_x,       &
                                 remapped_field_x,    &
                                 stencil_sizes_rx,    &
                                 stencil_max_rx,      &
                                 stencil_map_rx,      &
                                 dry_mass_for_x,      &
                                 stencil_sizes_mx,    &
                                 stencil_max_mx,      &
                                 stencil_map_mx,      &
                                 field_for_y,         &
                                 stencil_sizes_y,     &
                                 stencil_max_y,       &
                                 stencil_map_y,       &
                                 remapped_field_y,    &
                                 stencil_sizes_ry,    &
                                 stencil_max_ry,      &
                                 stencil_map_ry,      &
                                 dry_mass_for_y,      &
                                 stencil_sizes_my,    &
                                 stencil_max_my,      &
                                 stencil_map_my,      &
                                 dep_dist,            &
                                 frac_dry_flux,       &
                                 dep_lowest_k,        &
                                 dep_highest_k,       &
                                 panel_id,            &
                                 panel_edge_dist_W,   &
                                 panel_edge_dist_E,   &
                                 stencil_sizes_cx,    &
                                 stencil_max_cx,      &
                                 stencil_map_cx,      &
                                 panel_edge_dist_S,   &
                                 panel_edge_dist_N,   &
                                 stencil_sizes_cy,    &
                                 stencil_max_cy,      &
                                 stencil_map_cy,      &
                                 face_selector_ew,    &
                                 face_selector_ns,    &
                                 order,               &
                                 monotone,            &
                                 min_val,             &
                                 ndep,                &
                                 dt,                  &
                                 ndf_w2h,             &
                                 undf_w2h,            &
                                 map_w2h,             &
                                 ndf_w3,              &
                                 undf_w3,             &
                                 map_w3,              &
                                 ndf_depk,            &
                                 undf_depk,           &
                                 map_depk,            &
                                 ndf_pid,             &
                                 undf_pid,            &
                                 map_pid,             &
                                 ndf_w3_2d,           &
                                 undf_w3_2d,          &
                                 map_w3_2d )

    use ffsl_flux_xy_kernel_mod,             only: ffsl_flux_xy_1d
    use ffsl_flux_xy_sphere_kernel_mod,      only: ffsl_flux_xy_sphere_1d
    use ffsl_flux_xy_panel_remap_kernel_mod, only: ffsl_flux_xy_panel_remap_1d
    use panel_edge_support_mod,              only: crosses_panel_edge

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_pid
    integer(kind=i_def), intent(in) :: undf_pid
    integer(kind=i_def), intent(in) :: ndf_depk
    integer(kind=i_def), intent(in) :: undf_depk
    integer(kind=i_def), intent(in) :: ndep
    integer(kind=i_def), intent(in) :: stencil_max_x
    integer(kind=i_def), intent(in) :: stencil_max_rx
    integer(kind=i_def), intent(in) :: stencil_max_mx
    integer(kind=i_def), intent(in) :: stencil_max_y
    integer(kind=i_def), intent(in) :: stencil_max_ry
    integer(kind=i_def), intent(in) :: stencil_max_my
    integer(kind=i_def), intent(in) :: stencil_max_cx
    integer(kind=i_def), intent(in) :: stencil_max_cy
    integer(kind=i_def), intent(in) :: stencil_sizes_x(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_rx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_mx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_y(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_ry(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_my(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_cx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_cy(4)
    integer(kind=i_def), intent(in) :: order
    integer(kind=i_def), intent(in) :: monotone
    real(kind=r_tran),   intent(in) :: min_val
    real(kind=r_tran),   intent(in) :: dt

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def), intent(in) :: map_depk(ndf_depk)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_w3,stencil_max_x,4)
    integer(kind=i_def), intent(in) :: stencil_map_rx(ndf_w3,stencil_max_rx,4)
    integer(kind=i_def), intent(in) :: stencil_map_mx(ndf_w3,stencil_max_mx,4)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_w3,stencil_max_y,4)
    integer(kind=i_def), intent(in) :: stencil_map_ry(ndf_w3,stencil_max_ry,4)
    integer(kind=i_def), intent(in) :: stencil_map_my(ndf_w3,stencil_max_my,4)
    integer(kind=i_def), intent(in) :: stencil_map_cx(ndf_pid,stencil_max_cx,4)
    integer(kind=i_def), intent(in) :: stencil_map_cy(ndf_pid,stencil_max_cy,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux(undf_w2h)
    real(kind=r_tran),   intent(inout) :: flux_edge_downwind(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: field_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: remapped_field_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: remapped_field_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_tran),   intent(in)    :: frac_dry_flux(undf_w2h)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: dep_lowest_k(undf_depk)
    integer(kind=i_def), intent(in)    :: dep_highest_k(undf_depk)

    ! Internal arguments
    integer(kind=i_def) :: i, ipanel
    integer(kind=i_def) :: edge_dist_E, edge_dist_W, edge_dist_S, edge_dist_N
    integer(kind=i_def) :: near_dist
    logical(kind=l_def) :: near_edge, is_upwind, is_immediate
    integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
    integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
    integer(kind=i_def) :: stencil_map_x_1d(-stencil_max_x:stencil_max_x)
    integer(kind=i_def) :: stencil_map_y_1d(-stencil_max_y:stencil_max_y)
    integer(kind=i_def) :: stencil_map_cx_1d(-stencil_max_x:stencil_max_x)
    integer(kind=i_def) :: stencil_map_cy_1d(-stencil_max_y:stencil_max_y)
    integer(kind=i_def) :: face_selector_edge(undf_w3_2d)
    integer(kind=i_def) :: edge_face

    ! Form X and Y 1D stencils
    stencil_extent_xl = stencil_sizes_x(1) - 1
    stencil_extent_xr = stencil_sizes_x(3) - 1
    stencil_extent_yl = stencil_sizes_y(2) - 1
    stencil_extent_yr = stencil_sizes_y(4) - 1

    do i = -stencil_extent_xl, 0
      stencil_map_x_1d(i) = stencil_map_x(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_x_1d(i) = stencil_map_x(1, i+1, 3)
    end do

    do i = -stencil_extent_yl, 0
      stencil_map_y_1d(i) = stencil_map_y(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_y_1d(i) = stencil_map_y(1, i+1, 4)
    end do

    do i = -stencil_extent_xl, 0
      stencil_map_cx_1d(i) = stencil_map_cx(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_cx_1d(i) = stencil_map_cx(1, i+1, 3)
    end do

    do i = -stencil_extent_yl, 0
      stencil_map_cy_1d(i) = stencil_map_cy(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_cy_1d(i) = stencil_map_cy(1, i+1, 4)
    end do

    ipanel = INT(panel_id(map_pid(1)), i_def)
    edge_dist_E = panel_edge_dist_E(map_pid(1))
    edge_dist_W = panel_edge_dist_W(map_pid(1))
    edge_dist_S = panel_edge_dist_S(map_pid(1))
    edge_dist_N = panel_edge_dist_N(map_pid(1))

    ! X direction ==============================================================
    near_edge = crosses_panel_edge(                                            &
        edge_dist_W, edge_dist_E, MAX(stencil_extent_xl, stencil_extent_xr),   &
        order, face_selector_ew(map_w3_2d(1)), .false., ipanel, 1, W, E,       &
        dep_highest_k, ndep, ndf_depk, undf_depk, map_depk                     &
    )

    if (.not. near_edge) then
      ! (a) Not near a panel edge: default to standard FFSL
      call ffsl_flux_xy_1d( nlayers,             &
                            .true.,              &
                            flux,                &
                            field_for_x,         &
                            stencil_extent_xl,   &
                            stencil_extent_xr,   &
                            stencil_max_x,       &
                            stencil_map_x_1d,    &
                            dry_mass_for_x,      &
                            dep_dist,            &
                            frac_dry_flux,       &
                            dep_lowest_k,        &
                            dep_highest_k,       &
                            face_selector_ew,    &
                            order,               &
                            monotone,            &
                            min_val,             &
                            ndep,                &
                            dt,                  &
                            ndf_w2h,             &
                            undf_w2h,            &
                            map_w2h,             &
                            ndf_w3,              &
                            undf_w3,             &
                            map_w3,              &
                            ndf_depk,            &
                            undf_depk,           &
                            map_depk,            &
                            ndf_w3_2d,           &
                            undf_w3_2d,          &
                            map_w3_2d )
    else
      ! Work out which side is nearest, and whether this column is upwind or
      ! downwind of that edge (see module-level note on this assumption). A
      ! column is only treated as downwind if NONE of its levels are upwind
      ! of the edge; if ANY level is upwind, the whole column uses the
      ! upwind (panel-remap) treatment for `flux`, since a column may
      ! contain winds of different directions at different levels.
      if (edge_dist_E > 0 .and. (edge_dist_W <= 0 .or. edge_dist_E < edge_dist_W)) then
        near_dist = edge_dist_E
        edge_face = E
        is_upwind = ANY(dep_dist(map_w2h(E):map_w2h(E)+nlayers-1) > 0.0_r_tran)
      else
        near_dist = edge_dist_W
        edge_face = W
        is_upwind = ANY(dep_dist(map_w2h(W):map_w2h(W)+nlayers-1) < 0.0_r_tran)
      end if
      is_immediate = (near_dist == 1)

      if (is_upwind) then
        ! (b) Upwind side: use the eave-remapped fields in the fractional
        ! reconstruction, ensuring the panel edge flux is calculated
        call ffsl_flux_xy_panel_remap_1d( nlayers,             &
                                          .true.,              &
                                          flux,                &
                                          field_for_x,         &
                                          field_for_y,         &
                                          remapped_field_x,    &
                                          remapped_field_y,    &
                                          stencil_extent_xl,   &
                                          stencil_extent_xr,   &
                                          stencil_max_x,       &
                                          stencil_map_x_1d,    &
                                          dry_mass_for_x,      &
                                          dry_mass_for_y,      &
                                          dep_dist,            &
                                          frac_dry_flux,       &
                                          dep_lowest_k,        &
                                          dep_highest_k,       &
                                          ipanel,              &
                                          panel_edge_dist_W,   &
                                          panel_edge_dist_E,   &
                                          panel_edge_dist_S,   &
                                          panel_edge_dist_N,   &
                                          stencil_map_cx_1d,   &
                                          face_selector_ew,    &
                                          order,               &
                                          monotone,            &
                                          min_val,             &
                                          ndep,                &
                                          dt,                  &
                                          ndf_w2h,             &
                                          undf_w2h,            &
                                          map_w2h,             &
                                          ndf_w3,              &
                                          undf_w3,             &
                                          map_w3,              &
                                          ndf_depk,            &
                                          undf_depk,           &
                                          map_depk,            &
                                          ndf_pid,             &
                                          undf_pid,            &
                                          map_pid,             &
                                          ndf_w3_2d,           &
                                          undf_w3_2d,          &
                                          map_w3_2d )
      else
        ! (c) Downwind side: use the eave-remapped fields directly (via
        ! ffsl_flux_xy_sphere_1d), with the incoming face selector unchanged
        call ffsl_flux_xy_sphere_1d( nlayers,             &
                                     .true.,              &
                                     flux,                &
                                     field_for_x,         &
                                     field_for_y,         &
                                     stencil_extent_xl,   &
                                     stencil_extent_xr,   &
                                     stencil_max_x,       &
                                     stencil_map_x_1d,    &
                                     dry_mass_for_x,      &
                                     dry_mass_for_y,      &
                                     dep_dist,            &
                                     frac_dry_flux,       &
                                     dep_lowest_k,        &
                                     dep_highest_k,       &
                                     ipanel,              &
                                     edge_dist_W,         &
                                     edge_dist_E,         &
                                     face_selector_ew,    &
                                     order,               &
                                     monotone,            &
                                     min_val,             &
                                     ndep,                &
                                     dt,                  &
                                     ndf_w2h,             &
                                     undf_w2h,            &
                                     map_w2h,             &
                                     ndf_w3,              &
                                     undf_w3,             &
                                     map_w3,              &
                                     ndf_depk,            &
                                     undf_depk,           &
                                     map_depk,            &
                                     ndf_w3_2d,           &
                                     undf_w3_2d,          &
                                     map_w3_2d )
      end if

      ! (d) Immediately adjacent to a panel edge: also calculate the flux
      ! across the edge face itself into `flux_edge_downwind`, for later
      ! redistribution by peregrin_unification_kernel_mod.F90 (tracker step
      ! 3.4.4). This is done regardless of whether the column as a whole is
      ! upwind or downwind above, since individual levels within it may be
      ! downwind of the edge even when the column as a whole is treated as
      ! upwind. The local face selector only needs to differ from the
      ! incoming one in this edge-adjacent cell, to ensure the edge face is
      ! picked up here even if it would not otherwise be computed for this
      ! cell.
      if (is_immediate) then
        face_selector_edge = face_selector_ew
        face_selector_edge(map_w3_2d(1)) =                                     &
            restrict_face_selector(face_selector_ew(map_w3_2d(1)), edge_face)

        call ffsl_flux_xy_sphere_1d( nlayers,             &
                                     .true.,              &
                                     flux_edge_downwind,  &
                                     field_for_x,         &
                                     field_for_y,         &
                                     stencil_extent_xl,   &
                                     stencil_extent_xr,   &
                                     stencil_max_x,       &
                                     stencil_map_x_1d,    &
                                     dry_mass_for_x,      &
                                     dry_mass_for_y,      &
                                     dep_dist,            &
                                     frac_dry_flux,       &
                                     dep_lowest_k,        &
                                     dep_highest_k,       &
                                     ipanel,              &
                                     edge_dist_W,         &
                                     edge_dist_E,         &
                                     face_selector_edge,  &
                                     order,               &
                                     monotone,            &
                                     min_val,             &
                                     ndep,                &
                                     dt,                  &
                                     ndf_w2h,             &
                                     undf_w2h,            &
                                     map_w2h,             &
                                     ndf_w3,              &
                                     undf_w3,             &
                                     map_w3,              &
                                     ndf_depk,            &
                                     undf_depk,           &
                                     map_depk,            &
                                     ndf_w3_2d,           &
                                     undf_w3_2d,          &
                                     map_w3_2d )
      end if
    end if

    ! Y direction ==============================================================
    near_edge = crosses_panel_edge(                                            &
        edge_dist_S, edge_dist_N, MAX(stencil_extent_yl, stencil_extent_yr),   &
        order, face_selector_ns(map_w3_2d(1)), .false., ipanel, 2, S, N,       &
        dep_highest_k, ndep, ndf_depk, undf_depk, map_depk                     &
    )

    if (.not. near_edge) then
      ! (a) Not near a panel edge: default to standard FFSL
      call ffsl_flux_xy_1d( nlayers,             &
                            .false.,             &
                            flux,                &
                            field_for_y,         &
                            stencil_extent_yl,   &
                            stencil_extent_yr,   &
                            stencil_max_y,       &
                            stencil_map_y_1d,    &
                            dry_mass_for_y,      &
                            dep_dist,            &
                            frac_dry_flux,       &
                            dep_lowest_k,        &
                            dep_highest_k,       &
                            face_selector_ns,    &
                            order,               &
                            monotone,            &
                            min_val,             &
                            ndep,                &
                            dt,                  &
                            ndf_w2h,             &
                            undf_w2h,            &
                            map_w2h,             &
                            ndf_w3,              &
                            undf_w3,             &
                            map_w3,              &
                            ndf_depk,            &
                            undf_depk,           &
                            map_depk,            &
                            ndf_w3_2d,           &
                            undf_w3_2d,          &
                            map_w3_2d )
    else
      ! Work out which side is nearest, and whether this column is upwind or
      ! downwind of that edge. The y-direction sign convention is reversed
      ! relative to x (see module-level note on this assumption). A column
      ! is only treated as downwind if NONE of its levels are upwind of the
      ! edge; if ANY level is upwind, the whole column uses the upwind
      ! (panel-remap) treatment for `flux`, since a column may contain winds
      ! of different directions at different levels.
      if (edge_dist_N > 0 .and. (edge_dist_S <= 0 .or. edge_dist_N < edge_dist_S)) then
        near_dist = edge_dist_N
        edge_face = N
        is_upwind = ANY(dep_dist(map_w2h(N):map_w2h(N)+nlayers-1) < 0.0_r_tran)
      else
        near_dist = edge_dist_S
        edge_face = S
        is_upwind = ANY(dep_dist(map_w2h(S):map_w2h(S)+nlayers-1) > 0.0_r_tran)
      end if
      is_immediate = (near_dist == 1)

      if (is_upwind) then
        ! (b) Upwind side: use the eave-remapped fields in the fractional
        ! reconstruction, ensuring the panel edge flux is calculated
        call ffsl_flux_xy_panel_remap_1d( nlayers,             &
                                          .false.,             &
                                          flux,                &
                                          field_for_y,         &
                                          field_for_x,         &
                                          remapped_field_y,    &
                                          remapped_field_x,    &
                                          stencil_extent_yl,   &
                                          stencil_extent_yr,   &
                                          stencil_max_y,       &
                                          stencil_map_y_1d,    &
                                          dry_mass_for_y,      &
                                          dry_mass_for_x,      &
                                          dep_dist,            &
                                          frac_dry_flux,       &
                                          dep_lowest_k,        &
                                          dep_highest_k,       &
                                          ipanel,              &
                                          panel_edge_dist_S,   &
                                          panel_edge_dist_N,   &
                                          panel_edge_dist_W,   &
                                          panel_edge_dist_E,   &
                                          stencil_map_cy_1d,   &
                                          face_selector_ns,    &
                                          order,               &
                                          monotone,            &
                                          min_val,             &
                                          ndep,                &
                                          dt,                  &
                                          ndf_w2h,             &
                                          undf_w2h,            &
                                          map_w2h,             &
                                          ndf_w3,              &
                                          undf_w3,             &
                                          map_w3,              &
                                          ndf_depk,            &
                                          undf_depk,           &
                                          map_depk,            &
                                          ndf_pid,             &
                                          undf_pid,            &
                                          map_pid,             &
                                          ndf_w3_2d,           &
                                          undf_w3_2d,          &
                                          map_w3_2d )
      else
        ! (c) Downwind side: as in the x-direction above, with the incoming
        ! face selector unchanged
        call ffsl_flux_xy_sphere_1d( nlayers,             &
                                     .false.,             &
                                     flux,                &
                                     field_for_y,         &
                                     field_for_x,         &
                                     stencil_extent_yl,   &
                                     stencil_extent_yr,   &
                                     stencil_max_y,       &
                                     stencil_map_y_1d,    &
                                     dry_mass_for_y,      &
                                     dry_mass_for_x,      &
                                     dep_dist,            &
                                     frac_dry_flux,       &
                                     dep_lowest_k,        &
                                     dep_highest_k,       &
                                     ipanel,              &
                                     edge_dist_S,         &
                                     edge_dist_N,         &
                                     face_selector_ns,    &
                                     order,               &
                                     monotone,            &
                                     min_val,             &
                                     ndep,                &
                                     dt,                  &
                                     ndf_w2h,             &
                                     undf_w2h,            &
                                     map_w2h,             &
                                     ndf_w3,              &
                                     undf_w3,             &
                                     map_w3,              &
                                     ndf_depk,            &
                                     undf_depk,           &
                                     map_depk,            &
                                     ndf_w3_2d,           &
                                     undf_w3_2d,          &
                                     map_w3_2d )
      end if

      ! (d) Immediately adjacent to a panel edge: also calculate the flux
      ! across the edge face itself into `flux_edge_downwind`, as in the
      ! x-direction above
      if (is_immediate) then
        face_selector_edge = face_selector_ns
        face_selector_edge(map_w3_2d(1)) =                                     &
            restrict_face_selector(face_selector_ns(map_w3_2d(1)), edge_face)

        call ffsl_flux_xy_sphere_1d( nlayers,             &
                                     .false.,             &
                                     flux_edge_downwind,  &
                                     field_for_y,         &
                                     field_for_x,         &
                                     stencil_extent_yl,   &
                                     stencil_extent_yr,   &
                                     stencil_max_y,       &
                                     stencil_map_y_1d,    &
                                     dry_mass_for_y,      &
                                     dry_mass_for_x,      &
                                     dep_dist,            &
                                     frac_dry_flux,       &
                                     dep_lowest_k,        &
                                     dep_highest_k,       &
                                     ipanel,              &
                                     edge_dist_S,         &
                                     edge_dist_N,         &
                                     face_selector_edge,  &
                                     order,               &
                                     monotone,            &
                                     min_val,             &
                                     ndep,                &
                                     dt,                  &
                                     ndf_w2h,             &
                                     undf_w2h,            &
                                     map_w2h,             &
                                     ndf_w3,              &
                                     undf_w3,             &
                                     map_w3,              &
                                     ndf_depk,            &
                                     undf_depk,           &
                                     map_depk,            &
                                     ndf_w3_2d,           &
                                     undf_w3_2d,          &
                                     map_w3_2d )
      end if
    end if

  end subroutine peregrin_flux_code

  !> @brief Restrict a face selector value to a single face.
  !> @param[in] fs_val    The original face selector value (1 = left face
  !!                      only, 2 = both faces, -1 = right face only)
  !> @param[in] keep_face Which face to keep: W or S for the "left" face,
  !!                      E or N for the "right" face
  !> @return    The restricted face selector value (1, -1 or 0 if the
  !!            requested face was not present in the original selector)
  function restrict_face_selector(fs_val, keep_face) result(new_fs_val)

    implicit none

    integer(kind=i_def), intent(in) :: fs_val
    integer(kind=i_def), intent(in) :: keep_face
    integer(kind=i_def) :: new_fs_val

    if (keep_face == W .or. keep_face == S) then
      if (fs_val == 1 .or. fs_val == 2) then
        new_fs_val = 1
      else
        new_fs_val = 0
      end if
    else
      if (fs_val == -1 .or. fs_val == 2) then
        new_fs_val = -1
      else
        new_fs_val = 0
      end if
    end if

  end function restrict_face_selector

end module peregrin_flux_kernel_mod
