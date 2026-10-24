!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief   Calculates the horizontal mass flux for FFSL near cubed-sphere
!!          panel edges, using the PEREGRIN treatment.
!> @details Per column and direction (x using W/E faces, y using S/N faces),
!!          uses `ffsl_flux_xy_1d` away from panel edges and
!!          `ffsl_flux_xy_panel_remap_1d` for normal fluxes near an edge. For
!!          a face on the panel boundary, `peregrin_flux_1d` computes
!!          `flux_edge_downwind`, selecting remapped fields for stencil cells
!!          across the edge. That flux is later redistributed into `flux` by
!!          `peregrin_unification_kernel_mod.F90`.
!!
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
    integer(kind=i_def) :: face_idx, face
    integer(kind=i_def) :: face_dofs(2)
    logical(kind=l_def) :: near_edge
    logical(kind=l_def) :: is_panel_edge
    integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
    integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
    integer(kind=i_def) :: stencil_map_x_1d(-stencil_max_x:stencil_max_x)
    integer(kind=i_def) :: stencil_map_y_1d(-stencil_max_y:stencil_max_y)
    integer(kind=i_def) :: stencil_map_cx_1d(-stencil_max_x:stencil_max_x)
    integer(kind=i_def) :: stencil_map_cy_1d(-stencil_max_y:stencil_max_y)
    integer(kind=i_def) :: face_selector_edge(undf_w3_2d)

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
    face_dofs = (/ W, E /)
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
      do face_idx = 1, ABS(face_selector_ew(map_w3_2d(1)))
        face = face_dofs(face_idx - MIN(0, face_selector_ew(map_w3_2d(1))))
        is_panel_edge = (face == W .and. edge_dist_W == 1) .or.                &
                        (face == E .and. edge_dist_E == 1)
        if (is_panel_edge) then
          call peregrin_flux_1d(                                               &
              nlayers, .true., face, flux_edge_downwind, field_for_x,          &
              field_for_y, remapped_field_x, remapped_field_y,                 &
              stencil_extent_xl, stencil_extent_xr, stencil_max_x,             &
              stencil_map_x_1d, stencil_map_cx_1d, dry_mass_for_x,             &
              dry_mass_for_y, dep_dist, frac_dry_flux, dep_lowest_k,           &
              dep_highest_k, ipanel, panel_edge_dist_W, panel_edge_dist_E,     &
              panel_edge_dist_S, panel_edge_dist_N, order, monotone, min_val,  &
              ndep, dt, ndf_w2h, undf_w2h, map_w2h, ndf_w3, undf_w3, map_w3,   &
              ndf_depk, undf_depk, map_depk, ndf_pid, undf_pid, map_pid,       &
              ndf_w3_2d, undf_w3_2d, map_w3_2d)

          face_selector_edge = 0
          if (face == W) then
            face_selector_edge(map_w3_2d(1)) = 1
          else
            face_selector_edge(map_w3_2d(1)) = -1
          end if
          call ffsl_flux_xy_panel_remap_1d(                                    &
              nlayers, .true., flux, field_for_x, field_for_y,                 &
              remapped_field_x, remapped_field_y, stencil_extent_xl,           &
              stencil_extent_xr, stencil_max_x, stencil_map_x_1d,              &
              dry_mass_for_x, dry_mass_for_y, dep_dist, frac_dry_flux,         &
              dep_lowest_k, dep_highest_k, ipanel, panel_edge_dist_W,          &
              panel_edge_dist_E, panel_edge_dist_S, panel_edge_dist_N,         &
              stencil_map_cx_1d, face_selector_edge, order, monotone,          &
              min_val, ndep, dt, ndf_w2h, undf_w2h, map_w2h, ndf_w3, undf_w3,  &
              map_w3, ndf_depk, undf_depk, map_depk, ndf_pid, undf_pid,        &
              map_pid, ndf_w3_2d, undf_w3_2d, map_w3_2d)
        else
          call peregrin_flux_1d(                                               &
              nlayers, .true., face, flux, field_for_x, field_for_y,           &
              remapped_field_x, remapped_field_y, stencil_extent_xl,           &
              stencil_extent_xr, stencil_max_x, stencil_map_x_1d,              &
              stencil_map_cx_1d, dry_mass_for_x, dry_mass_for_y, dep_dist,     &
              frac_dry_flux, dep_lowest_k, dep_highest_k, ipanel,              &
              panel_edge_dist_W, panel_edge_dist_E, panel_edge_dist_S,         &
              panel_edge_dist_N, order, monotone, min_val, ndep, dt,           &
              ndf_w2h, undf_w2h, map_w2h, ndf_w3, undf_w3, map_w3,             &
              ndf_depk, undf_depk, map_depk, ndf_pid, undf_pid, map_pid,       &
              ndf_w3_2d, undf_w3_2d, map_w3_2d)
        end if
      end do
    end if

    ! Y direction ==============================================================
    face_dofs = (/ S, N /)
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
      do face_idx = 1, ABS(face_selector_ns(map_w3_2d(1)))
        face = face_dofs(face_idx - MIN(0, face_selector_ns(map_w3_2d(1))))
        is_panel_edge = (face == S .and. edge_dist_S == 1) .or.                &
                        (face == N .and. edge_dist_N == 1)
        if (is_panel_edge) then
          call peregrin_flux_1d(                                               &
              nlayers, .false., face, flux_edge_downwind, field_for_y,         &
              field_for_x, remapped_field_y, remapped_field_x,                 &
              stencil_extent_yl, stencil_extent_yr, stencil_max_y,             &
              stencil_map_y_1d, stencil_map_cy_1d, dry_mass_for_y,             &
              dry_mass_for_x, dep_dist, frac_dry_flux, dep_lowest_k,           &
              dep_highest_k, ipanel, panel_edge_dist_S, panel_edge_dist_N,     &
              panel_edge_dist_E, panel_edge_dist_W, order, monotone, min_val,  &
              ndep, dt, ndf_w2h, undf_w2h, map_w2h, ndf_w3, undf_w3, map_w3,   &
              ndf_depk, undf_depk, map_depk, ndf_pid, undf_pid, map_pid,       &
              ndf_w3_2d, undf_w3_2d, map_w3_2d)

          face_selector_edge = 0
          if (face == S) then
            face_selector_edge(map_w3_2d(1)) = 1
          else
            face_selector_edge(map_w3_2d(1)) = -1
          end if
          call ffsl_flux_xy_panel_remap_1d(                                    &
              nlayers, .false., flux, field_for_y, field_for_x,                &
              remapped_field_y, remapped_field_x, stencil_extent_yl,           &
              stencil_extent_yr, stencil_max_y, stencil_map_y_1d,              &
              dry_mass_for_y, dry_mass_for_x, dep_dist, frac_dry_flux,         &
              dep_lowest_k, dep_highest_k, ipanel, panel_edge_dist_S,          &
              panel_edge_dist_N, panel_edge_dist_E, panel_edge_dist_W,         &
              stencil_map_cy_1d, face_selector_edge, order, monotone,          &
              min_val, ndep, dt, ndf_w2h, undf_w2h, map_w2h, ndf_w3, undf_w3,  &
              map_w3, ndf_depk, undf_depk, map_depk, ndf_pid, undf_pid,        &
              map_pid, ndf_w3_2d, undf_w3_2d, map_w3_2d)
        else
          call peregrin_flux_1d(                                               &
              nlayers, .false., face, flux, field_for_y, field_for_x,          &
              remapped_field_y, remapped_field_x, stencil_extent_yl,           &
              stencil_extent_yr, stencil_max_y, stencil_map_y_1d,              &
              stencil_map_cy_1d, dry_mass_for_y, dry_mass_for_x, dep_dist,     &
              frac_dry_flux, dep_lowest_k, dep_highest_k, ipanel,              &
              panel_edge_dist_S, panel_edge_dist_N, panel_edge_dist_E,         &
              panel_edge_dist_W, order, monotone, min_val, ndep, dt,           &
              ndf_w2h, undf_w2h, map_w2h, ndf_w3, undf_w3, map_w3,             &
              ndf_depk, undf_depk, map_depk, ndf_pid, undf_pid, map_pid,       &
              ndf_w3_2d, undf_w3_2d, map_w3_2d)
        end if
      end do
    end if

  end subroutine peregrin_flux_code

  !> @brief Compute one flux face, including remapped integer contributions.
  !> @param[in]     face                 Face for this flux calculation
  subroutine peregrin_flux_1d(nlayers, x_direction, face, flux,                &
                              field, field_swapped, remapped_field,            &
                              remapped_field_swapped, stencil_extent_l,        &
                              stencil_extent_r, stencil_max, stencil_map,      &
                              stencil_map_pid, dry_mass, dry_mass_swapped,     &
                              dep_dist, frac_dry_flux, dep_lowest_k,           &
                              dep_highest_k, ipanel, panel_edge_dist_l,        &
                              panel_edge_dist_r, panel_edge_dist_swapped_l,    &
                              panel_edge_dist_swapped_r, order, monotone,      &
                              min_val, ndep, dt, ndf_w2h, undf_w2h, map_w2h,   &
                              ndf_w3, undf_w3, map_w3, ndf_depk, undf_depk,    &
                              map_depk, ndf_pid, undf_pid, map_pid,            &
                              ndf_w3_2d, undf_w3_2d, map_w3_2d)

    use panel_edge_support_mod,         only: rotated_panel_neighbour, FAR_AWAY
    use subgrid_horizontal_support_mod, only: fourth_order_horizontal_edge,    &
                                              nirvana_horizontal_edge
    use subgrid_common_support_mod,     only: monotonic_edge,                  &
                                              subgrid_quadratic_recon

    implicit none

    integer(kind=i_def),       intent(in)    :: nlayers, face
    integer(kind=i_def),       intent(in)    :: ndf_w2h, undf_w2h, ndf_w3, undf_w3
    integer(kind=i_def),       intent(in)    :: ndf_depk, undf_depk, ndf_pid, undf_pid
    integer(kind=i_def),       intent(in)    :: ndf_w3_2d, undf_w3_2d, ndep
    integer(kind=i_def),       intent(in)    :: stencil_extent_l, stencil_extent_r
    integer(kind=i_def),       intent(in)    :: stencil_max, ipanel, order, monotone
    integer(kind=i_def),       intent(in)    :: map_w2h(ndf_w2h), map_w3(ndf_w3)
    integer(kind=i_def),       intent(in)    :: map_depk(ndf_depk), map_pid(ndf_pid)
    integer(kind=i_def),       intent(in)    :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def),       intent(in)    :: stencil_map(-stencil_max:stencil_max)
    integer(kind=i_def),       intent(in)    :: stencil_map_pid(-stencil_max:stencil_max)
    integer(kind=i_def),       intent(in)    :: panel_edge_dist_l(undf_pid)
    integer(kind=i_def),       intent(in)    :: panel_edge_dist_r(undf_pid)
    integer(kind=i_def),       intent(in)    :: panel_edge_dist_swapped_l(undf_pid)
    integer(kind=i_def),       intent(in)    :: panel_edge_dist_swapped_r(undf_pid)
    integer(kind=i_def),       intent(in)    :: dep_lowest_k(undf_depk)
    integer(kind=i_def),       intent(in)    :: dep_highest_k(undf_depk)
    real(kind=r_tran),         intent(inout) :: flux(undf_w2h)
    real(kind=r_tran), target, intent(in)    :: field(undf_w3)
    real(kind=r_tran), target, intent(in)    :: field_swapped(undf_w3)
    real(kind=r_tran), target, intent(in)    :: remapped_field(undf_w3)
    real(kind=r_tran), target, intent(in)    :: remapped_field_swapped(undf_w3)
    real(kind=r_tran), target, intent(in)    :: dry_mass(undf_w3)
    real(kind=r_tran), target, intent(in)    :: dry_mass_swapped(undf_w3)
    real(kind=r_tran),         intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_tran),         intent(in)    :: frac_dry_flux(undf_w2h)
    real(kind=r_tran),         intent(in)    :: min_val, dt
    logical(kind=l_def),       intent(in)    :: x_direction

    integer(kind=i_def)         :: local_dofs(2), dof_offset
    integer(kind=i_def)         :: rel_idx_min, rel_idx_max, ndep_half
    integer(kind=i_def)         :: rot_edge_dist_l, rot_edge_dist_r
    integer(kind=i_def)         :: rotated_panel_l, rotated_panel_r
    integer(kind=i_def)         :: edge_dist_local(-stencil_max:stencil_max)
    integer(kind=i_def)         :: edge_dist_l, edge_dist_r, pid_idx
    integer(kind=i_def)         :: rel_idx_k, idx_dep, w2h_idx, col_idx
    integer(kind=i_def)         :: k_low, k_high, k, j, int_cell_idx, recon_size
    integer(kind=i_def)         :: sign_disp(nlayers), edge_dist(nlayers)
    integer(kind=i_def)         :: edge_dist_dep(nlayers)
    integer(kind=i_def)         :: rel_dep_cell_idx(nlayers)
    integer(kind=i_def)         :: rel_dep_cell_idx_bounded(nlayers)
    integer(kind=i_def)         :: rel_idx(nlayers), loc_idx(nlayers)
    real(kind=r_tran)           :: displacement(nlayers), frac_dist(nlayers)
    real(kind=r_tran)           :: recon_field(nlayers)
    real(kind=r_tran)           :: field_edge_left(nlayers), field_edge_right(nlayers)
    real(kind=r_tran)           :: field_local(nlayers, 1+2*order)
    real(kind=r_tran)           :: int_cell_switch(nlayers)
    real(kind=r_tran)           :: swap_switch(nlayers), remap_switch(nlayers)
    real(kind=r_tran)           :: direction, inv_dt
    real(kind=r_tran), pointer  :: field_ptr(:)
    real(kind=r_tran), pointer  :: mass_ptr(:)

    inv_dt = 1.0_r_tran / dt
    if (x_direction) then
      local_dofs = (/ W, E /)
      direction = 1.0_r_tran
    else
      local_dofs = (/ S, N /)
      direction = -1.0_r_tran
    end if
    if (face == local_dofs(1)) then
      dof_offset = 0
    else
      dof_offset = 1
    end if

    ndep_half = (ndep - 1_i_def) / 2_i_def
    rel_idx_min = MAX(-stencil_extent_l, -ndep_half)
    rel_idx_max = MIN(stencil_extent_r, ndep_half)
    pid_idx = map_pid(1)
    rotated_panel_l = rotated_panel_neighbour(ipanel, local_dofs(1))
    rotated_panel_r = rotated_panel_neighbour(ipanel, local_dofs(2))
    if (ABS(rotated_panel_l) > 0) then
      rot_edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
    else
      rot_edge_dist_l = -FAR_AWAY
    end if
    if (ABS(rotated_panel_r) > 0) then
      rot_edge_dist_r = ABS(panel_edge_dist_r(pid_idx))
    else
      rot_edge_dist_r = FAR_AWAY
    end if

    do rel_idx_k = -stencil_extent_l, -1
      pid_idx = stencil_map_pid(rel_idx_k)
      if (rotated_panel_l == 1 .and. rel_idx_k < rot_edge_dist_l + 1) then
        edge_dist_l = -ABS(panel_edge_dist_swapped_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_r(pid_idx))
      else if (rotated_panel_l == -1 .and. rel_idx_k < rot_edge_dist_l + 1) then
        edge_dist_l = -ABS(panel_edge_dist_swapped_r(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_l(pid_idx))
      else
        edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_r(pid_idx))
      end if
      if (ABS(edge_dist_l) < ABS(edge_dist_r)) then
        edge_dist_local(rel_idx_k) = edge_dist_l
      else
        edge_dist_local(rel_idx_k) = edge_dist_r
      end if
    end do

    pid_idx = map_pid(1)
    edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
    edge_dist_r = ABS(panel_edge_dist_r(pid_idx))
    if (ABS(edge_dist_l) < ABS(edge_dist_r)) then
      edge_dist_local(0) = edge_dist_l
    else
      edge_dist_local(0) = edge_dist_r
    end if

    do rel_idx_k = 1, stencil_extent_r
      pid_idx = stencil_map_pid(rel_idx_k)
      if (rotated_panel_r == 1 .and. rel_idx_k > rot_edge_dist_r - 1) then
        edge_dist_l = -ABS(panel_edge_dist_swapped_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_r(pid_idx))
      else if (rotated_panel_r == -1 .and. rel_idx_k > rot_edge_dist_r - 1) then
        edge_dist_l = -ABS(panel_edge_dist_swapped_r(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_swapped_l(pid_idx))
      else
        edge_dist_l = -ABS(panel_edge_dist_l(pid_idx))
        edge_dist_r = ABS(panel_edge_dist_r(pid_idx))
      end if
      if (ABS(edge_dist_l) < ABS(edge_dist_r)) then
        edge_dist_local(rel_idx_k) = edge_dist_l
      else
        edge_dist_local(rel_idx_k) = edge_dist_r
      end if
    end do

    w2h_idx = map_w2h(face)
    idx_dep = map_depk(face) + ndep_half
    displacement(:) = direction * dep_dist(w2h_idx:w2h_idx+nlayers-1)
    frac_dist(:) = ABS(displacement(:) - REAL(INT(displacement(:), i_def), r_tran))
    sign_disp(:) = INT(SIGN(1.0_r_tran, displacement(:)))
    rel_dep_cell_idx(:) = dof_offset - INT(displacement(:), i_def) +           &
                          (1 - sign_disp(:)) / 2 - 1
    field_local = 0.0_r_tran
    flux(w2h_idx:w2h_idx+nlayers-1) = 0.0_r_tran

    do rel_idx_k = rel_idx_min, rel_idx_max
      k_low = dep_lowest_k(idx_dep + rel_idx_k)
      k_high = dep_highest_k(idx_dep + rel_idx_k)
      if (k_high < 0_i_def) cycle

      col_idx = stencil_map(rel_idx_k)
      if (rel_idx_k < rot_edge_dist_l + 1 .or.                                 &
          rel_idx_k > rot_edge_dist_r - 1) then
        if (SIGN(1, edge_dist_local(rel_idx_k)) /=                           &
            SIGN(1, edge_dist_local(0))) then
          field_ptr => remapped_field_swapped
        else
          field_ptr => field_swapped
        end if
        mass_ptr => dry_mass_swapped
      else
        if (SIGN(1, edge_dist_local(rel_idx_k)) /=                           &
            SIGN(1, edge_dist_local(0))) then
          field_ptr => remapped_field
        else
          field_ptr => field
        end if
        mass_ptr => dry_mass
      end if

      do k = k_low, k_high
        int_cell_idx = 1 + sign_disp(k+1) * (dof_offset - rel_idx_k - 1 +      &
                                             (1 - sign_disp(k+1)) / 2)
        int_cell_switch(k+1) = REAL(                                           &
            (1 + SIGN(1, ABS(INT(displacement(k+1), i_def)) - int_cell_idx)    &
              * SIGN(1, int_cell_idx - 1)) / 2, r_tran)
        if (int_cell_switch(k+1) == 0.0_r_tran) cycle

        flux(w2h_idx+k) = flux(w2h_idx+k) +                                   &
            int_cell_switch(k+1) * field_ptr(col_idx+k) * mass_ptr(col_idx+k)
      end do
    end do

    do rel_idx_k = -stencil_extent_l, stencil_extent_r
      if (ABS(edge_dist_local(rel_idx_k)) > order) then
        edge_dist_local(rel_idx_k) = FAR_AWAY
      end if
    end do

    rel_dep_cell_idx_bounded(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l, &
                                                            rel_dep_cell_idx(:)))
    do k = 1, nlayers
      edge_dist_dep(k) = edge_dist_local(rel_dep_cell_idx_bounded(k))
    end do

    recon_size = 1 + 2*order
    do j = 1, recon_size
      rel_idx(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l,               &
          rel_dep_cell_idx(:) + j - order - 1))
      swap_switch(:) = REAL(1 -                                              &
          (1 + SIGN(1, rel_idx - rot_edge_dist_l - 1)) *                      &
          (1 + SIGN(1, rot_edge_dist_r - rel_idx - 1)) / 4, r_tran)
      loc_idx(:) = (1 + sign_disp(:))/2*j +                                   &
                   (1 - sign_disp(:))/2*(recon_size - j + 1)

      if (j == order + 1) then
        do k = 1, nlayers
          col_idx = stencil_map(rel_idx(k))
          field_local(k, loc_idx(k)) =                                        &
              (1.0_r_tran - swap_switch(k))*field(col_idx+k-1) +              &
              swap_switch(k)*field_swapped(col_idx+k-1)
        end do
      else
        do k = 1, nlayers
          edge_dist(k) = edge_dist_local(rel_idx(k))
        end do
        remap_switch(:) = REAL(                                               &
            (1 - SIGN(1, edge_dist)*SIGN(1, edge_dist_dep))                   &
            * (1 + SIGN(1, ABS(edge_dist_dep) - ABS(edge_dist)))              &
            * (1 + SIGN(1, order - ABS(edge_dist_dep))) / 8, r_tran)
        do k = 1, nlayers
          col_idx = stencil_map(rel_idx(k))
          field_local(k, loc_idx(k)) =                                        &
              (1.0_r_tran - remap_switch(k)) * (                              &
                  (1.0_r_tran - swap_switch(k))*field(col_idx+k-1) +          &
                  swap_switch(k)*field_swapped(col_idx+k-1)) +               &
              remap_switch(k) * (                                             &
                  (1.0_r_tran - swap_switch(k))*remapped_field(col_idx+k-1) + &
                  swap_switch(k)*remapped_field_swapped(col_idx+k-1))
        end do
      end if
    end do

    select case (order)
    case (1)
      call nirvana_horizontal_edge(                                           &
          field_local, field_edge_left, field_edge_right, nlayers)
    case (2)
      call fourth_order_horizontal_edge(                                      &
          field_local, field_edge_left, field_edge_right, nlayers)
    end select

    if (order == 0) then
      recon_field(:) = field_local(:,1)
    else
      call monotonic_edge(field_local, monotone, min_val,                     &
                          field_edge_left, field_edge_right, order, 1, nlayers)
      call subgrid_quadratic_recon(recon_field, frac_dist, field_local,        &
          field_edge_left, field_edge_right, monotone, order, nlayers)
    end if

    flux(w2h_idx:w2h_idx+nlayers-1) = inv_dt * (                              &
        recon_field(:) * frac_dry_flux(w2h_idx:w2h_idx+nlayers-1)              &
        + direction * sign_disp(:) * flux(w2h_idx:w2h_idx+nlayers-1))

  end subroutine peregrin_flux_1d

end module peregrin_flux_kernel_mod
