!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the fields in x and y at time n+1 using linear
!!          semi-Lagrangian transport, using panel-edge remapped field values
!!          when a departure point crosses a cubed-sphere panel edge.
!> @details This is the panel-edge-aware counterpart of
!!          `horizontal_linear_sl_kernel_mod.F90`: it solves the same
!!          one-dimensional advection equation in x and y using linear
!!          interpolation, but for each level (since departure points can
!!          vary with height) it picks, per interpolation point, whether to
!!          use the normal-mesh `field` or the panel-edge-remapped
!!          `remapped_field` -- the remapped field is used only for those
!!          points whose relative stencil offset lies beyond this column's
!!          distance to the panel edge in that direction (i.e. the point
!!          would otherwise be looked up on the wrong panel).
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest
!!       order.
!> @note This is a first draft (tracker step 3.5.2). The convention that
!!       `panel_edge_dist_*` equal to zero means the column is immediately
!!       adjacent to (i.e. on) the panel edge, and that a stencil point at
!!       relative offset `o` crosses the edge when `o` exceeds the local
!!       edge distance in that direction, has not yet been validated against
!!       the PEREGRIN design and should be reviewed.

module horizontal_linear_sl_panel_remap_kernel_mod

  use argument_mod,          only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   CELL_COLUMN, GH_WRITE,     &
                                   GH_READ, GH_SCALAR,        &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_3, &
                                   STENCIL, CROSS2D, GH_INTEGER
  use constants_mod,         only: i_def, r_tran, l_def, r_def
  use fs_continuity_mod,     only: W2H
  use kernel_mod,            only: kernel_type
  use reference_element_mod, only: W, E, S, N

  ! TODO: to remove
  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_DEBUG

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_linear_sl_panel_remap_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field_out_x
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field_out_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! remapped_field_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! remapped_field_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H),                       & ! dep_pts
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_id
        arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  & ! panel_edge_dist
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_linear_sl_panel_remap_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: horizontal_linear_sl_panel_remap_code
  public :: horizontal_linear_sl_panel_remap_1d

contains

  !> @brief Compute advective transport in x and y directions using 1D
  !!        Semi-Lagrangian schemes, with a linear reconstruction, using
  !!        panel-edge remapped field values where the stencil crosses a
  !!        cubed-sphere panel edge. This is the "inner" step of a COSMIC
  !!        splitting scheme.
  !> @param[in]     nlayers            Number of layers
  !> @param[in,out] field_x            Field at time n+1 in x direction
  !> @param[in,out] field_y            Field at time n+1 in y direction
  !> @param[in]     field              Field to transport
  !> @param[in]     stencil_sizes      Sizes of the branches of the cross
  !!                                   stencil
  !> @param[in]     stencil_max        Maximum size of a cross stencil branch
  !> @param[in]     stencil_map        Dofmap for the field stencil
  !> @param[in]     remapped_field_x   Panel-edge-remapped field_x, used for
  !!                                   stencil points that cross a panel edge
  !> @param[in]     stencil_sizes_rx   Sizes of the branches of the cross
  !!                                   stencil for remapped_field_x
  !> @param[in]     stencil_max_rx     Maximum size of a cross stencil branch
  !!                                   for remapped_field_x
  !> @param[in]     stencil_map_rx     Dofmap for the remapped_field_x stencil
  !> @param[in]     remapped_field_y   Panel-edge-remapped field_y, used for
  !!                                   stencil points that cross a panel edge
  !> @param[in]     stencil_sizes_ry   Sizes of the branches of the cross
  !!                                   stencil for remapped_field_y
  !> @param[in]     stencil_max_ry     Maximum size of a cross stencil branch
  !!                                   for remapped_field_y
  !> @param[in]     stencil_map_ry     Dofmap for the remapped_field_y stencil
  !> @param[in]     dep_pts            Departure points
  !> @param[in]     panel_id           Field containing IDs of mesh panels
  !> @param[in]     panel_edge_dist_W  2D field: distance from this column to
  !!                                   the panel edge to the West
  !> @param[in]     panel_edge_dist_E  2D field: distance from this column to
  !!                                   the panel edge to the East
  !> @param[in]     panel_edge_dist_S  2D field: distance from this column to
  !!                                   the panel edge to the South
  !> @param[in]     panel_edge_dist_N  2D field: distance from this column to
  !!                                   the panel edge to the North
  !> @param[in]     ndf_wf             Num of DoFs for field per cell
  !> @param[in]     undf_wf            Num of DoFs in this partition for field
  !> @param[in]     map_wf             Map for field
  !> @param[in]     ndf_w2h            Num of DoFs for W2H per cell
  !> @param[in]     undf_w2h           Num of DoFs in this partition for W2H
  !> @param[in]     map_w2h            Map for W2H
  !> @param[in]     ndf_2d             Num of DoFs for 2D fields per cell
  !> @param[in]     undf_2d            Num of DoFs in this partition for 2D
  !!                                   fields
  !> @param[in]     map_2d             Map for 2D fields
  subroutine horizontal_linear_sl_panel_remap_code( nlayers,             &
                                                    field_x,             &
                                                    field_y,             &
                                                    field,               &
                                                    stencil_sizes,       &
                                                    stencil_max,         &
                                                    stencil_map,         &
                                                    remapped_field_x,    &
                                                    stencil_sizes_rx,    &
                                                    stencil_max_rx,      &
                                                    stencil_map_rx,      &
                                                    remapped_field_y,    &
                                                    stencil_sizes_ry,    &
                                                    stencil_max_ry,      &
                                                    stencil_map_ry,      &
                                                    dep_pts,             &
                                                    panel_id,            &
                                                    panel_edge_dist_W,   &
                                                    panel_edge_dist_E,   &
                                                    panel_edge_dist_S,   &
                                                    panel_edge_dist_N,   &
                                                    ndf_wf,              &
                                                    undf_wf,             &
                                                    map_wf,              &
                                                    ndf_w2h,             &
                                                    undf_w2h,            &
                                                    map_w2h,             &
                                                    ndf_2d,              &
                                                    undf_2d,             &
                                                    map_2d )

    use horizontal_linear_sl_kernel_mod, only: horizontal_linear_sl_1d

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_2d
    integer(kind=i_def), intent(in) :: ndf_2d
    integer(kind=i_def), intent(in) :: stencil_max
    integer(kind=i_def), intent(in) :: stencil_max_rx
    integer(kind=i_def), intent(in) :: stencil_max_ry
    integer(kind=i_def), intent(in) :: stencil_sizes(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_rx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_ry(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
    integer(kind=i_def), intent(in) :: stencil_map(ndf_wf,stencil_max,4)
    integer(kind=i_def), intent(in) :: stencil_map_rx(ndf_wf,stencil_max_rx,4)
    integer(kind=i_def), intent(in) :: stencil_map_ry(ndf_wf,stencil_max_ry,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: field_x(undf_wf)
    real(kind=r_tran),   intent(inout) :: field_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field_x(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)
    real(kind=r_def),    intent(in)    :: panel_id(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_2d)

    ! Internal arguments
    integer(kind=i_def) :: i, ipanel
    integer(kind=i_def) :: edge_dist_E, edge_dist_W, edge_dist_S, edge_dist_N
    integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
    integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
    integer(kind=i_def) :: stencil_map_x_1d(-stencil_max:stencil_max)
    integer(kind=i_def) :: stencil_map_y_1d(-stencil_max:stencil_max)
    integer(kind=i_def) :: stencil_map_rx_1d(-stencil_max_rx:stencil_max_rx)
    integer(kind=i_def) :: stencil_map_ry_1d(-stencil_max_ry:stencil_max_ry)
    logical(kind=l_def) :: near_edge

    ipanel = INT(panel_id(map_2d(1)), i_def)
    edge_dist_E = panel_edge_dist_E(map_2d(1))
    edge_dist_W = panel_edge_dist_W(map_2d(1))
    edge_dist_S = panel_edge_dist_S(map_2d(1))
    edge_dist_N = panel_edge_dist_N(map_2d(1))

    ! Form X and Y 1D stencils (normal field)
    stencil_extent_xl = stencil_sizes(1) - 1
    stencil_extent_xr = stencil_sizes(3) - 1
    stencil_extent_yl = stencil_sizes(2) - 1
    stencil_extent_yr = stencil_sizes(4) - 1

    do i = -stencil_extent_xl, 0
      stencil_map_x_1d(i) = stencil_map(1, 1-i, 1)
      stencil_map_rx_1d(i) = stencil_map_rx(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_x_1d(i) = stencil_map(1, i+1, 3)
      stencil_map_rx_1d(i) = stencil_map_rx(1, i+1, 3)
    end do

    do i = -stencil_extent_yl, 0
      stencil_map_y_1d(i) = stencil_map(1, 1-i, 2)
      stencil_map_ry_1d(i) = stencil_map_ry(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_y_1d(i) = stencil_map(1, i+1, 4)
      stencil_map_ry_1d(i) = stencil_map_ry(1, i+1, 4)
    end do

    ! X-calculation ------------------------------------------------------------
    near_edge = (                                                              &
        (ABS(edge_dist_W) <= MAX(stencil_extent_xl, stencil_extent_xr))        &
        .or. (ABS(edge_dist_E) <= MAX(stencil_extent_xl, stencil_extent_xr))   &
    )

    if (near_edge) then
      call horizontal_linear_sl_panel_remap_1d( nlayers,           &
                                                .true.,            &
                                                field_x,           &
                                                field,             &
                                                remapped_field_x,  &
                                                remapped_field_y,  &
                                                stencil_extent_xl, &
                                                stencil_extent_xr, &
                                                stencil_max,       &
                                                stencil_map_x_1d,  &
                                                stencil_map_rx_1d, &
                                                edge_dist_W,       &
                                                edge_dist_E,       &
                                                ipanel,            &
                                                dep_pts,           &
                                                ndf_wf,            &
                                                undf_wf,           &
                                                map_wf,            &
                                                ndf_w2h,           &
                                                undf_w2h,          &
                                                map_w2h )
    else
      call horizontal_linear_sl_1d( nlayers,           &
                                    .true.,            &
                                    field_x,           &
                                    field,             &
                                    stencil_extent_xl, &
                                    stencil_extent_xr, &
                                    stencil_max,       &
                                    stencil_map_x_1d,  &
                                    dep_pts,           &
                                    ndf_wf,            &
                                    undf_wf,           &
                                    map_wf,            &
                                    ndf_w2h,           &
                                    undf_w2h,          &
                                    map_w2h )
    end if

    ! Y-calculation ------------------------------------------------------------
    near_edge = (                                                              &
        (ABS(edge_dist_S) <= MAX(stencil_extent_yl, stencil_extent_yr))        &
        .or. (ABS(edge_dist_N) <= MAX(stencil_extent_yl, stencil_extent_yr))   &
    )

    if (near_edge) then
      call horizontal_linear_sl_panel_remap_1d( nlayers,           &
                                                .false.,           &
                                                field_y,           &
                                                field,             &
                                                remapped_field_y,  &
                                                remapped_field_x,  &
                                                stencil_extent_yl, &
                                                stencil_extent_yr, &
                                                stencil_max,       &
                                                stencil_map_y_1d,  &
                                                stencil_map_ry_1d, &
                                                edge_dist_S,       &
                                                edge_dist_N,       &
                                                ipanel,            &
                                                dep_pts,           &
                                                ndf_wf,            &
                                                undf_wf,           &
                                                map_wf,            &
                                                ndf_w2h,           &
                                                undf_w2h,          &
                                                map_w2h )
    else
      call horizontal_linear_sl_1d( nlayers,           &
                                    .false.,           &
                                    field_y,           &
                                    field,             &
                                    stencil_extent_yl, &
                                    stencil_extent_yr, &
                                    stencil_max,       &
                                    stencil_map_y_1d,  &
                                    dep_pts,           &
                                    ndf_wf,            &
                                    undf_wf,           &
                                    map_wf,            &
                                    ndf_w2h,           &
                                    undf_w2h,          &
                                    map_w2h )
    end if

  end subroutine horizontal_linear_sl_panel_remap_code

! ============================================================================ !
! SINGLE UNDERLYING 1D ROUTINE
! ============================================================================ !

  !> @brief General 1D calculation of linear Semi-Lagrangian advected field,
  !!        using panel-edge remapped field values for stencil points that
  !!        cross a panel edge.
  !> @param[in]     edge_dist_l  Distance (number of columns) from this
  !!                             column to the panel edge on the "left"
  !!                             (W for x-direction, S for y-direction)
  !> @param[in]     edge_dist_r  Distance (number of columns) from this
  !!                             column to the panel edge on the "right"
  !!                             (E for x-direction, N for y-direction)
  subroutine horizontal_linear_sl_panel_remap_1d( nlayers,                &
                                                  x_direction,            &
                                                  field_out,              &
                                                  field,                  &
                                                  remapped_field,         &
                                                  remapped_field_swapped, &
                                                  stencil_extent_l,       &
                                                  stencil_extent_r,       &
                                                  stencil_max,            &
                                                  stencil_map,            &
                                                  stencil_map_r,          &
                                                  edge_dist_l_in,         &
                                                  edge_dist_r_in,         &
                                                  ipanel,                 &
                                                  dep_pts,                &
                                                  ndf_wf,                 &
                                                  undf_wf,                &
                                                  map_wf,                 &
                                                  ndf_w2h,                &
                                                  undf_w2h,               &
                                                  map_w2h )

    use panel_edge_support_mod,         only: rotated_panel_neighbour, FAR_AWAY

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_extent_l
    integer(kind=i_def), intent(in) :: stencil_extent_r
    integer(kind=i_def), intent(in) :: stencil_max
    integer(kind=i_def), intent(in) :: edge_dist_l_in
    integer(kind=i_def), intent(in) :: edge_dist_r_in
    integer(kind=i_def), intent(in) :: ipanel
    logical(kind=l_def), intent(in) :: x_direction

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: stencil_map(-stencil_max:stencil_max)
    integer(kind=i_def), intent(in) :: stencil_map_r(-stencil_max:stencil_max)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: field_out(undf_wf)
    real(kind=r_tran),   intent(in)    :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field_swapped(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)

    ! Local arrays
    integer(kind=i_def) :: int_disp(nlayers+ndf_wf-1)
    integer(kind=i_def) :: sign_disp(nlayers+ndf_wf-1)
    integer(kind=i_def) :: rel_idx_hi(nlayers+ndf_wf-1)
    integer(kind=i_def) :: rel_idx(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: displacement(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: field_local(nlayers+ndf_wf-1,2)
    real(kind=r_tran)   :: xx(nlayers+ndf_wf-1)
    logical(kind=l_def) :: swapped

    ! Local scalars
    integer(kind=i_def) :: j, k, nl
    integer(kind=i_def) :: w2h_df_l, w2h_df_r, f_idx
    real(kind=r_tran)   :: direction
    integer(kind=i_def) :: edge_dist_l, edge_dist_r
    integer(kind=i_def) :: rot_edge_dist_l, rot_edge_dist_r
    integer(kind=i_def) :: rotated_panel_l, rotated_panel_r

    ! nl = nlayers      for w3
    !    = nlayers+1    for wtheta
    nl = nlayers + ndf_wf - 1

    if (x_direction) then
      w2h_df_l = map_w2h(W)
      w2h_df_r = map_w2h(E)
      direction = 1.0_r_tran
      rotated_panel_l = rotated_panel_neighbour(ipanel, W)
      rotated_panel_r = rotated_panel_neighbour(ipanel, E)
    else
      ! y-direction
      w2h_df_l = map_w2h(S)
      w2h_df_r = map_w2h(N)
      direction = -1.0_r_tran
      rotated_panel_l = rotated_panel_neighbour(ipanel, S)
      rotated_panel_r = rotated_panel_neighbour(ipanel, N)
    end if

    ! ======================================================================== !
    ! Extract departure info
    ! ======================================================================== !

    ! Pre-determine aspects of crossing a rotated panel boundary, as in
    ! horizontal_cubic_sl_sphere_kernel_mod.F90
    if (ABS(rotated_panel_l) > 0) then
      rot_edge_dist_l = -ABS(edge_dist_l_in)
    else
      rot_edge_dist_l = -FAR_AWAY
    end if
    edge_dist_l = -ABS(edge_dist_l_in)
    if (ABS(rotated_panel_r) > 0) then
      rot_edge_dist_r = ABS(edge_dist_r_in)
    else
      rot_edge_dist_r = FAR_AWAY
    end if
    edge_dist_r = ABS(edge_dist_r_in)

    if (ndf_wf == 1) then
      ! Advecting W3 field: average the dep distances from this cell's faces
      displacement(:) = 0.5_r_tran * direction * (                             &
        dep_pts(w2h_df_l : w2h_df_l+nl-1)                                      &
        + dep_pts(w2h_df_r : w2h_df_r+nl-1)                                    &
      )
    else
      ! Advecting Wtheta field:
      ! In top and bottom layers, take the dep distances for top/bottom layer
      displacement(1) = 0.5_r_tran * direction * (                             &
        dep_pts(w2h_df_l) + dep_pts(w2h_df_r)                                  &
      )
      if (nlayers > 1) then
        ! NB: nl = nlayers + 1
        displacement(2:nl-1) = 0.25_r_tran * direction * (                     &
          dep_pts(w2h_df_l : w2h_df_l+nl-3)                                    &
          + dep_pts(w2h_df_l+1 : w2h_df_l+nl-2)                                &
          + dep_pts(w2h_df_r : w2h_df_r+nl-3)                                  &
          + dep_pts(w2h_df_r+1 : w2h_df_r+nl-2)                                &
        )
      end if
      ! Top layer
      displacement(nl) = 0.5_r_tran * direction * (                            &
        dep_pts(w2h_df_l+nl-2) + dep_pts(w2h_df_r+nl-2)                        &
      )
    end if

    int_disp(:) = INT(displacement(:), i_def)
    xx(:) = ABS(displacement(:) - REAL(int_disp, r_tran))
    sign_disp(:) = INT(SIGN(1.0_r_tran, displacement(:)))

    ! The relative index of the furthest cell to use in the stencil
    rel_idx_hi(:) = - sign_disp(:) - int_disp(:)

    ! ======================================================================== !
    ! Populate local arrays for interpolation, using the panel-edge remapped
    ! field for any stencil point that crosses the panel edge (i.e. whose
    ! relative offset exceeds this column's distance to the edge in that
    ! direction), and the normal field otherwise.
    ! ======================================================================== !

    ! Loop over points to use in reconstruction
    do j = 1, 2
      ! departure cell, between -stencil_extent_l and stencil_extent_r, e.g.
      ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
      rel_idx(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l,                &
          rel_idx_hi(:) + (2 - j)*sign_disp(:)                                 &
      ))

      ! Loop over layers
      do k = 1, nl
        swapped = (                                                            &
          rel_idx(k) <= rot_edge_dist_l .or. rel_idx(k) >= rot_edge_dist_r     &
        )
        f_idx = stencil_map(rel_idx(k))+k-1
        if (swapped) then
          if (rel_idx(k) <= edge_dist_l .or. rel_idx(k) >= edge_dist_r) then
            field_local(k,j) = remapped_field_swapped(stencil_map_r(rel_idx(k))+k-1)
          else
            field_local(k,j) = field(f_idx)
          end if
        else
          if (rel_idx(k) <= edge_dist_l .or. rel_idx(k) >= edge_dist_r) then
            field_local(k,j) = remapped_field(stencil_map_r(rel_idx(k))+k-1)
          else
            field_local(k,j) = field(f_idx)
          end if
        end if

        ! TODO: to remove
        if (field_local(k,j) < -100.0_r_tran) then
          write(log_scratch_space, *) 'LINEAR REMAP BAD VALUE: ', field_local(k,j), map_wf(1), &
            edge_dist_l, edge_dist_r, swapped, rel_idx(k)
          call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
        end if
      end do
    end do

    ! ======================================================================== !
    ! Perform linear interpolation
    ! ======================================================================== !

    ! Linear interpolation in x
    ! interp = (x-x1) / (x0-x1) f(x0) + (x-x0) / (x1-x0) f(x1)
    ! Set x0 = 0, x1 = 1, and 0 <= x <= 1
    ! interp = - (x-1) f(0) + x f(1)
    field_out(map_wf(1) : map_wf(1)+nl-1) = (                                  &
      -(xx(:)-1.0_r_tran) * field_local(:,1) + xx(:) * field_local(:,2)        &
    )

    select case (map_wf(1))
    case (862, 863, 864)
      write(log_scratch_space, *) 'SL LINEAR REMAP:', map_wf(1), direction, &
        field(map_wf(1)), field_out(map_wf(1)), displacement(1)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end select

  end subroutine horizontal_linear_sl_panel_remap_1d

end module horizontal_linear_sl_panel_remap_kernel_mod