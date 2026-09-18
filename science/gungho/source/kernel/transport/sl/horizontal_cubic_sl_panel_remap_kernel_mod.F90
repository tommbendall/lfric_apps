!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the advective increments in x and y at time n+1 using
!!          cubic semi-Lagrangian transport, using panel-edge remapped field
!!          values when a departure point crosses a cubed-sphere panel edge.
!> @details This is the panel-edge-aware counterpart of
!!          `horizontal_cubic_sl_kernel_mod.F90`: it solves the same
!!          one-dimensional advection equation in x and y using cubic
!!          interpolation, giving advective increments in both directions
!!          (this is the "outer" step of the COSMIC splitting, so the x
!!          increment works on the field previously advected in the
!!          y-direction, and vice versa), but for each level (since
!!          departure points can vary with height) it picks, per
!!          interpolation point, whether to use the normal-mesh `field` or
!!          the panel-edge-remapped `remapped_field` -- the remapped field
!!          is used only for those points whose relative stencil offset lies
!!          beyond this column's distance to the panel edge in that
!!          direction (i.e. the point would otherwise be looked up on the
!!          wrong panel).
!!
!!          When a stencil point crosses a *rotated* panel edge (i.e. one
!!          where the local x/y directions swap, as identified by
!!          `rotated_panel_neighbour` in `panel_edge_support_mod.F90`), the
!!          field/remapped_field pair from the *other* direction is used
!!          instead, exactly as in `horizontal_cubic_sl_sphere_kernel_mod.F90`.
!!          This gives, per point, up to four possible sources: `field`,
!!          `remapped_field`, `field_swapped` and `remapped_field_swapped`.
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest
!!       order.
!> @note This is a first draft (tracker step 3.5.2). The convention that
!!       `panel_edge_dist_*` equal to zero means the column is immediately
!!       adjacent to (i.e. on) the panel edge, and that a stencil point at
!!       relative offset `o` crosses the edge when `o` exceeds the local
!!       edge distance in that direction, has not yet been validated against
!!       the PEREGRIN design and should be reviewed.

module horizontal_cubic_sl_panel_remap_kernel_mod

  use argument_mod,          only: arg_type,                     &
                                   GH_FIELD, GH_REAL,            &
                                   CELL_COLUMN, GH_WRITE,        &
                                   GH_READ, GH_SCALAR,           &
                                   STENCIL, CROSS2D, GH_INTEGER, &
                                   ANY_DISCONTINUOUS_SPACE_1,    &
                                   ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,         only: r_tran, i_def, l_def, r_def
  use fs_continuity_mod,     only: W2H
  use kernel_mod,            only: kernel_type
  use reference_element_mod, only: W, E, S, N

  ! TODO: to remove
  use log_mod, only: log_event, log_level_debug, log_scratch_space

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_cubic_sl_panel_remap_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! increment_x
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! increment_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! remapped_field_x_in_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! remapped_field_y_in_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H),                       & ! dep_pts
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_id
        arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_edge_dist
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                              & ! monotone
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_cubic_sl_panel_remap_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: horizontal_cubic_sl_panel_remap_code
  public :: horizontal_cubic_sl_panel_remap_1d

contains

  !> @brief Compute advective transport in x and y directions using 1D
  !!        Semi-Lagrangian schemes, with a cubic reconstruction, using
  !!        panel-edge remapped field values where the stencil crosses a
  !!        cubed-sphere panel edge. This is the "outer" step of a COSMIC
  !!        splitting scheme.
  !> @param[in]     nlayers            Number of layers
  !> @param[in,out] increment_x        Advective increment in x direction
  !> @param[in,out] increment_y        Advective increment in y direction
  !> @param[in]     field_x            Field from x direction
  !> @param[in]     stencil_sizes_x    Sizes of the branches of the cross
  !!                                   stencil
  !> @param[in]     stencil_max_x      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_x      Dofmap for the field_x stencil
  !> @param[in]     field_y            Field from y direction
  !> @param[in]     stencil_sizes_y    Sizes of the branches of the cross
  !!                                   stencil
  !> @param[in]     stencil_max_y      Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_y      Dofmap for the field_y stencil
  !> @param[in]     remapped_field_x_in_y   field_x remapped in y, used for
  !!                                   stencil points that cross a panel edge
  !> @param[in]     stencil_sizes_rx   Sizes of the branches of the cross
  !!                                   stencil for remapped_field_x
  !> @param[in]     stencil_max_rx     Maximum size of a cross stencil branch
  !!                                   for remapped_field_x
  !> @param[in]     stencil_map_rx     Dofmap for the remapped_field_x stencil
  !> @param[in]     remapped_field_y_in_x  field_y remapped in X, used for
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
  !> @param[in]     monotone           Horizontal monotone option for cubic SL
  !> @param[in]     ndf_wf             Num of DoFs for field per cell
  !> @param[in]     undf_wf            Num of DoFs for this partition for
  !!                                   field
  !> @param[in]     map_wf             Map for Wf
  !> @param[in]     ndf_w2h            Num of DoFs for W2H per cell
  !> @param[in]     undf_w2h           Num of DoFs for this partition for
  !!                                   W2H
  !> @param[in]     map_w2h            Map for W2H
  !> @param[in]     ndf_2d             Num of DoFs for 2D fields per cell
  !> @param[in]     undf_2d            Num of DoFs in this partition for 2D
  !!                                   fields
  !> @param[in]     map_2d             Map for 2D fields
  subroutine horizontal_cubic_sl_panel_remap_code( nlayers,                &
                                                   increment_x,            &
                                                   increment_y,            &
                                                   field_x,                &
                                                   stencil_sizes_x,        &
                                                   stencil_max_x,          &
                                                   stencil_map_x,          &
                                                   field_y,                &
                                                   stencil_sizes_y,        &
                                                   stencil_max_y,          &
                                                   stencil_map_y,          &
                                                   remapped_field_x_in_y,  &
                                                   stencil_sizes_rx,       &
                                                   stencil_max_rx,         &
                                                   stencil_map_rx,         &
                                                   remapped_field_y_in_x,  &
                                                   stencil_sizes_ry,       &
                                                   stencil_max_ry,         &
                                                   stencil_map_ry,         &
                                                   dep_pts,                &
                                                   panel_id,               &
                                                   panel_edge_dist_W,      &
                                                   panel_edge_dist_E,      &
                                                   panel_edge_dist_S,      &
                                                   panel_edge_dist_N,      &
                                                   monotone,               &
                                                   ndf_wf,                 &
                                                   undf_wf,                &
                                                   map_wf,                 &
                                                   ndf_w2h,                &
                                                   undf_w2h,               &
                                                   map_w2h,                &
                                                   ndf_2d,                 &
                                                   undf_2d,                &
                                                   map_2d )

    use horizontal_cubic_sl_kernel_mod, only: horizontal_cubic_sl_1d

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_2d
    integer(kind=i_def), intent(in) :: ndf_2d
    integer(kind=i_def), intent(in) :: stencil_max_x
    integer(kind=i_def), intent(in) :: stencil_max_y
    integer(kind=i_def), intent(in) :: stencil_max_rx
    integer(kind=i_def), intent(in) :: stencil_max_ry
    integer(kind=i_def), intent(in) :: stencil_sizes_x(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_y(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_rx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_ry(4)
    integer(kind=i_def), intent(in) :: monotone

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_wf,stencil_max_x,4)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_wf,stencil_max_y,4)
    integer(kind=i_def), intent(in) :: stencil_map_rx(ndf_wf,stencil_max_rx,4)
    integer(kind=i_def), intent(in) :: stencil_map_ry(ndf_wf,stencil_max_ry,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: increment_x(undf_wf)
    real(kind=r_tran),   intent(inout) :: increment_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_x(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field_x_in_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field_y_in_x(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)
    real(kind=r_def),    intent(in)    :: panel_id(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_2d)
    integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_2d)

    ! Internal arguments
    integer(kind=i_def) :: i, ipanel
    integer(kind=i_def) :: edge_dist_E, edge_dist_W, edge_dist_S, edge_dist_N
    logical(kind=l_def) :: near_edge
    integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
    integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
    integer(kind=i_def) :: stencil_map_x_1d(-stencil_max_x:stencil_max_x)
    integer(kind=i_def) :: stencil_map_y_1d(-stencil_max_y:stencil_max_y)
    integer(kind=i_def) :: stencil_map_rx_1d(-stencil_max_rx:stencil_max_rx)
    integer(kind=i_def) :: stencil_map_ry_1d(-stencil_max_ry:stencil_max_ry)

    ipanel = INT(panel_id(map_2d(1)), i_def)
    edge_dist_E = panel_edge_dist_E(map_2d(1))
    edge_dist_W = panel_edge_dist_W(map_2d(1))
    edge_dist_S = panel_edge_dist_S(map_2d(1))
    edge_dist_N = panel_edge_dist_N(map_2d(1))

    ! Form X and Y 1D stencils
    stencil_extent_xl = stencil_sizes_x(1) - 1
    stencil_extent_xr = stencil_sizes_x(3) - 1
    stencil_extent_yl = stencil_sizes_y(2) - 1
    stencil_extent_yr = stencil_sizes_y(4) - 1

    do i = -stencil_extent_xl, 0
      stencil_map_x_1d(i) = stencil_map_x(1, 1-i, 1)
      stencil_map_rx_1d(i) = stencil_map_rx(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_x_1d(i) = stencil_map_x(1, i+1, 3)
      stencil_map_rx_1d(i) = stencil_map_rx(1, i+1, 3)
    end do

    do i = -stencil_extent_yl, 0
      stencil_map_y_1d(i) = stencil_map_y(1, 1-i, 2)
      stencil_map_ry_1d(i) = stencil_map_ry(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_y_1d(i) = stencil_map_y(1, i+1, 4)
      stencil_map_ry_1d(i) = stencil_map_ry(1, i+1, 4)
    end do

    ! X-calculation ------------------------------------------------------------
    near_edge = (                                                              &
        (ABS(edge_dist_W) <= MAX(stencil_extent_xl, stencil_extent_xr))        &
        .or. (ABS(edge_dist_E) <= MAX(stencil_extent_xl, stencil_extent_xr))   &
    )

    if (near_edge) then
      call horizontal_cubic_sl_panel_remap_1d( nlayers,               &
                                               .true.,                &
                                               increment_x,           &
                                               field_y,               &
                                               remapped_field_y_in_x, &
                                               field_x,               &
                                               remapped_field_x_in_y, &
                                               stencil_extent_xl,     &
                                               stencil_extent_xr,     &
                                               stencil_max_x,         &
                                               stencil_map_x_1d,      &
                                               stencil_map_rx_1d,     &
                                               edge_dist_W,           &
                                               edge_dist_E,           &
                                               ipanel,                &
                                               dep_pts,               &
                                               monotone,              &
                                               ndf_wf,                &
                                               undf_wf,               &
                                               map_wf,                &
                                               ndf_w2h,               &
                                               undf_w2h,              &
                                               map_w2h )
    else
      call horizontal_cubic_sl_1d( nlayers,           &
                                   .true.,            &
                                   increment_x,       &
                                   field_y,           &
                                   stencil_extent_xl, &
                                   stencil_extent_xr, &
                                   stencil_max_x,     &
                                   stencil_map_x_1d,  &
                                   dep_pts,           &
                                   monotone,          &
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
      call horizontal_cubic_sl_panel_remap_1d( nlayers,                &
                                               .false.,                &
                                               increment_y,            &
                                               field_x,                &
                                               remapped_field_x_in_y,  &
                                               field_y,                &
                                               remapped_field_y_in_x,  &
                                               stencil_extent_yl,      &
                                               stencil_extent_yr,      &
                                               stencil_max_y,          &
                                               stencil_map_y_1d,       &
                                               stencil_map_ry_1d,      &
                                               edge_dist_S,            &
                                               edge_dist_N,            &
                                               ipanel,                 &
                                               dep_pts,                &
                                               monotone,               &
                                               ndf_wf,                 &
                                               undf_wf,                &
                                               map_wf,                 &
                                               ndf_w2h,                &
                                               undf_w2h,               &
                                               map_w2h )
    else
      call horizontal_cubic_sl_1d( nlayers,           &
                                   .false.,           &
                                   increment_y,       &
                                   field_x,           &
                                   stencil_extent_yl, &
                                   stencil_extent_yr, &
                                   stencil_max_y,     &
                                   stencil_map_y_1d,  &
                                   dep_pts,           &
                                   monotone,          &
                                   ndf_wf,            &
                                   undf_wf,           &
                                   map_wf,            &
                                   ndf_w2h,           &
                                   undf_w2h,          &
                                   map_w2h )
    end if

  end subroutine horizontal_cubic_sl_panel_remap_code

! ============================================================================ !
! SINGLE UNDERLYING 1D ROUTINE
! ============================================================================ !

  !> @brief General 1D calculation of cubic Semi-Lagrangian advective
  !!        increment, using panel-edge remapped field values for stencil
  !!        points that cross a panel edge, and swapping to the other
  !!        direction's field/remapped_field pair for points that cross a
  !!        *rotated* panel edge (mirroring
  !!        `horizontal_cubic_sl_sphere_kernel_mod.F90`).
  !> @param[in]     field                  Field from this direction
  !> @param[in]     remapped_field         Panel-edge-remapped field for this
  !!                                       direction
  !> @param[in]     field_swapped          Field from the other direction,
  !!                                       used for points that cross a
  !!                                       rotated panel edge
  !> @param[in]     remapped_field_swapped Panel-edge-remapped field from the
  !!                                       other direction, used for points
  !!                                       that cross a rotated panel edge
  !!                                       and also cross the panel edge
  !> @param[in]     edge_dist_l  Distance (number of columns) from this
  !!                             column to the panel edge on the "left"
  !!                             (W for x-direction, S for y-direction)
  !> @param[in]     edge_dist_r  Distance (number of columns) from this
  !!                             column to the panel edge on the "right"
  !!                             (E for x-direction, N for y-direction)
  !> @param[in]     ipanel       ID of the mesh panel for this column
  subroutine horizontal_cubic_sl_panel_remap_1d( nlayers,                &
                                                 x_direction,            &
                                                 increment,              &
                                                 field,                  &
                                                 remapped_field,         &
                                                 field_swapped,          &
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
                                                 monotone,               &
                                                 ndf_wf,                 &
                                                 undf_wf,                &
                                                 map_wf,                 &
                                                 ndf_w2h,                &
                                                 undf_w2h,               &
                                                 map_w2h )

    use panel_edge_support_mod,         only: rotated_panel_neighbour, FAR_AWAY
    use transport_enumerated_types_mod, only: monotone_strict,                 &
                                              monotone_relaxed,                &
                                              monotone_positive

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
    integer(kind=i_def), intent(in) :: monotone
    logical(kind=l_def), intent(in) :: x_direction

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: stencil_map(-stencil_max:stencil_max)
    integer(kind=i_def), intent(in) :: stencil_map_r(-stencil_max:stencil_max)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: increment(undf_wf)
    real(kind=r_tran),   intent(in)    :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_swapped(undf_wf)
    real(kind=r_tran),   intent(in)    :: remapped_field_swapped(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)

    ! Local arrays
    integer(kind=i_def) :: int_disp(nlayers+ndf_wf-1)
    integer(kind=i_def) :: sign_disp(nlayers+ndf_wf-1)
    integer(kind=i_def) :: rel_idx_hi_p(nlayers+ndf_wf-1)
    integer(kind=i_def) :: rel_idx(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: displacement(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: field_out(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: field_local(nlayers+ndf_wf-1,4)
    real(kind=r_tran)   :: q_max(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: q_min(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: xx(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: remap_switch(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: rotated_switch(nlayers+ndf_wf-1)

    ! Local scalars
    integer(kind=i_def) :: j, k, nl, f_idx
    integer(kind=i_def) :: w2h_df_l, w2h_df_r
    real(kind=r_tran)   :: direction
    integer(kind=i_def) :: edge_dist_l, edge_dist_r
    integer(kind=i_def) :: rot_edge_dist_l, rot_edge_dist_r
    integer(kind=i_def) :: rotated_panel_l, rotated_panel_r

    ! Cubic interpolation weights
    real(kind=r_tran), parameter :: x0 = 0.0_r_tran
    real(kind=r_tran), parameter :: x1 = 1.0_r_tran
    real(kind=r_tran), parameter :: x2 = 2.0_r_tran
    real(kind=r_tran), parameter :: x3 = 3.0_r_tran
    real(kind=r_tran), parameter :: den0 = 1.0_r_tran/((x0-x1)*(x0-x2)*(x0-x3))
    real(kind=r_tran), parameter :: den1 = 1.0_r_tran/((x1-x0)*(x1-x2)*(x1-x3))
    real(kind=r_tran), parameter :: den2 = 1.0_r_tran/((x2-x0)*(x2-x1)*(x2-x3))
    real(kind=r_tran), parameter :: den3 = 1.0_r_tran/((x3-x0)*(x3-x1)*(x3-x2))

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

    ! ======================================================================== !
    ! Extract departure info
    ! ======================================================================== !

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
    xx(:) = 1.0_r_tran + ABS(displacement(:) - REAL(int_disp, r_tran))
    sign_disp(:) = INT(SIGN(1.0_r_tran, displacement(:)))

    ! The relative index of the most downwind cell to use in the stencil
    rel_idx_hi_p(:) = - 2*sign_disp(:) - int_disp(:)

    ! ======================================================================== !
    ! Populate local arrays for interpolation. Per point: if the point crosses
    ! a rotated panel edge, use the field/remapped_field pair from the other
    ! (swapped) direction; otherwise use the pair from this direction. Within
    ! whichever pair is selected, use the panel-edge remapped field for points
    ! that cross the (possibly non-rotated) panel edge, and the normal field
    ! otherwise.
    ! ======================================================================== !

    ! Loop over points to use in reconstruction
    do j = 1, 4
      ! If this column has idx 0, find relative index for the column of the
      ! departure cell, between -stencil_extent_l and stencil_extent_r, e.g.
      ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
      rel_idx(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l,                &
          rel_idx_hi_p(:) + (4 - j)*sign_disp(:)                               &
      ))

      ! Determine whether this column has crossed a rotated panel edge
      ! This factor is 1 if the column has crossed a rotated panel edge, and 0
      ! if it has not. It is made from (1 - switch_l*switch_r) where:
      ! switch_l is 0 when on a rotated left panel, 1 otherwise
      ! switch_r is 0 when on a rotated right panel, 1 otherwise
      rotated_switch(:) = REAL(1                                               &
          - (1 + SIGN(1, rel_idx - rot_edge_dist_l - 1))                       &
          ! 0 if crossed to rotated panel, 1 if not
          * (1 + SIGN(1, rot_edge_dist_r - rel_idx - 1)) / 4, r_tran           &
      )

      ! Determine whether to use the field that has been remapped to the
      ! other panel, instead of the field that has not been remapped
      ! This factor is 1 if the column has crossed a panel edge, and 0
      ! if it has not. It is made from (1 - switch_l*switch_r) where:
      ! switch_l is 0 when on a left panel, 1 otherwise
      ! switch_r is 0 when on a right panel, 1 otherwise
      remap_switch(:) = REAL(1                                                 &
          - (1 + SIGN(1, rel_idx - edge_dist_l - 1))                           &
          ! 0 if crossed to rotated panel, 1 if not
          * (1 + SIGN(1, edge_dist_r - rel_idx - 1)) / 4, r_tran               &
      )

      do k = 1, nl
        f_idx = stencil_map(rel_idx(k))
        field_local(k, j) = (                                                  &
          (1.0_r_tran - remap_switch(k))*(                                     &
            (1.0_r_tran - rotated_switch(k))*field(f_idx+k-1)                  &
            + rotated_switch(k)*field_swapped(f_idx+k-1)                       &
          ) + remap_switch(k)*(                                                &
            (1.0_r_tran - rotated_switch(k))*remapped_field(f_idx+k-1)         &
            + rotated_switch(k)*remapped_field_swapped(f_idx+k-1)              &
          )                                                                    &
        )
        ! TODO: to remove
        if (field_local(k,j) < -100.0_r_tran) then
          write(log_scratch_space, *) 'CUBIC REMAP BAD VALUE: ', field_local(k,j), map_wf(1), &
            j, rel_idx(k), edge_dist_l, edge_dist_r, rot_edge_dist_l, rot_edge_dist_r, &
            rotated_switch(k), remap_switch(k)
          call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
        end if
      end do
    end do

    ! ======================================================================== !
    ! Perform cubic interpolation
    ! ======================================================================== !

    field_out = (                                                              &
        (xx(:)-x1) * (xx(:)-x2) * (xx(:)-x3) * den0 * field_local(:,1)         &
        + (xx(:)-x0) * (xx(:)-x2) * (xx(:)-x3) * den1 * field_local(:,2)       &
        + (xx(:)-x0) * (xx(:)-x1) * (xx(:)-x3) * den2 * field_local(:,3)       &
        + (xx(:)-x0) * (xx(:)-x1) * (xx(:)-x2) * den3 * field_local(:,4)       &
    )

    select case (map_wf(1))
    case (50, 66, 82, 98)
      write(log_scratch_space, *) 'SL CUBIC REMAP', map_wf(1), direction, &
        field(map_wf(1)), field_out(1), displacement(1)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end select

    ! ======================================================================== !
    ! Apply monotonicity constraints
    ! ======================================================================== !

    select case (monotone)
    case (monotone_strict)
      ! Bound field by immediately neighbouring values
      q_min(:) = MIN(field_local(:,2), field_local(:,3))
      q_max(:) = MAX(field_local(:,2), field_local(:,3))
      field_out(:) = MIN(q_max(:), MAX(field_out(:), q_min(:)))

    case (monotone_relaxed)
      ! Bound field by all values in the stencil
      q_min(:) = MIN(                                                          &
          field_local(:,1), field_local(:,2),                                  &
          field_local(:,3), field_local(:,4)                                   &
      )
      q_max(:) = MAX(                                                          &
          field_local(:,1), field_local(:,2),                                  &
          field_local(:,3), field_local(:,4)                                   &
      )
      field_out(:) = MIN(q_max(:), MAX(field_out(:), q_min(:)))

    case (monotone_positive)
      ! Just make sure field out is positive
      field_out(:) = MAX(field_out(:), 0.0_r_tran)

    end select

    ! ======================================================================== !
    ! Compute increment
    ! ======================================================================== !

    increment(map_wf(1) : map_wf(1)+nl-1) = (                                  &
      field_out(:) - field(map_wf(1) : map_wf(1)+nl-1)                         &
    )

  end subroutine horizontal_cubic_sl_panel_remap_1d

end module horizontal_cubic_sl_panel_remap_kernel_mod
