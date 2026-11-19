!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the advective increments in x and y at time n+1 for the
!!          PEREGRIN treatment of cubic semi-Lagrangian transport near
!!          cubed-sphere panel edges.
!> @details This module has no reconstruction logic of its own: it dispatches,
!!          per column and per direction (x using W/E faces, y using S/N
!!          faces), to the existing 1D SL routines:
!!            (a) not near a panel edge (this column's distance to both the
!!                left and right panel edges in that direction is greater
!!                than the stencil extent, so no stencil point could reach
!!                the edge): the plain routine `horizontal_cubic_sl_1d`
!!                (`horizontal_cubic_sl_kernel_mod.F90`), using the normal
!!                field only.
!!            (b) near a panel edge (this column's distance to the left or
!!                right panel edge is within the stencil extent, so some
!!                stencil point may cross the edge): the panel-remap routine
!!                `horizontal_cubic_sl_panel_remap_1d`
!!                (`horizontal_cubic_sl_panel_remap_kernel_mod.F90`), using
!!                `field_px`/`field_py` (the PEREGRIN-transported fields
!!                computed by `peregrin_eave_sl_kernel_mod.F90` in the inner
!!                COSMIC-splitting step) for whichever individual stencil
!!                points (per level) have a departure distance that crosses
!!                the edge, and the normal field otherwise.
!!          This is the PEREGRIN counterpart to `peregrin_flux_kernel_mod.F90`
!!          for FFSL, containing the "outer" COSMIC-splitting step kernel
!!          `peregrin_cubic_sl_kernel_type`, analogous to
!!          `horizontal_cubic_sl_kernel_type`/`horizontal_cubic_sl_sphere_kernel_type`.
!!
!> @note This is a first draft (tracker step 3.5.3). Unlike the FFSL
!!       treatment in `peregrin_flux_kernel_mod.F90`, no distinction is made
!!       here between "upwind"/"downwind" columns: the panel-remap routine
!!       already falls back to the normal field, per level, for any stencil
!!       point that does not itself cross the edge, so a single "near
!!       edge"/"not near edge" dispatch (based only on this column's
!!       distance to the edge versus the stencil extent) is sufficient. This
!!       has not yet been validated against the PEREGRIN design and should
!!       be reviewed alongside 3.5.1/3.5.2.
!> @note This kernel only works when field is a W3/Wtheta field at lowest
!!       order.

module peregrin_cubic_sl_kernel_mod

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

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the "outer"/cubic PEREGRIN SL kernel. Contains
  !! the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: peregrin_cubic_sl_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! increment_x
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! increment_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_px
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_py
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H),                       & ! dep_pts
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_id
        arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_edge_dist
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                              & ! monotone
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: peregrin_cubic_sl_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: peregrin_cubic_sl_code

contains

  !> @brief Compute the PEREGRIN "outer"/cubic-reconstruction SL advective
  !!        increments in x and y directions, dispatching per direction to
  !!        the plain or panel-remap 1D routines depending on proximity to a
  !!        panel edge.
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
  !> @param[in]     field_px           PEREGRIN-transported field_x (from
  !!                                   `peregrin_eave_sl_kernel_mod.F90`),
  !!                                   used for stencil points whose
  !!                                   departure distance crosses a panel
  !!                                   edge
  !> @param[in]     stencil_sizes_px   Sizes of the branches of the cross
  !!                                   stencil for field_px
  !> @param[in]     stencil_max_px     Maximum size of a cross stencil branch
  !!                                   for field_px
  !> @param[in]     stencil_map_px     Dofmap for the field_px stencil
  !> @param[in]     field_py           PEREGRIN-transported field_y (from
  !!                                   `peregrin_eave_sl_kernel_mod.F90`),
  !!                                   used for stencil points whose
  !!                                   departure distance crosses a panel
  !!                                   edge
  !> @param[in]     stencil_sizes_py   Sizes of the branches of the cross
  !!                                   stencil for field_py
  !> @param[in]     stencil_max_py     Maximum size of a cross stencil branch
  !!                                   for field_py
  !> @param[in]     stencil_map_py     Dofmap for the field_py stencil
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
  subroutine peregrin_cubic_sl_code( nlayers,           &
                                     increment_x,       &
                                     increment_y,       &
                                     field_x,           &
                                     stencil_sizes_x,   &
                                     stencil_max_x,     &
                                     stencil_map_x,     &
                                     field_y,           &
                                     stencil_sizes_y,   &
                                     stencil_max_y,     &
                                     stencil_map_y,     &
                                     field_px,          &
                                     stencil_sizes_px,  &
                                     stencil_max_px,    &
                                     stencil_map_px,    &
                                     field_py,          &
                                     stencil_sizes_py,  &
                                     stencil_max_py,    &
                                     stencil_map_py,    &
                                     dep_pts,           &
                                     panel_id,          &
                                     panel_edge_dist_W, &
                                     panel_edge_dist_E, &
                                     panel_edge_dist_S, &
                                     panel_edge_dist_N, &
                                     monotone,          &
                                     ndf_wf,            &
                                     undf_wf,           &
                                     map_wf,            &
                                     ndf_w2h,           &
                                     undf_w2h,          &
                                     map_w2h,           &
                                     ndf_2d,            &
                                     undf_2d,           &
                                     map_2d )

    use horizontal_cubic_sl_kernel_mod,             only: horizontal_cubic_sl_1d
    use horizontal_cubic_sl_panel_remap_kernel_mod, only: horizontal_cubic_sl_panel_remap_1d

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
    integer(kind=i_def), intent(in) :: stencil_max_px
    integer(kind=i_def), intent(in) :: stencil_max_py
    integer(kind=i_def), intent(in) :: stencil_sizes_x(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_y(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_px(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_py(4)
    integer(kind=i_def), intent(in) :: monotone

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_wf,stencil_max_x,4)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_wf,stencil_max_y,4)
    integer(kind=i_def), intent(in) :: stencil_map_px(ndf_wf,stencil_max_px,4)
    integer(kind=i_def), intent(in) :: stencil_map_py(ndf_wf,stencil_max_py,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: increment_x(undf_wf)
    real(kind=r_tran),   intent(inout) :: increment_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_x(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_px(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_py(undf_wf)
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
    integer(kind=i_def) :: stencil_map_px_1d(-stencil_max_px:stencil_max_px)
    integer(kind=i_def) :: stencil_map_py_1d(-stencil_max_py:stencil_max_py)

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
      stencil_map_px_1d(i) = stencil_map_px(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_x_1d(i) = stencil_map_x(1, i+1, 3)
      stencil_map_px_1d(i) = stencil_map_px(1, i+1, 3)
    end do

    do i = -stencil_extent_yl, 0
      stencil_map_y_1d(i) = stencil_map_y(1, 1-i, 2)
      stencil_map_py_1d(i) = stencil_map_py(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_y_1d(i) = stencil_map_y(1, i+1, 4)
      stencil_map_py_1d(i) = stencil_map_py(1, i+1, 4)
    end do

    ! X-calculation ------------------------------------------------------------
    ! Uses field_y (already advected in y)
    near_edge = (                                                              &
        (ABS(edge_dist_W) <= MAX(stencil_extent_xl, stencil_extent_xr))        &
        .or. (ABS(edge_dist_E) <= MAX(stencil_extent_xl, stencil_extent_xr))   &
    )

    if (near_edge) then
      call horizontal_cubic_sl_panel_remap_1d( nlayers,           &
                                               .true.,            &
                                               increment_x,       &
                                               field_y,           &
                                               field_py,          &
                                               field_x,           &
                                               field_px,          &
                                               stencil_extent_xl, &
                                               stencil_extent_xr, &
                                               stencil_max_x,     &
                                               stencil_map_x_1d,  &
                                               stencil_map_py_1d, &
                                               edge_dist_W,       &
                                               edge_dist_E,       &
                                               ipanel,            &
                                               dep_pts,           &
                                               monotone,          &
                                               ndf_wf,            &
                                               undf_wf,           &
                                               map_wf,            &
                                               ndf_w2h,           &
                                               undf_w2h,          &
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
    ! Uses field_x (already advected in x)
    near_edge = (                                                              &
        (ABS(edge_dist_S) <= MAX(stencil_extent_yl, stencil_extent_yr))        &
        .or. (ABS(edge_dist_N) <= MAX(stencil_extent_yl, stencil_extent_yr))   &
    )

    if (near_edge) then
      call horizontal_cubic_sl_panel_remap_1d( nlayers,           &
                                               .false.,           &
                                               increment_y,       &
                                               field_x,           &
                                               field_px,          &
                                               field_y,           &
                                               field_py,          &
                                               stencil_extent_yl, &
                                               stencil_extent_yr, &
                                               stencil_max_y,     &
                                               stencil_map_y_1d,  &
                                               stencil_map_px_1d, &
                                               edge_dist_S,       &
                                               edge_dist_N,       &
                                               ipanel,            &
                                               dep_pts,           &
                                               monotone,          &
                                               ndf_wf,            &
                                               undf_wf,           &
                                               map_wf,            &
                                               ndf_w2h,           &
                                               undf_w2h,          &
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

  end subroutine peregrin_cubic_sl_code

end module peregrin_cubic_sl_kernel_mod
