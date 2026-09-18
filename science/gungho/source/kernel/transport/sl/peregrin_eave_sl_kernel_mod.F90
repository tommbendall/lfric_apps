!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates fields in x and y at time n+1 near cubed-sphere panel
!!          edges, using 1D semi-Lagrangian transport of the PEREGRIN wind.
!> @details This is the PEREGRIN-specific counterpart of
!!          `horizontal_linear_sl_kernel_mod.F90`: it solves the same
!!          one-dimensional advection equation in x and y using linear
!!          interpolation, but uses the pre-computed PEREGRIN departure
!!          distances (`dep_dist_peregrin`) rather than the raw departure
!!          points.
!!
!!          Unlike the general SL step, this kernel does not transport the
!!          plain advected fields `field_x`/`field_y`. Instead:
!!            - the x-direction correction (near an S/N panel edge, i.e. the
!!              "y-eave") transports `field_ry` (field_n remapped in y);
!!            - the y-direction correction (near a W/E panel edge, i.e. the
!!              "x-eave") transports `field_rx` (field_n remapped in x).
!!          If the 1D stencil for one of these corrections itself crosses a
!!          panel edge, that means the column is in fact at a panel corner,
!!          so the corner-remapped field (`corner_field_x`/`corner_field_y`)
!!          is used instead -- exactly analogous to the way the general,
!!          non-PEREGRIN transport picks up an edge-remapped field when its
!!          stencil crosses an edge.
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest
!!       order.
!> @note Transport in x is only performed near a panel edge in the
!!       y-direction (S/N); transport in y is only performed near a panel
!!       edge in the x-direction (W/E). Elsewhere the field is passed
!!       through unchanged (using field_x/field_y).

module peregrin_eave_sl_kernel_mod

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
  use log_mod,               only: log_event, log_scratch_space, LOG_LEVEL_DEBUG

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! PSy layer
  type, public, extends(kernel_type) :: peregrin_eave_sl_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                      &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field_px
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field_py
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_rx (field_n remapped in x)
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! field_ry (field_n remapped in y)
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! corner_field_x (field_n corner-remapped in x)
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                            STENCIL(CROSS2D)), & ! corner_field_y (field_n corner-remapped in y)
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H),                       & ! dep_dist_peregrin_x
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H),                       & ! dep_dist_peregrin_y
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_id
        arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  & ! panel_edge_dist
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: peregrin_eave_sl_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: peregrin_eave_sl_code

contains

  !> @brief Compute PEREGRIN advective transport in x and y directions using
  !!        a 1D semi-Lagrangian scheme with a linear reconstruction, near
  !!        cubed-sphere panel edges.
  !> @param[in]     nlayers               Number of layers
  !> @param[in,out] field_px              Field at time n+1 in x direction
  !> @param[in,out] field_py              Field at time n+1 in y direction
  !> @param[in]     field_rx              field_n remapped in x, transported
  !!                                      in y near a W/E edge (the "x-eave")
  !> @param[in]     field_ry              field_n remapped in y, transported
  !!                                      in x near an S/N edge (the "y-eave")
  !> @param[in]     corner_field_x        field_ry corner-remapped in x, used
  !!                                      when the x-direction transport of
  !!                                      field_ry itself crosses a W/E edge
  !!                                      (i.e. this column is at a corner)
  !> @param[in]     corner_field_y        field_rx corner-remapped in y, used
  !!                                      when the y-direction transport of
  !!                                      field_rx itself crosses an S/N edge
  !!                                      (i.e. this column is at a corner)
  !> @param[in]     dep_dist_peregrin_x   PEREGRIN horizontal departure
  !!                                      distances at cell faces, calculated
  !!                                      from remapped_flux_x. The corrected
  !!                                      DoFs of remapped_flux_x are its S/N
  !!                                      faces, which are read by the
  !!                                      y-transport, so this is used for the
  !!                                      y-transport calculation
  !> @param[in]     dep_dist_peregrin_y   PEREGRIN horizontal departure
  !!                                      distances at cell faces, calculated
  !!                                      from remapped_flux_y. The corrected
  !!                                      DoFs of remapped_flux_y are its W/E
  !!                                      faces, which are read by the
  !!                                      x-transport, so this is used for the
  !!                                      x-transport calculation
  !> @param[in]     panel_edge_dist_W     2D field containing the distance of
  !!                                      each column from the panel edge to
  !!                                      the West
  !> @param[in]     panel_edge_dist_E     2D field containing the distance of
  !!                                      each column from the panel edge to
  !!                                      the East
  !> @param[in]     panel_edge_dist_S     2D field containing the distance of
  !!                                      each column from the panel edge to
  !!                                      the South
  !> @param[in]     panel_edge_dist_N     2D field containing the distance of
  !!                                      each column from the panel edge to
  !!                                      the North
  !> @param[in]     ndf_wf                Num of DoFs for field per cell
  !> @param[in]     undf_wf               Num of DoFs in this partition for
  !!                                      field
  !> @param[in]     map_wf                Map for field
  !> @param[in]     ndf_w2h               Num of DoFs for W2H per cell
  !> @param[in]     undf_w2h              Num of DoFs in this partition for W2H
  !> @param[in]     map_w2h               Map for W2H
  !> @param[in]     ndf_2d                Num of DoFs for 2D fields per cell
  !> @param[in]     undf_2d               Num of DoFs in this partition for 2D
  !!                                      fields
  !> @param[in]     map_2d                Map for 2D fields
  subroutine peregrin_eave_sl_code( nlayers,                &
                                    field_px,               &
                                    field_py,               &
                                    field_rx,               &
                                    stencil_sizes_rx,       &
                                    stencil_max_rx,         &
                                    stencil_map_rx,         &
                                    field_ry,               &
                                    stencil_sizes_ry,       &
                                    stencil_max_ry,         &
                                    stencil_map_ry,         &
                                    corner_field_x,         &
                                    stencil_sizes_cx,       &
                                    stencil_max_cx,         &
                                    stencil_map_cx,         &
                                    corner_field_y,         &
                                    stencil_sizes_cy,       &
                                    stencil_max_cy,         &
                                    stencil_map_cy,         &
                                    dep_dist_peregrin_x,    &
                                    dep_dist_peregrin_y,    &
                                    panel_id,               &
                                    panel_edge_dist_W,      &
                                    panel_edge_dist_E,      &
                                    panel_edge_dist_S,      &
                                    panel_edge_dist_N,      &
                                    ndf_wf,                 &
                                    undf_wf,                &
                                    map_wf,                 &
                                    ndf_w2h,                &
                                    undf_w2h,               &
                                    map_w2h,                &
                                    ndf_2d,                 &
                                    undf_2d,                &
                                    map_2d )

    use horizontal_linear_sl_panel_remap_kernel_mod,                           &
        only: horizontal_linear_sl_panel_remap_1d
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
    integer(kind=i_def), intent(in) :: stencil_max_rx
    integer(kind=i_def), intent(in) :: stencil_max_ry
    integer(kind=i_def), intent(in) :: stencil_max_cx
    integer(kind=i_def), intent(in) :: stencil_max_cy
    integer(kind=i_def), intent(in) :: stencil_sizes_rx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_ry(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_cx(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_cy(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
    integer(kind=i_def), intent(in) :: stencil_map_rx(ndf_wf,stencil_max_rx,4)
    integer(kind=i_def), intent(in) :: stencil_map_ry(ndf_wf,stencil_max_ry,4)
    integer(kind=i_def), intent(in) :: stencil_map_cx(ndf_wf,stencil_max_cx,4)
    integer(kind=i_def), intent(in) :: stencil_map_cy(ndf_wf,stencil_max_cy,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: field_px(undf_wf)
    real(kind=r_tran),   intent(inout) :: field_py(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_rx(undf_wf)
    real(kind=r_tran),   intent(in)    :: field_ry(undf_wf)
    real(kind=r_tran),   intent(in)    :: corner_field_x(undf_wf)
    real(kind=r_tran),   intent(in)    :: corner_field_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_dist_peregrin_x(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_dist_peregrin_y(undf_w2h)
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
    integer(kind=i_def) :: stencil_map_ry_1d(-stencil_max_ry:stencil_max_ry)
    integer(kind=i_def) :: stencil_map_cx_1d(-stencil_max_cx:stencil_max_cx)
    integer(kind=i_def) :: stencil_map_rx_1d(-stencil_max_rx:stencil_max_rx)
    integer(kind=i_def) :: stencil_map_cy_1d(-stencil_max_cy:stencil_max_cy)
    logical(kind=l_def) :: do_peregrin, near_edge
    integer(kind=i_def) :: k, nl

    ipanel = INT(panel_id(map_2d(1)), i_def)
    edge_dist_E = panel_edge_dist_E(map_2d(1))
    edge_dist_W = panel_edge_dist_W(map_2d(1))
    edge_dist_S = panel_edge_dist_S(map_2d(1))
    edge_dist_N = panel_edge_dist_N(map_2d(1))

    ! nl = nlayers   for W3, nlayers+1 for Wtheta
    nl = nlayers + ndf_wf - 1

    ! Form 1D stencils for the x-direction transport of field_ry, and the
    ! corner-remapped field used when that transport crosses a W/E edge
    stencil_extent_xl = stencil_sizes_ry(1) - 1
    stencil_extent_xr = stencil_sizes_ry(3) - 1

    do i = -stencil_extent_xl, 0
      stencil_map_ry_1d(i) = stencil_map_ry(1, 1-i, 1)
      stencil_map_cx_1d(i) = stencil_map_cx(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_ry_1d(i) = stencil_map_ry(1, i+1, 3)
      stencil_map_cx_1d(i) = stencil_map_cx(1, i+1, 3)
    end do

    ! Form 1D stencils for the y-direction transport of field_rx, and the
    ! corner-remapped field used when that transport crosses an S/N edge
    stencil_extent_yl = stencil_sizes_rx(2) - 1
    stencil_extent_yr = stencil_sizes_rx(4) - 1

    do i = -stencil_extent_yl, 0
      stencil_map_rx_1d(i) = stencil_map_rx(1, 1-i, 2)
      stencil_map_cy_1d(i) = stencil_map_cy(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_rx_1d(i) = stencil_map_rx(1, i+1, 4)
      stencil_map_cy_1d(i) = stencil_map_cy(1, i+1, 4)
    end do

    ! X-calculation --------------------------------------------------------
    ! Only carry out the PEREGRIN transport in x when this column is near a
    ! panel edge in the y-direction (S/N); otherwise just copy field_x
    ! through unchanged. This transports field_ry (field_n remapped in y),
    ! switching to corner_field_x if the x-direction stencil also crosses a
    ! W/E edge (i.e. this column is at a panel corner).
    do_peregrin = (                                                            &
        (ABS(edge_dist_S) <= MAX(stencil_extent_yl, stencil_extent_yr))        &
        .or. (ABS(edge_dist_N) <= MAX(stencil_extent_yl, stencil_extent_yr))   &
    )

    if (do_peregrin) then
      near_edge = (                                                            &
          (ABS(edge_dist_W) <= MAX(stencil_extent_xl, stencil_extent_xr))      &
          .or. (ABS(edge_dist_E) <= MAX(stencil_extent_xl, stencil_extent_xr)) &
      )

      if (near_edge) then
        call horizontal_linear_sl_panel_remap_1d( nlayers,                &
                                                  .true.,                 &
                                                  field_px,               &
                                                  field_ry,               &
                                                  corner_field_x,         &
                                                  corner_field_y,         &
                                                  stencil_extent_xl,      &
                                                  stencil_extent_xr,      &
                                                  stencil_max_ry,         &
                                                  stencil_map_ry_1d,      &
                                                  stencil_map_cx_1d,      &
                                                  edge_dist_W,            &
                                                  edge_dist_E,            &
                                                  ipanel,                 &
                                      ! The x-transport reads W/E DoFs, which
                                      ! are the corrected ones in _peregrin_y
                                                  dep_dist_peregrin_y,    &
                                                  ndf_wf,                 &
                                                  undf_wf,                &
                                                  map_wf,                 &
                                                  ndf_w2h,                &
                                                  undf_w2h,               &
                                                  map_w2h )
      else
        call horizontal_linear_sl_1d( nlayers,             &
                                      .true.,              &
                                      field_px,            &
                                      field_ry,            &
                                      stencil_extent_xl,   &
                                      stencil_extent_xr,   &
                                      stencil_max_ry,      &
                                      stencil_map_ry_1d,   &
                                      ! The x-transport reads W/E DoFs, which
                                      ! are the corrected ones in _peregrin_y
                                      dep_dist_peregrin_y, &
                                      ndf_wf,              &
                                      undf_wf,             &
                                      map_wf,              &
                                      ndf_w2h,             &
                                      undf_w2h,            &
                                      map_w2h )
      end if

      ! TODO: to remove
      select case (map_wf(1))
      case (862, 863, 864)
        write(log_scratch_space,*)                                              &
            'PEREGRIN_EAVE_X_RESULT: cell=', map_wf(1),                         &
            ' field_px(0)=', field_px(map_wf(1))
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
      end select
    else
      do k = 0, nl-1
        field_px(map_wf(1)+k) = -2900.0_r_tran
      end do
    end if

    ! Y-calculation --------------------------------------------------------
    ! Only carry out the PEREGRIN transport in y when this column is near a
    ! panel edge in the x-direction (W/E); otherwise just copy field_y
    ! through unchanged. This transports field_rx (field_n remapped in x),
    ! switching to corner_field_y if the y-direction stencil also crosses an
    ! S/N edge (i.e. this column is at a panel corner).
    do_peregrin = (                                                            &
        (ABS(edge_dist_W) <= MAX(stencil_extent_xl, stencil_extent_xr))        &
        .or. (ABS(edge_dist_E) <= MAX(stencil_extent_xl, stencil_extent_xr))   &
    )

    if (do_peregrin) then
      near_edge = (                                                            &
          (ABS(edge_dist_S) <= MAX(stencil_extent_yl, stencil_extent_yr))      &
          .or. (ABS(edge_dist_N) <= MAX(stencil_extent_yl, stencil_extent_yr)) &
      )

      if (near_edge) then
        call horizontal_linear_sl_panel_remap_1d( nlayers,                &
                                                  .false.,                &
                                                  field_py,               &
                                                  field_rx,               &
                                                  corner_field_y,         &
                                                  corner_field_x,         &
                                                  stencil_extent_yl,      &
                                                  stencil_extent_yr,      &
                                                  stencil_max_rx,         &
                                                  stencil_map_rx_1d,      &
                                                  stencil_map_cy_1d,      &
                                                  edge_dist_S,            &
                                                  edge_dist_N,            &
                                                  ipanel,                 &
                                      ! The y-transport reads S/N DoFs, which
                                      ! are the corrected ones in _peregrin_x
                                                  dep_dist_peregrin_x,    &
                                                  ndf_wf,                 &
                                                  undf_wf,                &
                                                  map_wf,                 &
                                                  ndf_w2h,                &
                                                  undf_w2h,               &
                                                  map_w2h )
      else
        call horizontal_linear_sl_1d( nlayers,             &
                                      .false.,             &
                                      field_py,            &
                                      field_rx,            &
                                      stencil_extent_yl,   &
                                      stencil_extent_yr,   &
                                      stencil_max_rx,      &
                                      stencil_map_rx_1d,   &
                                      ! The y-transport reads S/N DoFs, which
                                      ! are the corrected ones in _peregrin_x
                                      dep_dist_peregrin_x, &
                                      ndf_wf,              &
                                      undf_wf,             &
                                      map_wf,              &
                                      ndf_w2h,             &
                                      undf_w2h,            &
                                      map_w2h )
      end if

    else
      do k = 0, nl-1
        field_py(map_wf(1)+k) = -2900.0_r_tran
      end do
    end if

  end subroutine peregrin_eave_sl_code

end module peregrin_eave_sl_kernel_mod
