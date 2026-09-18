!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Computes the coordinates corresponding to extending a cubed-sphere
!!        panel "around" a corner, using the two edge-neighbouring panels that
!!        meet at that corner.
!! @details At each cubed-sphere corner three panels meet: the owned panel and
!!          its x-crossing (W/E) and y-crossing (S/N) neighbours. The
!!          equiangular gnomonic mapping for a single panel, xi=tan(alpha),
!!          eta=tan(beta), is a smooth, single-valued, invertible function of
!!          (alpha,beta) over the whole range (-pi/2,pi/2), not merely that
!!          panel's own quarter (-pi/4,pi/4). So the corner-extended
!!          coordinate for a neighbouring panel is simply the direct
!!          reprojection of this column's true physical location into that
!!          neighbour's own frame - no intermediate "hop" through a third
!!          panel is needed. This kernel computes two such reprojections:
!!            - "anticlockwise": reprojected into the x-crossing neighbour
!!            - "clockwise": reprojected into the y-crossing neighbour
!!          These are genuinely different values, since they are each
!!          expressed in a different panel's own equiangular convention.
module panel_corner_coords_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER,                          &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_7,                    &
                                 OWNED_AND_HALO_CELL_COLUMN
use constants_mod,         only: r_def, i_def, l_def


use log_mod,               only: log_event, log_level_info, log_scratch_space

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_corner_coords_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                          &
      arg_type(GH_FIELD*2, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_7),   &
      arg_type(GH_FIELD*2, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_7),   &
      arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7),   &
      arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7),   &
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   &
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   &
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                                &
  /)
  integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
contains
  procedure, nopass :: panel_corner_coords_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_corner_coords_code

contains

!> @brief Extend the equiangular coordinate fields around cubed-sphere corners
!> @param[in]     nlayers             Number of layers
!> @param[in,out] alpha_anti          Anticlockwise alpha extension
!> @param[in,out] beta_anti           Anticlockwise beta extension
!> @param[in,out] alpha_clock         Clockwise alpha extension
!> @param[in,out] beta_clock          Clockwise beta extension
!> @param[in]     alpha_x             alpha in near panel edge in x direction
!> @param[in]     beta_x              beta in near panel edge in x direction
!> @param[in]     alpha_y             alpha in near panel edge in y direction
!> @param[in]     beta_y              beta in near panel edge in y direction
!> @param[in]     panel_id            ID of the panel for each column
!> @param[in]     panel_edge_dist_W   Distance to the West panel edge
!> @param[in]     panel_edge_dist_E   Distance to the East panel edge
!> @param[in]     panel_edge_dist_S   Distance to the South panel edge
!> @param[in]     panel_edge_dist_N   Distance to the North panel edge
!> @param[in]     stencil_extent      Max stencil extent to be used
!> @param[in]     ndf_wx_2d           Num DoFs per cell for 2D coords
!> @param[in]     undf_wx_2d          Num DoFs for this partition
!> @param[in]     map_wx_2d           DoFmap for 2D coord fields
!> @param[in]     ndf_pid             Num DoFs per cell for panel ID
!> @param[in]     undf_pid            Num DoFs for this partition
!> @param[in]     map_pid             DoFmap for panel ID
subroutine panel_corner_coords_code( nlayers,                                  &
                                     alpha_anti, beta_anti,                    &
                                     alpha_clock, beta_clock,                  &
                                     alpha_x, beta_x, alpha_y, beta_y,         &
                                     panel_id,                                 &
                                     panel_edge_dist_W, panel_edge_dist_E,     &
                                     panel_edge_dist_S, panel_edge_dist_N,     &
                                     stencil_extent,                           &
                                     ndf_wx_2d, undf_wx_2d, map_wx_2d,         &
                                     ndf_pid, undf_pid, map_pid                &
                                   )

  use panel_edge_support_mod,       only: panel_neighbour
  use reference_element_mod,        only: W, S, N, E

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_pid, ndf_wx_2d
  integer(kind=i_def), intent(in)    :: undf_pid, undf_wx_2d
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: stencil_extent
  real(kind=r_def),    intent(inout) :: alpha_anti(undf_wx_2d)
  real(kind=r_def),    intent(inout) :: beta_anti(undf_wx_2d)
  real(kind=r_def),    intent(inout) :: alpha_clock(undf_wx_2d)
  real(kind=r_def),    intent(inout) :: beta_clock(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: alpha_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: alpha_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)

  integer(kind=i_def) :: owned_panel, swapped_panel_x, swapped_panel_y
  integer(kind=i_def) :: panel_W, panel_E, panel_S, panel_N, panel_corner, df

  ! Output
  real(kind=r_def) :: alpha_anti_c, beta_anti_c, alpha_clock_c, beta_clock_c

  ! Panel id for this column
  owned_panel = int(panel_id(map_pid(1)), i_def)

  panel_W = panel_neighbour(owned_panel, W)
  panel_E = panel_neighbour(owned_panel, E)
  panel_S = panel_neighbour(owned_panel, S)
  panel_N = panel_neighbour(owned_panel, N)

  ! Determine if we are near a corner of the panel ----------------------------
  ! Initialise swapped panels
  swapped_panel_x = 0
  swapped_panel_y = 0

  ! At the moment we can only set one neighbour for the x/y directions, so use
  ! whichever edge is closest (mirrors panel_edge_coords_code)
  if (ABS(panel_edge_dist_W(map_pid(1))) <= stencil_extent) then
    swapped_panel_x = panel_W
  end if
  if (ABS(panel_edge_dist_E(map_pid(1))) <= stencil_extent .and.               &
      ABS(panel_edge_dist_E(map_pid(1))) < ABS(panel_edge_dist_W(map_pid(1)))) then
    swapped_panel_x = panel_E
  end if

  if (ABS(panel_edge_dist_S(map_pid(1))) <= stencil_extent) then
    swapped_panel_y = panel_S
  end if
  if (ABS(panel_edge_dist_N(map_pid(1))) <= stencil_extent .and.               &
      ABS(panel_edge_dist_N(map_pid(1))) < ABS(panel_edge_dist_S(map_pid(1)))) then
    swapped_panel_y = panel_N
  end if

  if ( swapped_panel_x /= 0_i_def .and. swapped_panel_y /= 0_i_def ) then
    ! Determine which coordinate is the extended coordinate
    panel_corner = 100*owned_panel + 10*swapped_panel_x + swapped_panel_y
    select case (panel_corner)
    ! ------------------------------------------------------------------------
    ! W/E-neighbour edge unrotated, S/N-neighbour edge unrotated
    case (125, 216, 364, 453, 541, 632)
      do df = 1, ndf_wx_2d
        alpha_anti(map_wx_2d(df)) = -beta_y(map_wx_2d(df))
        beta_anti(map_wx_2d(df)) = alpha_x(map_wx_2d(df))
        alpha_clock(map_wx_2d(df)) = beta_y(map_wx_2d(df))
        beta_clock(map_wx_2d(df)) = -alpha_x(map_wx_2d(df))
      end do

    ! ------------------------------------------------------------------------
    ! W/E-neighbour edge unrotated, S/N-neighbour edge rotated
    case (126, 215, 362, 451, 543, 634)
      do df = 1, ndf_wx_2d
        alpha_anti(map_wx_2d(df)) = alpha_x(map_wx_2d(df))
        beta_anti(map_wx_2d(df)) = alpha_y(map_wx_2d(df))
        alpha_clock(map_wx_2d(df)) = alpha_y(map_wx_2d(df))
        beta_clock(map_wx_2d(df)) = -alpha_x(map_wx_2d(df))
      end do

    ! ------------------------------------------------------------------------
    ! W/E-neighbour edge rotated, S/N-neighbour edge unrotated
    case (145, 236, 354, 463, 521, 612)
      do df = 1, ndf_wx_2d
        alpha_anti(map_wx_2d(df)) = -beta_x(map_wx_2d(df))
        beta_anti(map_wx_2d(df)) = beta_y(map_wx_2d(df))
        alpha_clock(map_wx_2d(df)) = beta_y(map_wx_2d(df))
        beta_clock(map_wx_2d(df)) = beta_x(map_wx_2d(df))
      end do

    ! ------------------------------------------------------------------------
    ! W/E-neighbour edge rotated, S/N-neighbour edge rotated
    case (146, 235, 352, 461, 523, 614)
      do df = 1, ndf_wx_2d
        alpha_anti(map_wx_2d(df)) = alpha_y(map_wx_2d(df))
        beta_anti(map_wx_2d(df)) = beta_x(map_wx_2d(df))
        alpha_clock(map_wx_2d(df)) = alpha_y(map_wx_2d(df))
        beta_clock(map_wx_2d(df)) = beta_x(map_wx_2d(df))
      end do
    end select

  else
    ! Not near a corner: fall back to the native coordinates
    do df = 1, ndf_wx_2d
      alpha_anti(map_wx_2d(df)) = alpha_x(map_wx_2d(df))
      beta_anti(map_wx_2d(df)) = beta_x(map_wx_2d(df))
      alpha_clock(map_wx_2d(df)) = alpha_y(map_wx_2d(df))
      beta_clock(map_wx_2d(df)) = beta_y(map_wx_2d(df))
    end do
  end if

  if ( swapped_panel_x /= 0_i_def .and. swapped_panel_y /= 0_i_def ) then
    alpha_anti_c = 0.0_r_def
    beta_anti_c = 0.0_r_def
    alpha_clock_c = 0.0_r_def
    beta_clock_c = 0.0_r_def
    do df = 1, ndf_wx_2d
      alpha_anti_c = alpha_anti_c +                                      &
          alpha_anti(map_wx_2d(df)) / real(ndf_wx_2d, r_def)
      beta_anti_c = beta_anti_c +                                         &
          beta_anti(map_wx_2d(df)) / real(ndf_wx_2d, r_def)
      alpha_clock_c = alpha_clock_c +                                       &
          alpha_clock(map_wx_2d(df)) / real(ndf_wx_2d, r_def)
      beta_clock_c = beta_clock_c +                                         &
          beta_clock(map_wx_2d(df)) / real(ndf_wx_2d, r_def)
    end do
  end if

end subroutine panel_corner_coords_code

end module panel_corner_coords_kernel_mod
