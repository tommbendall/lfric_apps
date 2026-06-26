!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Process Jules sea and sea-ice ancillaries to apply an SST anomaly
!> @details Kernel used without proper PSyclone support of multi-data fields
!>          see https://github.com/stfc/PSyclone/issues/868
module process_sstanom_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_REAL,         &
                           GH_READ, GH_READWRITE,     &
                           CELL_COLUMN,               &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  use jules_control_init_mod, only: n_sea_ice_tile, &
       first_sea_tile, first_sea_ice_tile
  use surface_config_mod,     only : sst_anomaly, sst_anomaly_single, &
                                     sst_anomaly_single_value

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: process_sstanom_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: process_sstanom_code
  end type process_sstanom_kernel_type

  public :: process_sstanom_code

contains

  !> @param[in]     nlayers            The number of layers
  !> @param[in]     tile_fraction      Surface tile fractions
  !> @param[in,out] tile_temperature   Tile temperature
  !> @param[in]     ndf_tile           Number of DOFs per cell for tiles
  !> @param[in]     undf_tile          Number of total DOFs for tiles
  !> @param[in]     map_tile           Dofmap for cell for surface tiles
  subroutine process_sstanom_code(nlayers,                       &
                                  tile_fraction,                 &
                                  tile_temperature,              &
                                  ndf_tile, undf_tile, map_tile)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)

    real(kind=r_def), intent(in)     :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout)  :: tile_temperature(undf_tile)

    ! Internal variables
    integer(kind=i_def) :: i
    real(kind=r_def) :: tot_ice, sstanom

    ! Calculate the current ice fraction for use below
    tot_ice = 0.0_r_def
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      tot_ice = tot_ice + tile_fraction(map_tile(1)+i-1)
    end do

    ! Apply the SST anomaly
    if (sst_anomaly == sst_anomaly_single) then
      if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def .and. &
           tot_ice == 0.0_r_def) then
          tile_temperature(map_tile(1)+first_sea_tile-1) = &
            tile_temperature(map_tile(1)+first_sea_tile-1) + sst_anomaly_single_value
      end if
    end if

  end subroutine process_sstanom_code

end module process_sstanom_kernel_mod
