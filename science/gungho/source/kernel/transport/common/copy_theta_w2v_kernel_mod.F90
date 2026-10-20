!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Builds a W2V flux field for the transport of theta variables, which
!!        has no vertical transport into the top layer of the mesh
module copy_theta_w2v_kernel_mod
use argument_mod,            only : arg_type, GH_FIELD, GH_READ,               &
                                    CELL_COLUMN, GH_REAL, GH_WRITE,            &
                                    ANY_DISCONTINUOUS_SPACE_3
use fs_continuity_mod,       only : W2V
use constants_mod,           only : r_tran, i_def
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: copy_theta_w2v_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                          &
      arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2V),                          &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2V)                           &
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: copy_theta_w2v_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public copy_theta_w2v_code

contains

!> @brief Builds a W2V flux field for the transport of theta variables, which
!!        has no vertical transport into the top layer of the mesh
!> @param[in]     nlayers_shifted   Number of layers in the shifted/theta mesh
!> @param[in,out] flux_theta        The flux field to compute
!> @param[in]     flux_shifted      Flux on the shifted mesh
!> @param[in]     ndf_shifted       Num of DoFs per cell for shifted W2H
!> @param[in]     undf_shifted      Num of DoFs per partition for shifted W2H
!> @param[in]     map_shifted       Base cell DoF-map for shifted W2H
subroutine copy_theta_w2v_code( nlayers_shifted,  &
                                flux_theta,       &
                                flux_shifted,     &
                                ndf_shifted,      &
                                undf_shifted,     &
                                map_shifted )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers_shifted
  integer(kind=i_def), intent(in)    :: undf_shifted, ndf_shifted
  integer(kind=i_def), intent(in)    :: map_shifted(ndf_shifted)
  real(kind=r_tran),   intent(in)    :: flux_shifted(undf_shifted)
  real(kind=r_tran),   intent(inout) :: flux_theta(undf_shifted)

  ! Internal variables
  integer(kind=i_def) :: idx_sh

  idx_sh = map_shifted(1)

  ! For all but bottom two layers, copy the shifted flux
  flux_theta(idx_sh : idx_sh+nlayers_shifted-2) =                              &
      flux_shifted(idx_sh : idx_sh+nlayers_shifted-2)

  ! Top two values are set to zero
  flux_theta(idx_sh+nlayers_shifted-1) = 0.0_r_tran
  flux_theta(idx_sh+nlayers_shifted)   = 0.0_r_tran

end subroutine copy_theta_w2v_code

end module copy_theta_w2v_kernel_mod