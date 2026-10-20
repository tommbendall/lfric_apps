!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Builds a mass field for the transport of theta variables, which
!!        differs in the top two layers of the mesh

module copy_theta_mass_kernel_mod
use argument_mod,            only : arg_type, GH_FIELD, GH_READ,               &
                                    CELL_COLUMN, GH_REAL, GH_WRITE,            &
                                    ANY_DISCONTINUOUS_SPACE_3
use fs_continuity_mod,       only : W3
use constants_mod,           only : r_tran, i_def
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: copy_theta_mass_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                          &
      arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_3),    &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3),    &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3)                            &
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: copy_theta_mass_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public copy_theta_mass_code

contains

!> @brief Builds a mass field for the transport of theta variables, which
!!        differs in the top two layers of the mesh
!> @param[in]     nlayers_shifted  Number of layers in the shifted/theta mesh
!> @param[in,out] mass_theta       The mass field to compute
!> @param[in]     mass_shifted     Mass on the shifted mesh
!> @param[in]     mass_prime       Mass on the prime extrusion
!> @param[in]     ndf_shifted      Num of DoFs per cell for shifted W3
!> @param[in]     undf_shifted     Num of DoFs per partition for shifted W3
!> @param[in]     map_shifted      Base cell DoF-map for shifted W3
!> @param[in]     ndf_prime        Num of DoFs per cell for prime W3
!> @param[in]     undf_prime       Num of DoFs per partition for prime W3
!> @param[in]     map_prime        Base cell DoF-map for prime W3
subroutine copy_theta_mass_code( nlayers_shifted, &
                                 mass_theta,      &
                                 mass_shifted,    &
                                 mass_prime,      &
                                 ndf_shifted,     &
                                 undf_shifted,    &
                                 map_shifted,     &
                                 ndf_prime,       &
                                 undf_prime,      &
                                 map_prime )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers_shifted
  integer(kind=i_def), intent(in)    :: undf_shifted, ndf_shifted
  integer(kind=i_def), intent(in)    :: undf_prime, ndf_prime
  integer(kind=i_def), intent(in)    :: map_shifted(ndf_shifted)
  integer(kind=i_def), intent(in)    :: map_prime(ndf_prime)
  real(kind=r_tran),   intent(in)    :: mass_prime(undf_prime)
  real(kind=r_tran),   intent(in)    :: mass_shifted(undf_shifted)
  real(kind=r_tran),   intent(inout) :: mass_theta(undf_shifted)

  ! Internal variables
  integer(kind=i_def) :: idx_sh, idx_prime, nlayers_prime

  idx_sh = map_shifted(1)
  idx_prime = map_prime(1)
  nlayers_prime = nlayers_shifted - 1

  ! For all but top two layers, copy the shifted mass
  if (nlayers_shifted > 2) then
    mass_theta(idx_sh : idx_sh+nlayers_shifted-3) =                            &
        mass_shifted(idx_sh : idx_sh+nlayers_shifted-3)
  end if

  ! For the next-to-top value, combine the top two shifted mass values
  mass_theta(idx_sh+nlayers_shifted-2) =                                       &
      mass_shifted(idx_sh+nlayers_shifted-2) + mass_shifted(idx_sh+nlayers_shifted-1)

  ! For very top value, copy the *prime* mass
  mass_theta(idx_sh+nlayers_shifted-1) = mass_prime(idx_prime+nlayers_prime-1)

end subroutine copy_theta_mass_code

end module copy_theta_mass_kernel_mod