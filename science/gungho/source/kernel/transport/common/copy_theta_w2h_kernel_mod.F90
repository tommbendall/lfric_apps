!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Builds a W2H flux field for the transport of theta variables, which
!!        differs in the top two layers of the mesh

module copy_theta_w2h_kernel_mod
use argument_mod,            only : arg_type, GH_FIELD, GH_READ,               &
                                    CELL_COLUMN, GH_REAL, GH_WRITE,            &
                                    ANY_DISCONTINUOUS_SPACE_2,                 &
                                    ANY_DISCONTINUOUS_SPACE_3, GH_INTEGER
use fs_continuity_mod,       only : W2H
use constants_mod,           only : r_tran, i_def
use kernel_mod,              only : kernel_type
use reference_element_mod,   only : N

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: copy_theta_w2h_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                          &
      arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2),    &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2),    &
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2H),                          &
      arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),    &
      arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)     &
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: copy_theta_w2h_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public copy_theta_w2h_code

contains

!> @brief Builds a W2H flux field for the transport of theta variables, which
!!        differs in the top two layers of the mesh
!> @param[in]     nlayers_shifted   Number of layers in the shifted/theta mesh
!> @param[in,out] flux_theta        The flux field to compute
!> @param[in]     flux_shifted      Flux on the shifted mesh
!> @param[in]     flux_prime        Flux on the prime extrusion
!> @param[in]     face_selector_ew  Determines faces to loop over in x direction
!> @param[in]     face_selector_ns  Determines faces to loop over in y direction
!> @param[in]     ndf_shifted       Num of DoFs per cell for shifted W2H
!> @param[in]     undf_shifted      Num of DoFs per partition for shifted W2H
!> @param[in]     map_shifted       Base cell DoF-map for shifted W2H
!> @param[in]     ndf_prime         Num of DoFs per cell for prime W2H
!> @param[in]     undf_prime        Num of DoFs per partition for prime W2H
!> @param[in]     map_prime         Base cell DoF-map for prime W2H
!> @param[in]     ndf_w3_2d         Num of DoFs per cell for W3 2D field
!> @param[in]     undf_w3_2d        Num of DoFs per partition for W3 2D field
!> @param[in]     map_w3_2d         Base cell DoF-map for W3 2D field
subroutine copy_theta_w2h_code( nlayers_shifted,  &
                                flux_theta,       &
                                flux_shifted,     &
                                flux_prime,       &
                                face_selector_ew, &
                                face_selector_ns, &
                                ndf_shifted,      &
                                undf_shifted,     &
                                map_shifted,      &
                                ndf_prime,        &
                                undf_prime,       &
                                map_prime,        &
                                ndf_w3_2d,        &
                                undf_w3_2d,       &
                                map_w3_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers_shifted
  integer(kind=i_def), intent(in)    :: undf_shifted, ndf_shifted
  integer(kind=i_def), intent(in)    :: undf_prime, ndf_prime
  integer(kind=i_def), intent(in)    :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), intent(in)    :: map_shifted(ndf_shifted)
  integer(kind=i_def), intent(in)    :: map_prime(ndf_prime)
  integer(kind=i_def), intent(in)    :: map_w3_2d(ndf_w3_2d)
  real(kind=r_tran),   intent(in)    :: flux_prime(undf_prime)
  real(kind=r_tran),   intent(in)    :: flux_shifted(undf_shifted)
  real(kind=r_tran),   intent(inout) :: flux_theta(undf_shifted)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)

  ! Internal variables
  integer(kind=i_def) :: idx_sh, idx_prime, nlayers_prime
  integer(kind=i_def) :: face, j

  nlayers_prime = nlayers_shifted - 1

  ! Loop over horizontal W2 DoFs
  do j = 1, face_selector_ew(map_w3_2d(1)) + face_selector_ns(map_w3_2d(1))
    face = j
    if (j == 3 .and. face_selector_ns(map_w3_2d(1)) == 2                       &
               .and. face_selector_ew(map_w3_2d(1)) == 1) face = N

    idx_sh = map_shifted(face)
    idx_prime = map_prime(face)
    nlayers_prime = nlayers_shifted - 1

    ! For all but top two layers, copy the shifted flux
    if (nlayers_shifted > 2) then
      flux_theta(idx_sh : idx_sh+nlayers_shifted-3) =                          &
          flux_shifted(idx_sh : idx_sh+nlayers_shifted-3)
    end if

    ! For the next-to-top value, combine the top two shifted flux values
    flux_theta(idx_sh+nlayers_shifted-2) =                                     &
        flux_shifted(idx_sh+nlayers_shifted-2) + flux_shifted(idx_sh+nlayers_shifted-1)

    ! For very top value, copy the *prime* flux
    flux_theta(idx_sh+nlayers_shifted-1) = flux_prime(idx_prime+nlayers_prime-1)
  end do

end subroutine copy_theta_w2h_code

end module copy_theta_w2h_kernel_mod