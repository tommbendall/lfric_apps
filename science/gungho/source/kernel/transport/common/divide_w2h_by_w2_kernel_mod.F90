!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! GitHub Copilot.
!-------------------------------------------------------------------------------
!> @brief Divide a W2H field by the horizontal DoFs of a co-located W2 field.
!!
!> @details W2's first 4 (horizontal) DoFs per cell are ordered identically to
!!          W2H's 4 DoFs (W, S, N, E), so the W2H DoF loop counter can be used
!!          directly to index into the W2 field's DoF map (see the same trick
!!          used by panel_edge_remap_w2h_weights_kernel_mod for dA_at_w2).
!!          This kernel exists so that area-normalisation fields such as
!!          dA_at_w2_native (a W2 field) can be applied to a W2H flux field
!!          without needing to pass a W2H-space copy of that field into other
!!          kernels, reducing their memory accesses.
!!          W2H faces are shared (continuous) between neighbouring columns at
!!          lowest order, so face_selector_ew/ns (as used by
!!          split_w2_field_kernel_mod) are used to visit each shared face only
!!          once, avoiding duplicate divide work.
module divide_w2h_by_w2_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD,                                     &
                                 GH_REAL, GH_INTEGER,                          &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 CELL_COLUMN
use fs_continuity_mod,      only: W2h, W2
use constants_mod,          only: r_def, r_tran, i_def, l_def
use sci_face_selector_support_mod, only: face_from_face_selector

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: divide_w2h_by_w2_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                          &
       arg_type(GH_FIELD, GH_REAL,    GH_WRITE, W2h),                         & ! field_out
       arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2h),                         & ! field_in
       arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                          & ! w2_divisor
       arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   & ! face_selector_ew
       arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)    & ! face_selector_ns
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: divide_w2h_by_w2_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: divide_w2h_by_w2_code

contains

!> @brief Divide a W2H field by the horizontal DoFs of a W2 field.
!> @param[in]     nlayers          Number of layers
!> @param[in,out] field_out        Result of field_in / w2_divisor (W2H)
!> @param[in]     field_in         Field to be divided (W2H)
!> @param[in]     w2_divisor       Divisor field (W2)
!> @param[in]     face_selector_ew 2D field indicating which W/E faces
!!                                 to loop over in this column
!> @param[in]     face_selector_ns 2D field indicating which N/S faces
!!                                 to loop over in this column
!> @param[in]     ndf_w2h    Num DoFs per cell for the W2H fields
!> @param[in]     undf_w2h   Num DoFs for this partition for the W2H fields
!> @param[in]     map_w2h    DoF map for the W2H fields
!> @param[in]     ndf_w2     Num DoFs per cell for the W2 field
!> @param[in]     undf_w2    Num DoFs for this partition for the W2 field
!> @param[in]     map_w2     DoF map for the W2 field
!> @param[in]     ndf_fs2d   Num DoFs per cell for the 2D face selector fields
!> @param[in]     undf_fs2d  Num DoFs for this partition for the face selector
!!                          fields
!> @param[in]     map_fs2d   DoF map for the face selector fields
subroutine divide_w2h_by_w2_code( nlayers,                                     &
                                  field_out,                                  &
                                  field_in,                                   &
                                  w2_divisor,                                 &
                                  face_selector_ew,                           &
                                  face_selector_ns,                          &
                                  ndf_w2h, undf_w2h, map_w2h,                 &
                                  ndf_w2, undf_w2, map_w2,                    &
                                  ndf_fs2d, undf_fs2d, map_fs2d               &
  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in)    :: ndf_fs2d, undf_fs2d

  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_fs2d(ndf_fs2d)

  real(kind=r_tran),   intent(inout) :: field_out(undf_w2h)
  real(kind=r_tran),   intent(in)    :: field_in(undf_w2h)
  real(kind=r_def),    intent(in)    :: w2_divisor(undf_w2)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_fs2d)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_fs2d)

  ! Internal variables
  integer(kind=i_def) :: df, j, w2h_idx, w2_idx, nvert
  integer(kind=i_def) :: hori_dofs_to_do
  logical(kind=l_def) :: lowest_order

  ! Number of vertical levels to loop over (W2H is layer-based, like W3)
  nvert = nlayers - 1

  ! At lowest order (ndf_w2 == 6) W2H faces are shared between neighbouring
  ! columns, so use face_selector_ew/ns to visit each face only once (as
  ! split_w2_field_kernel_mod does). At higher order there is no such
  ! duplication, so loop over all W2H DoFs directly.
  lowest_order = (ndf_w2 == 6)

  if (lowest_order) then
    hori_dofs_to_do = ABS(face_selector_ew(map_fs2d(1))) + ABS(face_selector_ns(map_fs2d(1)))
  else
    hori_dofs_to_do = ndf_w2h
  end if

  do j = 1, hori_dofs_to_do
    if (lowest_order) then
      df = face_from_face_selector(j, face_selector_ew(map_fs2d(1)), face_selector_ns(map_fs2d(1)))
    else
      df = j
    end if

    w2h_idx = map_w2h(df)
    w2_idx = map_w2(df)

    field_out(w2h_idx : w2h_idx+nvert) = field_in(w2h_idx : w2h_idx+nvert)     &
                                         / w2_divisor(w2_idx : w2_idx+nvert)
  end do

end subroutine divide_w2h_by_w2_code

end module divide_w2h_by_w2_kernel_mod
