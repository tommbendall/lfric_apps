!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! GitHub Copilot.
!-------------------------------------------------------------------------------
!> @brief Computes the "pure PEREGRIN wind" at cubed-sphere panel edges.
!> @details The PEREGRIN wind is the fixed-in-time, purely geometric velocity
!!          that maps a native equiangular coordinate onto its corresponding
!!          extended-mesh coordinate across a cubed-sphere panel edge. Following
!!          the W2H convention used throughout the transport code, a horizontal
!!          wind field is stored as \f$ \Delta t \, \mathrm{d}A \, u \f$, i.e. a
!!          volume swept through a face during a time step. This kernel therefore
!!          computes and stores the volume swept on a row of cells between a
!!          coordinate and its corresponding extended-mesh coordinate.
!!
!!          The swept volume is accumulated cell-by-cell, walking outwards from
!!          the face towards the panel edge, in exactly the same way that
!!          horizontal FFSL departure distances accumulate mass in
!!          hori_dep_dist_ffsl_1d, using a linear interpolation in the final
!!          ("departure") cell.
!!
!!          Away from panel edges (where the panel edge distance is FAR_AWAY)
!!          no computation is performed and the wind is left as zero.
module peregrin_wind_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type, func_type,                          &
                                 GH_FIELD,                                     &
                                 GH_REAL, GH_INTEGER,                          &
                                 GH_READ, GH_WRITE,                            &
                                 GH_BASIS, GH_EVALUATOR,                       &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_7,                    &
                                 STENCIL, CROSS2D, CELL_COLUMN
use fs_continuity_mod,     only: W2h, W3
use constants_mod,         only: r_def, r_tran, i_def, l_def

! Configuration modules
use base_mesh_config_mod,      only: geometry, topology
use finite_element_config_mod, only: coord_system
use planet_config_mod,         only: scaled_radius

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!! PSy layer
type, public, extends(kernel_type) :: peregrin_wind_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                          &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2h),                        & ! extended_panel_swept_volume
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3, STENCIL(CROSS2D)),       & ! cell_volume
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5,   &
                                                      STENCIL(CROSS2D)),       & ! chi (native alpha,beta,r)
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7),  & ! panel_edge_coords_x
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7),  & ! panel_edge_coords_y
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,   &
                                                        STENCIL(CROSS2D)),     & ! panel_id
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  & ! panel_edge_dist
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  & ! face_selector_ew
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)   & ! face_selector_ns
  /)
  !> Basis functions are evaluated at the W2H nodes (the write field) so that
  !! chi and the extended coordinates can be interpolated onto the faces
  type(func_type) :: meta_funcs(2) = (/                                        &
       func_type(ANY_DISCONTINUOUS_SPACE_5, GH_BASIS),                         & ! chi
       func_type(ANY_DISCONTINUOUS_SPACE_7, GH_BASIS)                          & ! extended coords
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: peregrin_wind_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: peregrin_wind_code
public :: peregrin_wind_1d

contains

!> @brief Computes the extended-panel swept volume at cubed-sphere panel edges.
!> @param[in]     nlayers            Number of layers in the mesh
!> @param[in,out] extended_panel_swept_volume W2H field storing the swept
!!                    volume between the native and extended-mesh coordinates
!> @param[in]     cell_volume        W3 field of cell volumes (the swept volume
!!                                   accumulated across the row of cells)
!> @param[in]     cv_stencil_size    Sizes of the branches of the cell-volume
!!                                   cross stencil
!> @param[in]     cv_max_length      Maximum cell-volume cross-stencil branch
!> @param[in]     cv_stencil_map     DoF map for the cell-volume cross stencil
!> @param[in]     chi_1              First (alpha) native coordinate field
!> @param[in]     chi_2              Second (beta) native coordinate field
!> @param[in]     chi_3              Third (radius) native coordinate field
!> @param[in]     chi_stencil_size   Sizes of the branches of the chi stencil
!> @param[in]     chi_max_length     Maximum chi cross-stencil branch length
!> @param[in]     chi_stencil_map    DoF map for the chi cross stencil
!> @param[in]     alpha_x            Extended alpha coordinate across the x edge
!> @param[in]     beta_x             Extended beta coordinate across the x edge
!> @param[in]     alpha_y            Extended alpha coordinate across the y edge
!> @param[in]     beta_y             Extended beta coordinate across the y edge
!> @param[in]     panel_id           Field of cubed-sphere panel IDs
!> @param[in]     pid_stencil_size   Sizes of the branches of the panel_id stencil
!> @param[in]     pid_max_length     Maximum panel_id cross-stencil branch length
!> @param[in]     pid_stencil_map    DoF map for the panel_id cross stencil
!> @param[in]     panel_edge_dist_W  Distance of each column from the W panel edge
!> @param[in]     panel_edge_dist_E  Distance of each column from the E panel edge
!> @param[in]     panel_edge_dist_S  Distance of each column from the S panel edge
!> @param[in]     panel_edge_dist_N  Distance of each column from the N panel edge
!> @param[in]     face_selector_ew   Indicates which x faces to loop over
!> @param[in]     face_selector_ns   Indicates which y faces to loop over
!> @param[in]     ndf_w2h            Num DoFs per cell for W2H
!> @param[in]     undf_w2h           Num W2H DoFs in memory for this partition
!> @param[in]     map_w2h            Map of lowest-cell W2H DoFs
!> @param[in]     ndf_w3             Num DoFs per cell for W3
!> @param[in]     undf_w3            Num W3 DoFs in memory for this partition
!> @param[in]     map_w3             Map of lowest-cell W3 DoFs
!> @param[in]     ndf_wx             Num DoFs per cell for coordinate fields
!> @param[in]     undf_wx            Num coordinate DoFs for this partition
!> @param[in]     map_wx             Map of lowest-cell coordinate DoFs
!> @param[in]     basis_wx           Basis functions of the chi space evaluated
!!                                   at the W2H nodes
!> @param[in]     ndf_wx_2d          Num DoFs per cell for 2D coord fields
!> @param[in]     undf_wx_2d         Num 2D coord DoFs for this partition
!> @param[in]     map_wx_2d          Map for 2D coordinate fields
!> @param[in]     basis_wx_2d        Basis functions of the extended-coordinate
!!                                   space evaluated at the W2H nodes
!> @param[in]     ndf_pid            Num DoFs per cell for panel ID
!> @param[in]     undf_pid           Num panel ID DoFs for this partition
!> @param[in]     map_pid            Map for panel ID field
subroutine peregrin_wind_code( nlayers,                                        &
                               extended_panel_swept_volume,                    &
                               cell_volume,                                    &
                               cv_stencil_size, cv_max_length, cv_stencil_map, &
                               chi_1, chi_2, chi_3,                            &
                               chi_stencil_size, chi_max_length,               &
                               chi_stencil_map,                                &
                               alpha_x, beta_x, alpha_y, beta_y,               &
                               panel_id,                                       &
                               pid_stencil_size, pid_max_length,               &
                               pid_stencil_map,                                &
                               panel_edge_dist_W, panel_edge_dist_E,           &
                               panel_edge_dist_S, panel_edge_dist_N,           &
                               face_selector_ew, face_selector_ns,             &
                               ndf_w2h, undf_w2h, map_w2h,                     &
                               ndf_w3, undf_w3, map_w3,                        &
                               ndf_wx, undf_wx, map_wx, basis_wx,              &
                               ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,  &
                               ndf_pid, undf_pid, map_pid                      &
                             )

  use panel_edge_support_mod, only: panel_neighbour, FAR_AWAY
  use reference_element_mod,  only: W, S, N, E

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: cv_stencil_size(4)
  integer(kind=i_def), intent(in)    :: chi_stencil_size(4)
  integer(kind=i_def), intent(in)    :: pid_stencil_size(4)
  integer(kind=i_def), intent(in)    :: cv_max_length, chi_max_length
  integer(kind=i_def), intent(in)    :: pid_max_length
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: cv_stencil_map(ndf_w3, cv_max_length, 4)
  integer(kind=i_def), intent(in)    :: chi_stencil_map(ndf_wx, chi_max_length, 4)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid, pid_max_length, 4)
  real(kind=r_tran),   intent(inout) :: extended_panel_swept_volume(undf_w2h)
  real(kind=r_tran),   intent(in)    :: cell_volume(undf_w3)
  real(kind=r_def),    intent(in)    :: chi_1(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_2(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_3(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: alpha_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_pid)
  real(kind=r_def),    intent(in)    :: basis_wx(1, ndf_wx, ndf_w2h)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_w2h)

  ! Internal variables
  integer(kind=i_def) :: pid_idx, owned_panel
  integer(kind=i_def) :: swapped_panel_x, swapped_panel_y

  pid_idx = map_pid(1)

  swapped_panel_x = 0
  swapped_panel_y = 0

  ! Nothing to do away from panel edges: the panel edge distance is FAR_AWAY
  if (ABS(panel_edge_dist_W(pid_idx)) >= FAR_AWAY .and.                        &
      ABS(panel_edge_dist_E(pid_idx)) >= FAR_AWAY .and.                        &
      ABS(panel_edge_dist_S(pid_idx)) >= FAR_AWAY .and.                        &
      ABS(panel_edge_dist_N(pid_idx)) >= FAR_AWAY) then
    return
  end if

  owned_panel = int(panel_id(pid_idx), i_def)

  ! Panel edge in the x direction: compute the change in the beta coordinate,
  ! walking stencils in the y (S-N) direction and storing at the S/N faces ---- !
  if (ABS(panel_edge_dist_W(pid_idx)) < FAR_AWAY .or.                          &
      ABS(panel_edge_dist_E(pid_idx)) < FAR_AWAY) then
    if (ABS(panel_edge_dist_W(pid_idx)) < ABS(panel_edge_dist_E(pid_idx))) then
      swapped_panel_x = panel_neighbour(owned_panel, W)
    else
      swapped_panel_x = panel_neighbour(owned_panel, E)
    end if

    if (swapped_panel_x /= 0_i_def) then
      call peregrin_wind_1d( nlayers, .true.,                                  &
                             extended_panel_swept_volume, cell_volume,         &
                             cv_stencil_size, cv_max_length, cv_stencil_map,   &
                             chi_1, chi_2, chi_3,                              &
                             chi_stencil_size, chi_max_length, chi_stencil_map,&
                             alpha_x, beta_x,                                  &
                             panel_id, pid_stencil_size, pid_max_length,       &
                             pid_stencil_map,                                  &
                             owned_panel, swapped_panel_x,                     &
                             face_selector_ns,                                 &
                             ndf_w2h, undf_w2h, map_w2h,                       &
                             ndf_w3, undf_w3, map_w3,                          &
                             ndf_wx, undf_wx, map_wx, basis_wx,                &
                             ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,    &
                             ndf_pid, undf_pid, map_pid )
    end if
  end if

  ! Panel edge in the y direction: compute the change in the alpha coordinate,
  ! walking stencils in the x (W-E) direction and storing at the W/E faces ---- !
  if (ABS(panel_edge_dist_S(pid_idx)) < FAR_AWAY .or.                          &
      ABS(panel_edge_dist_N(pid_idx)) < FAR_AWAY) then
    if (ABS(panel_edge_dist_S(pid_idx)) < ABS(panel_edge_dist_N(pid_idx))) then
      swapped_panel_y = panel_neighbour(owned_panel, S)
    else
      swapped_panel_y = panel_neighbour(owned_panel, N)
    end if

    if (swapped_panel_y /= 0_i_def) then
      call peregrin_wind_1d( nlayers, .false.,                                 &
                             extended_panel_swept_volume, cell_volume,         &
                             cv_stencil_size, cv_max_length, cv_stencil_map,   &
                             chi_1, chi_2, chi_3,                              &
                             chi_stencil_size, chi_max_length, chi_stencil_map,&
                             alpha_y, beta_y,                                  &
                             panel_id, pid_stencil_size, pid_max_length,       &
                             pid_stencil_map,                                  &
                             owned_panel, swapped_panel_y,                     &
                             face_selector_ew,                                 &
                             ndf_w2h, undf_w2h, map_w2h,                       &
                             ndf_w3, undf_w3, map_w3,                          &
                             ndf_wx, undf_wx, map_wx, basis_wx,                &
                             ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,    &
                             ndf_pid, undf_pid, map_pid )
    end if
  end if

end subroutine peregrin_wind_code

! ============================================================================ !
! SINGLE UNDERLYING 1D ROUTINE
! ============================================================================ !

!> @brief Computes the extended-panel swept volume in one direction (x or y).
!> @details For each selected face, chi and the extended-mesh coordinate are
!!          evaluated at the W2H point using the basis functions, and their
!!          horizontal displacement (in alpha or beta) is found. The branch of
!!          the cross stencil to walk is chosen from the sign of the
!!          displacement, so the row of cells is always swept towards the
!!          extended-mesh coordinate. Whole cell volumes are accumulated until
!!          the native coordinate spans this displacement, with a linear
!!          interpolation in the final ("departure") cell. The resulting swept
!!          volume is stored per layer in the W2H field.
subroutine peregrin_wind_1d( nlayers, beta_wind,                               &
                             extended_panel_swept_volume, cell_volume,         &
                             cv_stencil_size, cv_max_length, cv_stencil_map,   &
                             chi_1, chi_2, chi_3,                              &
                             chi_stencil_size, chi_max_length, chi_stencil_map,&
                             alpha_ext, beta_ext,                              &
                             panel_id, pid_stencil_size, pid_max_length,       &
                             pid_stencil_map,                                  &
                             owned_panel, swapped_panel,                       &
                             face_selector,                                    &
                             ndf_w2h, undf_w2h, map_w2h,                       &
                             ndf_w3, undf_w3, map_w3,                          &
                             ndf_wx, undf_wx, map_wx, basis_wx,                &
                             ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,    &
                             ndf_pid, undf_pid, map_pid )

  use sci_chi_transform_mod, only: chi2xyz, get_to_rotate, get_to_stretch,     &
                                   get_inverse_mesh_rotation_matrix,           &
                                   get_stretch_factor
  use coord_transform_mod,   only: alphabetar2xyz, xyz2alphabetar,             &
                                   inverse_schmidt_transform_xyz
  use reference_element_mod, only: W, S, N, E

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  logical(kind=l_def), intent(in)    :: beta_wind
  integer(kind=i_def), intent(in)    :: cv_stencil_size(4)
  integer(kind=i_def), intent(in)    :: chi_stencil_size(4)
  integer(kind=i_def), intent(in)    :: pid_stencil_size(4)
  integer(kind=i_def), intent(in)    :: cv_max_length, chi_max_length
  integer(kind=i_def), intent(in)    :: pid_max_length
  integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: cv_stencil_map(ndf_w3, cv_max_length, 4)
  integer(kind=i_def), intent(in)    :: chi_stencil_map(ndf_wx, chi_max_length, 4)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid, pid_max_length, 4)
  integer(kind=i_def), intent(in)    :: owned_panel, swapped_panel
  real(kind=r_tran),   intent(inout) :: extended_panel_swept_volume(undf_w2h)
  real(kind=r_tran),   intent(in)    :: cell_volume(undf_w3)
  real(kind=r_def),    intent(in)    :: chi_1(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_2(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_3(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: face_selector(undf_pid)
  real(kind=r_def),    intent(in)    :: basis_wx(1, ndf_wx, ndf_w2h)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_w2h)

  ! Local parameters
  integer(kind=i_def), parameter :: interp_dir_alpha = 1
  integer(kind=i_def), parameter :: interp_dir_beta  = 2
  real(kind=r_def),    parameter :: unit_radius = 1.0_r_def

  ! Local variables
  integer(kind=i_def) :: local_dofs(2)
  integer(kind=i_def) :: df, df_idx, dfc, branch
  integer(kind=i_def) :: interp_dir, panel_edge
  integer(kind=i_def) :: j, k, m, nrow, w2h_idx, cell_idx
  real(kind=r_def)    :: xyz(3), abh(3), h_dummy
  real(kind=r_def)    :: alpha_owned, beta_owned
  real(kind=r_def)    :: x_ext, displacement
  real(kind=r_def)    :: cell_width, remaining, frac
  real(kind=r_def)    :: direction
  real(kind=r_def)    :: stretch_factor
  real(kind=r_def)    :: inverse_rot_matrix(3,3)
  real(kind=r_def)    :: x_native_stencil(chi_max_length)
  real(kind=r_def)    :: x_native_this_face
  integer(kind=i_def) :: ncell
  integer(kind=i_def) :: swept_cell_idx(cv_max_length)
  real(kind=r_def)    :: swept_cell_width(cv_max_length)
  logical(kind=l_def) :: to_rotate, to_stretch
  real(kind=r_tran)   :: swept_volume(nlayers)

  to_rotate = get_to_rotate()
  to_stretch = get_to_stretch()
  if (to_rotate) then
    inverse_rot_matrix = get_inverse_mesh_rotation_matrix()
  end if
  if (to_stretch) then
    stretch_factor = get_stretch_factor()
  end if

  ! Select the faces to store at and the stencil branches to walk ------------ !
  ! A panel edge in the x direction gives a change in the beta coordinate, so
  ! the row of cells runs in the y (S-N) direction and the result is stored at
  ! the S/N faces. A panel edge in the y direction gives a change in alpha, so
  ! the row runs in the x (W-E) direction, stored at the W/E faces.
  if (beta_wind) then
    local_dofs = (/ S, N /)
    ! The stored y-direction (beta) wind is negative: a positive value means
    ! wind blowing N -> S
    direction = -1.0_r_def
  else
    local_dofs = (/ W, E /)
    ! The stored x-direction (alpha) wind follows the usual sign convention
    direction = 1.0_r_def
  end if

  ! Determine whether the displacement is measured in alpha or beta ---------- !
  panel_edge = 10*owned_panel + swapped_panel
  select case (panel_edge)
  case (15, 26, 34, 43, 51, 62)
    interp_dir = interp_dir_alpha
  case (12, 21, 36, 45, 54, 63)
    interp_dir = interp_dir_beta
  case (16, 41, 32, 25, 64, 53)
    interp_dir = interp_dir_alpha
  case (61, 14, 23, 52, 46, 35)
    interp_dir = interp_dir_beta
  case default
    interp_dir = interp_dir_beta
  end select

  ! Loop over the selected faces --------------------------------------------- !
  do df_idx = 1, ABS(face_selector(map_pid(1)))
    df = local_dofs(df_idx - MIN(0, face_selector(map_pid(1))))

    w2h_idx = map_w2h(df)

    ! ------------------------------------------------------------------------ !
    ! Extended-mesh coordinate evaluated at the stored face of the owned cell
    ! ------------------------------------------------------------------------ !
    ! The extended coordinate is interpolated to the W2H point using its basis
    ! functions, then continued from the neighbouring panel onto the owned one.
    abh(:) = 0.0_r_def
    do dfc = 1, ndf_wx_2d
      abh(interp_dir_alpha) = abh(interp_dir_alpha)                            &
                            + alpha_ext(map_wx_2d(dfc)) * basis_wx_2d(1, dfc, df)
      abh(interp_dir_beta)  = abh(interp_dir_beta)                             &
                            + beta_ext(map_wx_2d(dfc)) * basis_wx_2d(1, dfc, df)
    end do
    call alphabetar2xyz( abh(interp_dir_alpha), abh(interp_dir_beta),          &
                         unit_radius, swapped_panel, xyz(1), xyz(2), xyz(3) )
    call xyz2alphabetar( xyz(1), xyz(2), xyz(3), owned_panel,                  &
                         abh(1), abh(2), h_dummy )
    x_ext = abh(interp_dir)

    ! ------------------------------------------------------------------------ !
    ! Native coordinate at the stored face of the owned cell
    ! ------------------------------------------------------------------------ !
    ! The chi field is evaluated at the stored face's W2H point, after each DoF
    ! has been transformed into the owned panel's native equiangular system.
    abh(:) = 0.0_r_def
    do dfc = 1, ndf_wx
      call chi2xyz( chi_1(chi_stencil_map(dfc, 1, df)),                        &
                    chi_2(chi_stencil_map(dfc, 1, df)),                        &
                    unit_radius, owned_panel,                                  &
                    geometry, topology, coord_system, scaled_radius,           &
                    xyz(1), xyz(2), xyz(3) )
      if (to_rotate) xyz = matmul(inverse_rot_matrix, xyz)
      if (to_stretch) xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
      call xyz2alphabetar( xyz(1), xyz(2), xyz(3), owned_panel,                &
                           alpha_owned, beta_owned, h_dummy )
      abh(interp_dir_alpha) = abh(interp_dir_alpha)                            &
                            + alpha_owned * basis_wx(1, dfc, df)
      abh(interp_dir_beta)  = abh(interp_dir_beta)                             &
                            + beta_owned * basis_wx(1, dfc, df)
    end do
    x_native_this_face = abh(interp_dir)

    ! Displacement to be swept (native minus extended), at the stored face
    displacement = x_native_this_face - x_ext

    ! ------------------------------------------------------------------------ !
    ! Select the branch of the cross stencil to walk, from the sign of the
    ! displacement. In the owned panel's native frame the coordinate increases
    ! towards local_dofs(2) (N or E), so a positive displacement (native larger
    ! than extended) is swept towards local_dofs(1) (S or W) and vice versa.
    ! The swept row can therefore run in either direction, independently of
    ! which face the wind is stored at.
    ! ------------------------------------------------------------------------ !
    if (displacement >= 0.0_r_def) then
      branch = local_dofs(1)
    else
      branch = local_dofs(2)
    end if
    nrow = cv_stencil_size(branch)

    ! ------------------------------------------------------------------------ !
    ! Native coordinate at the branch face of each cell along the swept row
    ! ------------------------------------------------------------------------ !
    do j = 1, nrow
      abh(:) = 0.0_r_def
      do dfc = 1, ndf_wx
        call chi2xyz( chi_1(chi_stencil_map(dfc, j, branch)),                  &
                      chi_2(chi_stencil_map(dfc, j, branch)),                  &
                      unit_radius, owned_panel,                                &
                      geometry, topology, coord_system, scaled_radius,         &
                      xyz(1), xyz(2), xyz(3) )
        if (to_rotate) xyz = matmul(inverse_rot_matrix, xyz)
        if (to_stretch) xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
        call xyz2alphabetar( xyz(1), xyz(2), xyz(3), owned_panel,              &
                             alpha_owned, beta_owned, h_dummy )
        abh(interp_dir_alpha) = abh(interp_dir_alpha)                          &
                              + alpha_owned * basis_wx(1, dfc, branch)
        abh(interp_dir_beta)  = abh(interp_dir_beta)                           &
                              + beta_owned * basis_wx(1, dfc, branch)
      end do
      x_native_stencil(j) = abh(interp_dir)
    end do

    ! ------------------------------------------------------------------------ !
    ! Assemble the row of cells to sweep, and their coordinate widths
    ! ------------------------------------------------------------------------ !
    ! If the branch runs in the stored face's own direction the row is simply
    ! the branch's neighbours. If it runs the opposite way, the owned cell is
    ! swept first (it lies between the stored face and the branch face of the
    ! owned cell), followed by the branch's neighbours.
    if (branch == df) then
      ncell = nrow - 1
      do m = 1, ncell
        swept_cell_idx(m) = cv_stencil_map(1, m+1, branch)
        swept_cell_width(m) = ABS(x_native_stencil(m+1) - x_native_stencil(m))
      end do
    else
      ncell = nrow
      swept_cell_idx(1) = cv_stencil_map(1, 1, branch)
      swept_cell_width(1) = ABS(x_native_stencil(1) - x_native_this_face)
      do m = 1, nrow - 1
        swept_cell_idx(m+1) = cv_stencil_map(1, m+1, branch)
        swept_cell_width(m+1) = ABS(x_native_stencil(m+1) - x_native_stencil(m))
      end do
    end if

    ! ------------------------------------------------------------------------ !
    ! Accumulate the swept volume across the row of cells
    ! ------------------------------------------------------------------------ !
    swept_volume(:) = 0.0_r_tran
    remaining = ABS(displacement)

    do m = 1, ncell

      cell_idx = swept_cell_idx(m)
      cell_width = swept_cell_width(m)
      if (cell_width <= 0.0_r_def) EXIT

      if (remaining >= cell_width) then
        ! Whole cell is swept
        do k = 1, nlayers
          swept_volume(k) = swept_volume(k) + cell_volume(cell_idx + k - 1)
        end do
        remaining = remaining - cell_width
      else
        ! Linear interpolation in the final ("departure") cell
        frac = remaining / cell_width
        do k = 1, nlayers
          swept_volume(k) = swept_volume(k)                                    &
                          + REAL(frac, r_tran) * cell_volume(cell_idx + k - 1)
        end do
        remaining = 0.0_r_def
        EXIT
      end if
    end do

    ! Store the swept volume, signed by the displacement and the direction
    ! convention (positive beta wind means flow from N to S)
    do k = 1, nlayers
      extended_panel_swept_volume(w2h_idx + k - 1) =                           &
          REAL(SIGN(1.0_r_def, direction*displacement), r_tran) * swept_volume(k)
    end do

  end do  ! faces

end subroutine peregrin_wind_1d

end module peregrin_wind_kernel_mod
