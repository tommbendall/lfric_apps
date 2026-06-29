!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the energy flux associated with moisture entering/leaving
!>        the atmosphere at the surface.
!>
!> @details When moisture evaporates into the atmosphere, it adds mass that
!>   carries kinetic and internal energy. When precipitation leaves, it removes
!>   mass that was carrying kinetic and internal energy. These energy changes
!>   are not captured by the latent heat terms (Lvr*rain, Lsr*snow) already
!>   tracked in sum_fluxes_alg, so this kernel computes the additional
!>   contribution to the energy budget.
!>
!>   For evaporation (vapour entering):
!>     flux += [(cpv - Rv)*T_surf + 0.5*(u^2 + v^2)] * evap_rate
!>
!>   For rain leaving:
!>     flux -= [cl*T_surf + 0.5*(u^2 + v^2)] * rain_rate
!>
!>   For snow/graupel leaving:
!>     flux -= [ci*T_surf + 0.5*(u^2 + v^2)] * ice_rate
!>
!>   SIMPLIFICATION: The temperature used for the heat capacity with
!>   precipitation is the surface temperature and NOT the temperature at
!>   which the precipitation "vanished" from the model grid. This is a
!>   significant approximation for rain that forms aloft and falls to the
!>   surface, but avoids the complexity of tracking where precipitation
!>   originates.
!>
module moisture_energy_flux_kernel_mod

  use argument_mod,  only: arg_type, CELL_COLUMN,                              &
                           GH_FIELD, GH_SCALAR,                                &
                           GH_REAL, GH_READ, GH_READWRITE,                     &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: moisture_energy_flux_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                       &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ),                                 &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: moisture_energy_flux_code
  end type moisture_energy_flux_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: moisture_energy_flux_code

contains

  !> @brief Compute the energy flux from moisture entering/leaving atmosphere
  !! @param[in]     nlayers           Number of layers
  !! @param[in,out] accumulated_fluxes Accumulated energy flux field
  !! @param[in]     theta_surf        Potential temperature at surface
  !! @param[in]     exner_surf        Exner pressure at surface
  !! @param[in]     u_surf            Zonal wind at lowest level
  !! @param[in]     v_surf            Meridional wind at lowest level
  !! @param[in]     vap_in            Surface moisture flux (evaporation, kg/m2/s)
  !! @param[in]     liq_out           Total liquid precipitation rate (kg/m2/s)
  !! @param[in]     ice_out           Total ice precipitation rate (kg/m2/s)
  !! @param[in]     cpv_local         Heat capacity of water vapour at const. p
  !! @param[in]     cl_local          Heat capacity of liquid water
  !! @param[in]     ci_local          Heat capacity of ice
  !! @param[in]     Rv_local          Gas constant of water vapour
  !! @param[in]     ndf               Number of DOFs per cell
  !! @param[in]     undf              Total number of unique DOFs
  !! @param[in]     map               DOF map for the cell
  subroutine moisture_energy_flux_code( nlayers,                               &
                                        accumulated_fluxes,                    &
                                        theta_surf, exner_surf,                &
                                        u_surf, v_surf,                        &
                                        vap_in, liq_out, ice_out,              &
                                        cpv_local, cl_local, ci_local,         &
                                        Rv_local,                              &
                                        ndf, undf, map )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf, undf
    integer(kind=i_def), intent(in) :: map(ndf)

    real(kind=r_def), intent(inout) :: accumulated_fluxes(undf)
    real(kind=r_def), intent(in)    :: theta_surf(undf)
    real(kind=r_def), intent(in)    :: exner_surf(undf)
    real(kind=r_def), intent(in)    :: u_surf(undf)
    real(kind=r_def), intent(in)    :: v_surf(undf)
    real(kind=r_def), intent(in)    :: vap_in(undf)
    real(kind=r_def), intent(in)    :: liq_out(undf)
    real(kind=r_def), intent(in)    :: ice_out(undf)
    real(kind=r_def), intent(in)    :: cpv_local
    real(kind=r_def), intent(in)    :: cl_local
    real(kind=r_def), intent(in)    :: ci_local
    real(kind=r_def), intent(in)    :: Rv_local

    ! Local variables
    integer(kind=i_def) :: k
    real(kind=r_def)    :: T_surf, ke_h, cv_vap

    ! Heat capacity at constant volume for vapour
    cv_vap = cpv_local - Rv_local

    k = map(1)

    ! Compute surface temperature
    T_surf = theta_surf(k) * exner_surf(k)

    ! Compute horizontal kinetic energy at lowest level
    ke_h = 0.5_r_def * ( u_surf(k)**2 + v_surf(k)**2 )

    ! Add energy flux contribution from moisture exchange:
    !
    ! Evaporation adds vapour carrying internal energy + KE:
    !   + (cv_vap * T_surf + ke_h) * evap_rate
    !
    ! Precipitation removes liquid/ice carrying internal energy + KE:
    !   - (cl * T_surf + ke_h) * liq_precip_rate
    !   - (ci * T_surf + ke_h) * ice_precip_rate
    !
    ! SIMPLIFICATION: The temperature used for the heat capacity of
    ! precipitation is the surface temperature, not the temperature at
    ! which the precipitation "vanished" from the model grid.
    accumulated_fluxes(k) = accumulated_fluxes(k)                              &
                          + ( cv_vap * T_surf + ke_h ) * vap_in(k)             &
                          - ( cl_local * T_surf + ke_h ) * liq_out(k)          &
                          - ( ci_local * T_surf + ke_h ) * ice_out(k)

  end subroutine moisture_energy_flux_code

end module moisture_energy_flux_kernel_mod
