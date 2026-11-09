!----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of UM clock

module um_clock_init_mod

  ! LFRic modules
  use clock_mod,     only : clock_type
  use constants_mod, only : r_def, r_um

  implicit none

  private
  public :: um_clock_init
  public :: set_um_timestep

contains

  !> @brief Initialise UM clock
  !> @param[in] clock  The LFRic model clock object
  subroutine um_clock_init( clock )

    ! UM modules containing things that need setting
    use timestep_mod, only: timestep, recip_timestep

    implicit none

    class(clock_type), intent(in) :: clock

    ! Timestep used in UM code - contained in UM timestep_mod.
    ! Set from LFRic input timestep.
    timestep = real(clock%get_seconds_per_step(), r_um)
    recip_timestep = 1.0_r_um / timestep

  end subroutine um_clock_init

  !> @brief Set the timestep used in the UM code to a different length
  !> @details The microphysics constants depend on the timestep, so are
  !!          recalculated at the next call of the microphysics. Call this from
  !!          the algorithm layer, outside any OpenMP region, and restore the
  !!          timestep of the model afterwards.
  !> @param[in] dt  The timestep length (s)
  subroutine set_um_timestep( dt )

    use timestep_mod, only: timestep, recip_timestep
    use mphys_constants_mod, only: l_calc_mp_const

    implicit none

    real(r_def), intent(in) :: dt

    timestep = real(dt, r_um)
    recip_timestep = 1.0_r_um / timestep
    l_calc_mp_const = .true.

  end subroutine set_um_timestep

end module um_clock_init_mod
