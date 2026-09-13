module fusion_fortran
  !! ISO_C_BINDING wrapper for the stateless event-rate particle network.
  !!
  !! Event rates are supplied in m^-3 s^-1.  The wrapper exposes separate
  !! nonnegative reactant losses, product births, net sources, and neutron
  !! births.  Product births are nuclear births and do not imply thermalization.
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int
  integer(c_int), parameter, public :: PB11_METHOD_INTEGRAL = 0_c_int
  integer(c_int), parameter, public :: PB11_METHOD_FAST = 1_c_int

  integer(c_int), parameter, public :: FUSION_PROTON = 0_c_int
  integer(c_int), parameter, public :: FUSION_DEUTERON = 1_c_int
  integer(c_int), parameter, public :: FUSION_TRITON = 2_c_int
  integer(c_int), parameter, public :: FUSION_HELIUM3 = 3_c_int
  integer(c_int), parameter, public :: FUSION_HELIUM4 = 4_c_int
  integer(c_int), parameter, public :: FUSION_BORON11 = 5_c_int
  integer(c_int), parameter, public :: FUSION_SPECIES_COUNT = 6_c_int

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = 0_c_int
  integer(c_int), parameter, public :: FUSION_DD_TP = 1_c_int
  integer(c_int), parameter, public :: FUSION_DD_HE3N = 2_c_int
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = 3_c_int
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = 4_c_int
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = 5_c_int

  ! Public IDs follow the C ABI and start at zero.  Fortran array positions
  ! corresponding to an ID are therefore ID + 1.

  type, bind(C), public :: fusion_particle_sources_v1
    real(c_double) :: reactant_loss(6)
    real(c_double) :: product_birth(6)
    real(c_double) :: net_source(6)
    real(c_double) :: neutron_birth
  end type fusion_particle_sources_v1

  type, bind(C), public :: fusion_thermal_rates_v1
    real(c_double) :: reactivity_m3_s(5)
    real(c_double) :: event_rate_m3_s(5)
    type(fusion_particle_sources_v1) :: particles
    real(c_double) :: nuclear_power_W_m3(5)
    real(c_double) :: total_nuclear_power_W_m3
  end type fusion_thermal_rates_v1

  public :: fusion_particle_sources
  public :: fusion_channel_q
  public :: fusion_thermal_reactivity
  public :: fusion_thermal_rates
  public :: fusion_cross_section

  interface
     function c_fusion_particle_sources(event_rates, source) &
          bind(C, name="fusion_c_particle_sources") result(status)
       import :: c_double, c_int, fusion_particle_sources_v1
       real(c_double), intent(in) :: event_rates(5)
       type(fusion_particle_sources_v1), intent(out) :: source
       integer(c_int) :: status
     end function c_fusion_particle_sources

     function c_fusion_channel_q(channel, q_j) &
          bind(C, name="fusion_c_channel_q") result(status)
       import :: c_double, c_int
       integer(c_int), value :: channel
       real(c_double) :: q_j
       integer(c_int) :: status
     end function c_fusion_channel_q

     function c_fusion_thermal_reactivity(channel, kt_j, pb_method, value) &
          bind(C, name="fusion_c_thermal_reactivity") result(status)
       import :: c_double, c_int
       integer(c_int), value :: channel
       real(c_double), value :: kt_j
       integer(c_int), value :: pb_method
       real(c_double) :: value
       integer(c_int) :: status
     end function c_fusion_thermal_reactivity

     function c_fusion_thermal_rates(kt_j, density, channel_mask, pb_method, &
          rates) bind(C, name="fusion_c_thermal_rates") result(status)
       import :: c_double, c_int, fusion_thermal_rates_v1
       real(c_double), value :: kt_j
       real(c_double), intent(in) :: density(6)
       integer(c_int), value :: channel_mask
       integer(c_int), value :: pb_method
       type(fusion_thermal_rates_v1), intent(out) :: rates
       integer(c_int) :: status
     end function c_fusion_thermal_rates

     function c_fusion_cross_section(channel, relative_energy_j, value) &
          bind(C, name="fusion_c_cross_section") result(status)
       import :: c_double, c_int
       integer(c_int), value :: channel
       real(c_double), value :: relative_energy_j
       real(c_double) :: value
       integer(c_int) :: status
     end function c_fusion_cross_section
  end interface

contains

  subroutine fusion_particle_sources(event_rates, source, status)
    real(c_double), intent(in) :: event_rates(5)
    type(fusion_particle_sources_v1), intent(out) :: source
    integer(c_int), intent(out) :: status

    source%reactant_loss = 0.0_c_double
    source%product_birth = 0.0_c_double
    source%net_source = 0.0_c_double
    source%neutron_birth = 0.0_c_double
    status = c_fusion_particle_sources(event_rates, source)
  end subroutine fusion_particle_sources

  subroutine fusion_channel_q(channel, q_j, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(out) :: q_j
    integer(c_int), intent(out) :: status

    q_j = 0.0_c_double
    status = c_fusion_channel_q(channel, q_j)
  end subroutine fusion_channel_q

  subroutine fusion_thermal_reactivity(channel, kt_j, pb_method, value, &
       status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: kt_j
    integer(c_int), intent(in) :: pb_method
    real(c_double), intent(out) :: value
    integer(c_int), intent(out) :: status

    value = 0.0_c_double
    status = c_fusion_thermal_reactivity(channel, kt_j, pb_method, value)
  end subroutine fusion_thermal_reactivity

  subroutine fusion_thermal_rates(kt_j, density, channel_mask, pb_method, &
       rates, status)
    real(c_double), intent(in) :: kt_j
    real(c_double), intent(in) :: density(6)
    integer(c_int), intent(in) :: channel_mask
    integer(c_int), intent(in) :: pb_method
    type(fusion_thermal_rates_v1), intent(out) :: rates
    integer(c_int), intent(out) :: status

    call clear_thermal_rates(rates)
    ! The frozen C ABI uses int for this bit mask; accept only 0..31.
    if (channel_mask < 0_c_int .or. channel_mask > 31_c_int) then
       status = PB11_STATUS_INVALID_ARGUMENT
       return
    end if
    status = c_fusion_thermal_rates(kt_j, density, channel_mask, pb_method, &
         rates)
  end subroutine fusion_thermal_rates

  subroutine fusion_cross_section(channel, relative_energy_j, value, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: relative_energy_j
    real(c_double), intent(out) :: value
    integer(c_int), intent(out) :: status

    value = 0.0_c_double
    status = c_fusion_cross_section(channel, relative_energy_j, value)
  end subroutine fusion_cross_section

  subroutine clear_particle_sources(source)
    type(fusion_particle_sources_v1), intent(out) :: source

    source%reactant_loss = 0.0_c_double
    source%product_birth = 0.0_c_double
    source%net_source = 0.0_c_double
    source%neutron_birth = 0.0_c_double
  end subroutine clear_particle_sources

  subroutine clear_thermal_rates(rates)
    type(fusion_thermal_rates_v1), intent(out) :: rates

    rates%reactivity_m3_s = 0.0_c_double
    rates%event_rate_m3_s = 0.0_c_double
    call clear_particle_sources(rates%particles)
    rates%nuclear_power_W_m3 = 0.0_c_double
    rates%total_nuclear_power_W_m3 = 0.0_c_double
  end subroutine clear_thermal_rates

end module fusion_fortran
