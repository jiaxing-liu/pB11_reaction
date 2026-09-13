program test_fusion_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_positive_inf, &
       ieee_quiet_nan, ieee_value
  use fusion_fortran
  implicit none

  real(c_double), parameter :: mass_numbers(6) = &
       [1.0_c_double, 2.0_c_double, 3.0_c_double, 3.0_c_double, &
        4.0_c_double, 11.0_c_double]
  real(c_double), parameter :: charges(6) = &
       [1.0_c_double, 1.0_c_double, 1.0_c_double, 2.0_c_double, &
        2.0_c_double, 5.0_c_double]
  real(c_double), parameter :: joules_per_keV = 1.602176634e-16_c_double
  real(c_double), parameter :: joules_per_mb = 1.0e-31_c_double
  integer :: failures

  failures = 0

  call verify_channel(FUSION_PB11_3ALPHA)
  call verify_channel(FUSION_DD_TP)
  call verify_channel(FUSION_DD_HE3N)
  call verify_channel(FUSION_DT_ALPHAN)
  call verify_channel(FUSION_DHE3_ALPHAP)
  call verify_mixed_channels()
  call verify_zero_rates()
  call verify_errors()
  call verify_thermal_rates()
  call verify_scalar_wrappers()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran fusion-network test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran fusion-network tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function source_is_zero(source)
    type(fusion_particle_sources_v1), intent(in) :: source

    source_is_zero = all(source%reactant_loss == 0.0_c_double) .and. &
         all(source%product_birth == 0.0_c_double) .and. &
         all(source%net_source == 0.0_c_double) .and. &
         source%neutron_birth == 0.0_c_double
  end function source_is_zero

  logical function thermal_is_zero(rates)
    type(fusion_thermal_rates_v1), intent(in) :: rates

    thermal_is_zero = all(rates%reactivity_m3_s == 0.0_c_double) .and. &
         all(rates%event_rate_m3_s == 0.0_c_double) .and. &
         source_is_zero(rates%particles) .and. &
         all(rates%nuclear_power_W_m3 == 0.0_c_double) .and. &
         rates%total_nuclear_power_W_m3 == 0.0_c_double
  end function thermal_is_zero

  logical function close_relative(actual, expected, tolerance)
    real(c_double), intent(in) :: actual, expected, tolerance
    real(c_double) :: scale

    if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
       close_relative = .false.
       return
    end if
    if (actual == expected) then
       close_relative = .true.
       return
    end if
    scale = max(abs(actual), abs(expected))
    close_relative = scale > 0.0_c_double .and. &
         abs(actual - expected) <= tolerance * scale
  end function close_relative

  subroutine verify_channel(channel)
    integer(c_int), intent(in) :: channel
    real(c_double) :: rates(5), expected_loss(6), expected_birth(6)
    real(c_double) :: expected_neutron
    real(c_double) :: loss_nucleons, birth_nucleons
    real(c_double) :: loss_charge, birth_charge
    type(fusion_particle_sources_v1) :: source
    integer(c_int) :: status
    integer :: species

    rates = 0.0_c_double
    rates(channel + 1) = 1.0_c_double
    expected_loss = 0.0_c_double
    expected_birth = 0.0_c_double
    expected_neutron = 0.0_c_double

    select case (channel)
    case (FUSION_PB11_3ALPHA)
       expected_loss(FUSION_PROTON + 1) = 1.0_c_double
       expected_loss(FUSION_BORON11 + 1) = 1.0_c_double
       expected_birth(FUSION_HELIUM4 + 1) = 3.0_c_double
    case (FUSION_DD_TP)
       expected_loss(FUSION_DEUTERON + 1) = 2.0_c_double
       expected_birth(FUSION_TRITON + 1) = 1.0_c_double
       expected_birth(FUSION_PROTON + 1) = 1.0_c_double
    case (FUSION_DD_HE3N)
       expected_loss(FUSION_DEUTERON + 1) = 2.0_c_double
       expected_birth(FUSION_HELIUM3 + 1) = 1.0_c_double
       expected_neutron = 1.0_c_double
    case (FUSION_DT_ALPHAN)
       expected_loss(FUSION_DEUTERON + 1) = 1.0_c_double
       expected_loss(FUSION_TRITON + 1) = 1.0_c_double
       expected_birth(FUSION_HELIUM4 + 1) = 1.0_c_double
       expected_neutron = 1.0_c_double
    case (FUSION_DHE3_ALPHAP)
       expected_loss(FUSION_DEUTERON + 1) = 1.0_c_double
       expected_loss(FUSION_HELIUM3 + 1) = 1.0_c_double
       expected_birth(FUSION_HELIUM4 + 1) = 1.0_c_double
       expected_birth(FUSION_PROTON + 1) = 1.0_c_double
    end select

    call fusion_particle_sources(rates, source, status)
    call check(status == PB11_STATUS_OK, 'unit channel returns OK')
    do species = 1, 6
       call check(source%reactant_loss(species) == expected_loss(species), &
            'unit channel has expected reactant loss')
       call check(source%product_birth(species) == expected_birth(species), &
            'unit channel has expected product birth')
       call check(source%net_source(species) == &
            expected_birth(species) - expected_loss(species), &
            'unit channel separates birth and loss')
       call check(ieee_is_finite(source%reactant_loss(species)) .and. &
            ieee_is_finite(source%product_birth(species)) .and. &
            ieee_is_finite(source%net_source(species)), &
            'unit channel output is finite')
    end do
    call check(source%neutron_birth == expected_neutron, &
         'unit channel has expected neutron birth')

    loss_nucleons = sum(mass_numbers * expected_loss)
    birth_nucleons = sum(mass_numbers * expected_birth) + expected_neutron
    loss_charge = sum(charges * expected_loss)
    birth_charge = sum(charges * expected_birth)
    call check(loss_nucleons == birth_nucleons, &
         'unit channel conserves nucleon number including neutrons')
    call check(loss_charge == birth_charge, &
         'unit channel conserves electric charge')
  end subroutine verify_channel

  subroutine verify_mixed_channels()
    real(c_double) :: rates(5), expected_loss(6), expected_birth(6)
    type(fusion_particle_sources_v1) :: source
    integer(c_int) :: status
    integer :: species

    ! DD entries are event rates and already include their identical-pair 1/2.
    rates = [2.0_c_double, 7.0_c_double, 11.0_c_double, 13.0_c_double, &
         17.0_c_double]
    expected_loss = [2.0_c_double, 66.0_c_double, 13.0_c_double, &
         17.0_c_double, 0.0_c_double, 2.0_c_double]
    expected_birth = [24.0_c_double, 0.0_c_double, 7.0_c_double, &
         11.0_c_double, 36.0_c_double, 0.0_c_double]

    call fusion_particle_sources(rates, source, status)
    call check(status == PB11_STATUS_OK, 'mixed integer channels return OK')
    do species = 1, 6
       call check(source%reactant_loss(species) == expected_loss(species), &
            'mixed channels have known losses')
       call check(source%product_birth(species) == expected_birth(species), &
            'mixed channels have known births')
       call check(source%net_source(species) == &
            expected_birth(species) - expected_loss(species), &
            'mixed channels keep births and losses separate')
    end do
    call check(source%reactant_loss(FUSION_DEUTERON + 1) == 66.0_c_double, &
         'DD event rates are not halved a second time')
    call check(source%neutron_birth == 24.0_c_double, &
         'mixed channels count neutron births')
  end subroutine verify_mixed_channels

  subroutine verify_zero_rates()
    real(c_double) :: rates(5)
    type(fusion_particle_sources_v1) :: source
    integer(c_int) :: status

    rates = 0.0_c_double
    call fusion_particle_sources(rates, source, status)
    call check(status == PB11_STATUS_OK, 'zero rates return OK')
    call check(source_is_zero(source), 'zero rates produce zero outputs')
  end subroutine verify_zero_rates

  subroutine verify_thermal_rates()
    real(c_double) :: kt_j, density(6)
    real(c_double) :: q_tp, q_he3n
    real(c_double) :: expected_tp, expected_he3n, expected_total
    real(c_double) :: rate_tp, rate_he3n
    type(fusion_thermal_rates_v1) :: rates
    integer(c_int) :: mask, status

    kt_j = 10.0_c_double * joules_per_keV
    density = 0.0_c_double
    density(FUSION_DEUTERON + 1) = 1.0e20_c_double
    mask = ishft(1_c_int, FUSION_DD_TP) + &
         ishft(1_c_int, FUSION_DD_HE3N)

    call fusion_thermal_rates(kt_j, density, mask, PB11_METHOD_INTEGRAL, &
         rates, status)
    call check(status == PB11_STATUS_OK, 'Fortran DD-only rates return OK')
    call check(rates%reactivity_m3_s(FUSION_DD_TP + 1) > 0.0_c_double .and. &
         rates%reactivity_m3_s(FUSION_DD_HE3N + 1) > 0.0_c_double, &
         '10 keV DD branches have nonzero reactivities')

    rate_tp = rates%event_rate_m3_s(FUSION_DD_TP + 1)
    rate_he3n = rates%event_rate_m3_s(FUSION_DD_HE3N + 1)
    expected_tp = 0.5_c_double * density(FUSION_DEUTERON + 1)**2 * &
         rates%reactivity_m3_s(FUSION_DD_TP + 1)
    expected_he3n = 0.5_c_double * density(FUSION_DEUTERON + 1)**2 * &
         rates%reactivity_m3_s(FUSION_DD_HE3N + 1)
    call check(rate_tp > 0.0_c_double .and. rate_he3n > 0.0_c_double, &
         '10 keV DD branches have nonzero event rates')
    call check(close_relative(rate_tp, expected_tp, 1.0e-12_c_double), &
         'DD Tp event rate includes one half')
    call check(close_relative(rate_he3n, expected_he3n, 1.0e-12_c_double), &
         'DD He3n event rate includes one half')

    call fusion_channel_q(FUSION_DD_TP, q_tp, status)
    call check(status == PB11_STATUS_OK .and. q_tp > 0.0_c_double, &
         'DD Tp channel Q wrapper returns positive Q')
    call fusion_channel_q(FUSION_DD_HE3N, q_he3n, status)
    call check(status == PB11_STATUS_OK .and. q_he3n > 0.0_c_double, &
         'DD He3n channel Q wrapper returns positive Q')
    call check(close_relative(rates%nuclear_power_W_m3(FUSION_DD_TP + 1), &
         rate_tp * q_tp, 1.0e-12_c_double), &
         'DD Tp nuclear power uses channel Q')
    call check(close_relative(rates%nuclear_power_W_m3(FUSION_DD_HE3N + 1), &
         rate_he3n * q_he3n, 1.0e-12_c_double), &
         'DD He3n nuclear power uses channel Q')
    expected_total = rates%nuclear_power_W_m3(FUSION_DD_TP + 1) + &
         rates%nuclear_power_W_m3(FUSION_DD_HE3N + 1)
    call check(close_relative(rates%total_nuclear_power_W_m3, expected_total, &
         1.0e-12_c_double), 'DD total nuclear power is the channel sum')

    call check(close_relative(rates%particles%reactant_loss( &
         FUSION_DEUTERON + 1), 2.0_c_double * (rate_tp + rate_he3n), &
         1.0e-12_c_double), 'DD network consumes two deuterons per event')
    call check(close_relative(rates%particles%product_birth(FUSION_PROTON + 1), &
         rate_tp, 1.0e-12_c_double), 'DD Tp births one proton')
    call check(close_relative(rates%particles%product_birth(FUSION_TRITON + 1), &
         rate_tp, 1.0e-12_c_double), 'DD Tp births one triton')
    call check(close_relative(rates%particles%product_birth( &
         FUSION_HELIUM3 + 1), rate_he3n, 1.0e-12_c_double), &
         'DD He3n births one helium-3')
    call check(close_relative(rates%particles%neutron_birth, rate_he3n, &
         1.0e-12_c_double), 'DD He3n records one neutron per event')
    call check(rates%event_rate_m3_s(FUSION_PB11_3ALPHA + 1) == 0.0_c_double .and. &
         rates%event_rate_m3_s(FUSION_DT_ALPHAN + 1) == 0.0_c_double .and. &
         rates%event_rate_m3_s(FUSION_DHE3_ALPHAP + 1) == 0.0_c_double, &
         'unselected channels remain disabled')

    mask = -1_c_int
    call fusion_thermal_rates(kt_j, density, mask, PB11_METHOD_INTEGRAL, &
         rates, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'negative Fortran channel mask is rejected')
    call check(thermal_is_zero(rates), 'negative mask clears rates output')
    mask = 32_c_int
    call fusion_thermal_rates(kt_j, density, mask, PB11_METHOD_INTEGRAL, &
         rates, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'channel mask above 31 is rejected')
    call check(thermal_is_zero(rates), 'large mask clears rates output')
  end subroutine verify_thermal_rates

  subroutine verify_scalar_wrappers()
    real(c_double) :: kt_j, reactivity, sigma_m2, q_j, sigma_mb
    integer(c_int) :: status

    kt_j = 10.0_c_double * joules_per_keV
    call fusion_thermal_reactivity(FUSION_DD_TP, kt_j, PB11_METHOD_FAST, &
         reactivity, status)
    call check(status == PB11_STATUS_OK .and. reactivity > 0.0_c_double .and. &
         ieee_is_finite(reactivity), &
         'scalar thermal reactivity wrapper returns finite value')

    call fusion_channel_q(FUSION_DD_TP, q_j, status)
    call check(status == PB11_STATUS_OK .and. q_j > 0.0_c_double, &
         'scalar channel-Q wrapper returns finite positive value')

    call fusion_cross_section(FUSION_DT_ALPHAN, kt_j, sigma_m2, status)
    sigma_mb = sigma_m2 / joules_per_mb
    call check(status == PB11_STATUS_OK .and. sigma_m2 > 0.0_c_double .and. &
         ieee_is_finite(sigma_m2), &
         'scalar DT cross-section wrapper returns finite value')
    call check(sigma_mb > 27.0_c_double .and. sigma_mb < 27.2_c_double, &
         '10 keV DT cross section is about 27.02 millibarn')

    call fusion_channel_q(-1_c_int, q_j, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         q_j == 0.0_c_double, 'bad scalar Q channel clears output')
    call fusion_thermal_reactivity(FUSION_CHANNEL_COUNT, kt_j, &
         PB11_METHOD_INTEGRAL, reactivity, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         reactivity == 0.0_c_double, 'bad scalar rate channel clears output')
    call fusion_cross_section(FUSION_DT_ALPHAN, -kt_j, sigma_m2, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         sigma_m2 == 0.0_c_double, 'negative scalar energy clears output')
  end subroutine verify_scalar_wrappers

  subroutine verify_errors()
    real(c_double) :: rates(5)
    type(fusion_particle_sources_v1) :: source
    integer(c_int) :: status

    rates = 0.0_c_double
    rates(1) = ieee_value(0.0_c_double, ieee_quiet_nan)
    call fusion_particle_sources(rates, source, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'NaN rate returns invalid-argument status')
    call check(source_is_zero(source), 'NaN rate clears outputs')

    rates = 0.0_c_double
    rates(2) = ieee_value(0.0_c_double, ieee_positive_inf)
    call fusion_particle_sources(rates, source, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'infinite rate returns invalid-argument status')
    call check(source_is_zero(source), 'infinite rate clears outputs')

    rates = 0.0_c_double
    rates(4) = -1.0_c_double
    call fusion_particle_sources(rates, source, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE, &
         'negative rate returns out-of-range status')
    call check(source_is_zero(source), 'negative rate clears outputs')

    rates = 0.0_c_double
    rates(1) = huge(0.0_c_double)
    call fusion_particle_sources(rates, source, status)
    call check(status == PB11_STATUS_NUMERICAL_FAILURE, &
         'unrepresentable birth returns numerical-failure status')
    call check(source_is_zero(source), 'overflow clears outputs')
  end subroutine verify_errors

end program test_fusion_fortran
