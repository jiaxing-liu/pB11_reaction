program test_fusion_beam_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, &
       ieee_value
  use fusion_beam_fortran
  implicit none

  real(c_double), parameter :: joules_per_keV = 1.602176634e-16_c_double
  real(c_double), parameter :: deuteron_mass = 3.3435837724e-27_c_double
  real(c_double), parameter :: triton_mass = 5.0073567446e-27_c_double
  integer :: failures

  failures = 0
  call verify_domain_query()
  call verify_cold_target_rate()
  call verify_cold_outside_window()
  call verify_warm_domain_flag_and_identity()
  call verify_equal_temperature_pair()
  call verify_unequal_temperature_pair()
  call verify_errors()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran fusion-beam test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran fusion-beam tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

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

  logical function window_is_zero(out)
    type(fusion_beam_window_v1), intent(in) :: out

    window_is_zero = out%resolved_reactivity_m3_s == 0.0_c_double .and. &
         out%projectile_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%target_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%relative_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%cm_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%resolved_pair_probability == 0.0_c_double .and. &
         out%unresolved_pair_probability == 0.0_c_double .and. &
         out%unresolved_relative_speed_m_s == 0.0_c_double .and. &
         out%quadrature_error_m3_s == 0.0_c_double .and. &
         out%energy_identity_error_J_m3_s == 0.0_c_double .and. &
         out%domain_incomplete == 0_c_int
  end function window_is_zero

  subroutine verify_domain_query()
    real(c_double), parameter :: minimum_keV(5) = &
         [0.0_c_double, 0.5_c_double, 0.5_c_double, 0.5_c_double, &
          0.3_c_double]
    real(c_double), parameter :: maximum_keV(5) = &
         [9760.0_c_double, 5000.0_c_double, 4900.0_c_double, &
          4700.0_c_double, 4800.0_c_double]
    real(c_double) :: minimum, maximum
    integer(c_int) :: channel, status

    do channel = 0_c_int, FUSION_CHANNEL_COUNT - 1_c_int
       minimum = -1.0_c_double
       maximum = -1.0_c_double
       call fusion_cross_section_domain(channel, minimum, maximum, status)
       call check(status == PB11_STATUS_OK, 'beam domain query returns OK')
       call check(close_relative(minimum, minimum_keV(channel + 1) * &
            joules_per_keV, 1.0e-14_c_double) .and. &
            close_relative(maximum, maximum_keV(channel + 1) * &
            joules_per_keV, 1.0e-14_c_double), &
            'beam domain query exposes the channel SI window')
       call check(ieee_is_finite(minimum) .and. ieee_is_finite(maximum) .and. &
            minimum >= 0.0_c_double .and. maximum > minimum, &
            'beam domain endpoints are ordered and finite')
    end do
  end subroutine verify_domain_query

  subroutine verify_cold_target_rate()
    real(c_double), parameter :: expected_sigma_m2 = &
         27.02e-3_c_double * 1.0e-28_c_double
    real(c_double) :: relative_energy, projectile_energy, speed
    real(c_double) :: reduced_mass, expected_rate
    type(fusion_beam_window_v1) :: out
    integer(c_int) :: status

    relative_energy = 10.0_c_double * joules_per_keV
    projectile_energy = relative_energy * (deuteron_mass + triton_mass) / &
         triton_mass
    speed = sqrt(2.0_c_double * projectile_energy / deuteron_mass)
    expected_rate = speed * expected_sigma_m2
    reduced_mass = deuteron_mass / &
         (1.0_c_double + deuteron_mass / triton_mass)

    call fusion_beam_maxwellian_window(FUSION_DT_ALPHAN, deuteron_mass, &
         triton_mass, projectile_energy, 0.0_c_double, out, status)
    call check(status == PB11_STATUS_OK, 'cold-target beam call returns OK')
    call check(out%domain_incomplete == 0_c_int, &
         'cold in-window target is complete')
    call check(out%resolved_pair_probability == 1.0_c_double .and. &
         out%unresolved_pair_probability == 0.0_c_double, &
         'cold in-window probabilities are resolved')
    call check(close_relative(out%resolved_reactivity_m3_s, expected_rate, &
         3.0e-3_c_double), 'cold-target reactivity matches sigma times speed')
    call check(out%resolved_reactivity_m3_s > 0.0_c_double, &
         'cold-target reactivity is positive')
    call check(close_relative(out%projectile_energy_reactivity_J_m3_s, &
         projectile_energy * out%resolved_reactivity_m3_s, &
         1.0e-12_c_double), 'cold projectile energy moment is conditional')
    call check(close_relative(out%relative_energy_reactivity_J_m3_s, &
         reduced_mass * speed * speed / 2.0_c_double * &
         out%resolved_reactivity_m3_s, 1.0e-12_c_double), &
         'cold relative energy moment is consistent')
    call check(close_relative(out%cm_energy_reactivity_J_m3_s, &
         (projectile_energy - relative_energy) * &
         out%resolved_reactivity_m3_s, 1.0e-12_c_double), &
         'cold CM energy moment is consistent')
    call check(out%target_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%energy_identity_error_J_m3_s == 0.0_c_double, &
         'cold target has no target moment or identity residual')
  end subroutine verify_cold_target_rate

  subroutine verify_cold_outside_window()
    real(c_double) :: relative_energy, projectile_energy, speed
    type(fusion_beam_window_v1) :: out
    integer(c_int) :: status

    relative_energy = 6000.0_c_double * joules_per_keV
    projectile_energy = relative_energy * (deuteron_mass + triton_mass) / &
         triton_mass
    speed = sqrt(2.0_c_double * projectile_energy / deuteron_mass)

    call fusion_beam_maxwellian_window(FUSION_DT_ALPHAN, deuteron_mass, &
         triton_mass, projectile_energy, 0.0_c_double, out, status)
    call check(status == PB11_STATUS_OK, &
         'cold outside-window call returns unresolved OK')
    call check(out%domain_incomplete == 1_c_int .and. &
         out%resolved_pair_probability == 0.0_c_double .and. &
         out%unresolved_pair_probability == 1.0_c_double, &
         'cold outside-window pair is marked unresolved')
    call check(close_relative(out%unresolved_relative_speed_m_s, speed, &
         1.0e-12_c_double), 'cold unresolved speed is returned')
    call check(window_is_zero_except_unresolved(out), &
         'cold unresolved result has no resolved moments')
  end subroutine verify_cold_outside_window

  logical function window_is_zero_except_unresolved(out)
    type(fusion_beam_window_v1), intent(in) :: out

    window_is_zero_except_unresolved = &
         out%resolved_reactivity_m3_s == 0.0_c_double .and. &
         out%projectile_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%target_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%relative_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%cm_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%quadrature_error_m3_s == 0.0_c_double .and. &
         out%energy_identity_error_J_m3_s == 0.0_c_double
  end function window_is_zero_except_unresolved

  subroutine verify_warm_domain_flag_and_identity()
    real(c_double) :: relative_energy, projectile_energy
    real(c_double) :: lhs, rhs, probability_sum
    type(fusion_beam_window_v1) :: out
    integer(c_int) :: status

    relative_energy = 10.0_c_double * joules_per_keV
    projectile_energy = relative_energy * (deuteron_mass + triton_mass) / &
         triton_mass
    call fusion_beam_maxwellian_window(FUSION_DT_ALPHAN, deuteron_mass, &
         triton_mass, projectile_energy, 5.0_c_double * joules_per_keV, &
         out, status)
    call check(status == PB11_STATUS_OK, 'warm-target beam call returns OK')
    call check(out%domain_incomplete == 1_c_int, &
         'positive-temperature target reports incomplete domain')
    call check(out%resolved_reactivity_m3_s > 0.0_c_double .and. &
         out%resolved_pair_probability >= 0.0_c_double .and. &
         out%unresolved_pair_probability >= 0.0_c_double, &
         'warm resolved result is finite and nonnegative')
    probability_sum = out%resolved_pair_probability + &
         out%unresolved_pair_probability
    call check(close_relative(probability_sum, 1.0_c_double, 1.0e-9_c_double), &
         'warm resolved and unresolved probabilities close')
    lhs = out%projectile_energy_reactivity_J_m3_s + &
         out%target_energy_reactivity_J_m3_s
    rhs = out%relative_energy_reactivity_J_m3_s + &
         out%cm_energy_reactivity_J_m3_s
    call check(close_relative(lhs, rhs, 1.0e-9_c_double), &
         'warm energy moments satisfy projectile plus target identity')
    call check(abs(out%energy_identity_error_J_m3_s) <= &
         1.0e-9_c_double * max(1.0_c_double, abs(lhs) + abs(rhs)), &
         'warm reported energy identity residual is small')
  end subroutine verify_warm_domain_flag_and_identity

  subroutine verify_equal_temperature_pair()
    real(c_double) :: temperature, expected_cm
    type(fusion_beam_window_v1) :: out
    integer(c_int) :: status

    temperature = 10.0_c_double * joules_per_keV
    call fusion_thermal_pair_maxwellian_window(FUSION_DT_ALPHAN, &
         deuteron_mass, triton_mass, temperature, temperature, out, status)
    call check(status == PB11_STATUS_OK, &
         'equal-temperature pair call returns OK')
    call check(out%domain_incomplete == 1_c_int .and. &
         out%resolved_reactivity_m3_s > 0.0_c_double, &
         'equal-temperature pair retains domain flag and rate')
    expected_cm = 1.5_c_double * temperature * &
         out%resolved_reactivity_m3_s
    call check(close_relative(out%cm_energy_reactivity_J_m3_s, expected_cm, &
         1.0e-12_c_double), &
         'equal-temperature CM energy is 3/2 kT times the rate')
    call check(close_relative(out%projectile_energy_reactivity_J_m3_s + &
         out%target_energy_reactivity_J_m3_s, &
         out%relative_energy_reactivity_J_m3_s + &
         out%cm_energy_reactivity_J_m3_s, 1.0e-9_c_double), &
         'equal-temperature pair satisfies energy identity')
  end subroutine verify_equal_temperature_pair

  subroutine verify_unequal_temperature_pair()
    real(c_double) :: temperature_a, temperature_b, lhs, rhs
    type(fusion_beam_window_v1) :: out
    integer(c_int) :: status

    temperature_a = 4.0_c_double * joules_per_keV
    temperature_b = 12.0_c_double * joules_per_keV
    call fusion_thermal_pair_maxwellian_window(FUSION_DT_ALPHAN, &
         deuteron_mass, triton_mass, temperature_a, temperature_b, out, &
         status)
    call check(status == PB11_STATUS_OK, &
         'unequal-temperature pair call returns OK')
    call check(out%domain_incomplete == 1_c_int .and. &
         out%resolved_reactivity_m3_s > 0.0_c_double, &
         'unequal-temperature pair retains domain flag and rate')
    lhs = out%projectile_energy_reactivity_J_m3_s + &
         out%target_energy_reactivity_J_m3_s
    rhs = out%relative_energy_reactivity_J_m3_s + &
         out%cm_energy_reactivity_J_m3_s
    call check(close_relative(lhs, rhs, 1.0e-9_c_double), &
         'unequal-temperature pair satisfies energy identity')
    call check(abs(out%energy_identity_error_J_m3_s) <= &
         1.0e-9_c_double * max(1.0_c_double, abs(lhs) + abs(rhs)), &
         'unequal-temperature pair reports small identity residual')
  end subroutine verify_unequal_temperature_pair

  subroutine verify_errors()
    type(fusion_beam_window_v1) :: out
    real(c_double) :: minimum, maximum
    real(c_double) :: nan_value
    integer(c_int) :: status

    out = fusion_beam_window_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1_c_int)
    call fusion_beam_maxwellian_window(-1_c_int, deuteron_mass, triton_mass, &
         10.0_c_double * joules_per_keV, 0.0_c_double, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'invalid beam channel returns invalid-argument status')
    call check(window_is_zero(out), 'invalid beam channel clears output')

    out = fusion_beam_window_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1_c_int)
    call fusion_beam_maxwellian_window(FUSION_DT_ALPHAN, -deuteron_mass, &
         triton_mass, 10.0_c_double * joules_per_keV, 0.0_c_double, out, &
         status)
    call check(status == PB11_STATUS_OUT_OF_RANGE, &
         'invalid beam mass returns out-of-range status')
    call check(window_is_zero(out), 'invalid beam mass clears output')

    nan_value = ieee_value(0.0_c_double, ieee_quiet_nan)
    out = fusion_beam_window_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1_c_int)
    call fusion_beam_maxwellian_window(FUSION_DT_ALPHAN, deuteron_mass, &
         triton_mass, 10.0_c_double * joules_per_keV, nan_value, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'NaN target temperature returns invalid-argument status')
    call check(window_is_zero(out), 'NaN target temperature clears output')

    out = fusion_beam_window_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1_c_int)
    call fusion_thermal_pair_maxwellian_window(FUSION_DT_ALPHAN, &
         deuteron_mass, triton_mass, -1.0_c_double * joules_per_keV, &
         10.0_c_double * joules_per_keV, out, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE, &
         'negative pair temperature returns out-of-range status')
    call check(window_is_zero(out), 'negative pair temperature clears output')

    minimum = 2.0_c_double
    maximum = 3.0_c_double
    call fusion_cross_section_domain(-1_c_int, minimum, maximum, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         minimum == 0.0_c_double .and. maximum == 0.0_c_double, &
         'invalid domain query clears both outputs')
  end subroutine verify_errors

end program test_fusion_beam_fortran
