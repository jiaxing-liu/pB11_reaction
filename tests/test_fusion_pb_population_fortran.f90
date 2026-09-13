program test_fusion_pb_population_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_pb_population_fortran
  use fusion_beam_fortran, only : fusion_beam_window_v1
  use fusion_nuclear_data_fortran, only : fusion_nuclear_mass_v1, &
       fusion_nuclear_mass, FUSION_PROTON, FUSION_BORON11
  implicit none

  real(c_double), parameter :: joules_per_keV = 1.602176634e-16_c_double
  real(c_double), parameter :: rel_tol = 5.0e-8_c_double
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_population_point()
  call verify_invalid_population_clears()
  call verify_thermal_population_sums()
  call verify_beam_population_invalid_mass_clears()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran pB-population test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran pB-population tests passed'

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

  logical function population_is_zero(out)
    type(fusion_pb_population_v1), intent(in) :: out

    population_is_zero = out%total_cross_section_m2 == 0.0_c_double .and. &
         out%narrow_fit_cross_section_m2 == 0.0_c_double .and. &
         out%narrow_fit_fraction == 0.0_c_double .and. &
         out%alpha0_peak_fraction == 0.0_c_double .and. &
         out%narrow_remainder_fraction == 0.0_c_double .and. &
         out%other_remainder_fraction == 0.0_c_double .and. &
         out%continuum_peak_fraction == 0.0_c_double .and. &
         out%continuum_extrapolated_fraction == 0.0_c_double .and. &
         out%equivalent_proton_lab_energy_J == 0.0_c_double
  end function population_is_zero

  subroutine clear_window_for_test(out)
    type(fusion_beam_window_v1), intent(out) :: out

    out%resolved_reactivity_m3_s = -1.0_c_double
    out%projectile_energy_reactivity_J_m3_s = -1.0_c_double
    out%target_energy_reactivity_J_m3_s = -1.0_c_double
    out%relative_energy_reactivity_J_m3_s = -1.0_c_double
    out%cm_energy_reactivity_J_m3_s = -1.0_c_double
    out%resolved_pair_probability = -1.0_c_double
    out%unresolved_pair_probability = -1.0_c_double
    out%unresolved_relative_speed_m_s = -1.0_c_double
    out%quadrature_error_m3_s = -1.0_c_double
    out%energy_identity_error_J_m3_s = -1.0_c_double
    out%domain_incomplete = -1_c_int
  end subroutine clear_window_for_test

  subroutine fill_rates_for_test(out)
    type(fusion_pb_population_rates_v1), intent(out) :: out

    call clear_window_for_test(out%total)
    call clear_window_for_test(out%alpha0_peak)
    call clear_window_for_test(out%narrow_remainder)
    call clear_window_for_test(out%other_remainder)
    call clear_window_for_test(out%continuum_extrapolated)
  end subroutine fill_rates_for_test

  logical function rates_are_zero(out)
    type(fusion_pb_population_rates_v1), intent(in) :: out

    rates_are_zero = window_is_zero(out%total) .and. &
         window_is_zero(out%alpha0_peak) .and. &
         window_is_zero(out%narrow_remainder) .and. &
         window_is_zero(out%other_remainder) .and. &
         window_is_zero(out%continuum_extrapolated)
  end function rates_are_zero

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

  logical function window_is_finite(out)
    type(fusion_beam_window_v1), intent(in) :: out

    window_is_finite = ieee_is_finite(out%resolved_reactivity_m3_s) .and. &
         ieee_is_finite(out%projectile_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%target_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%relative_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%cm_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%resolved_pair_probability) .and. &
         ieee_is_finite(out%unresolved_pair_probability) .and. &
         ieee_is_finite(out%unresolved_relative_speed_m_s) .and. &
         ieee_is_finite(out%quadrature_error_m3_s) .and. &
         ieee_is_finite(out%energy_identity_error_J_m3_s)
  end function window_is_finite

  subroutine verify_layout()
    type(fusion_pb_population_v1) :: point
    type(fusion_pb_population_rates_v1) :: rates
    type(fusion_beam_window_v1) :: window
    real(c_double) :: scalar

    call check(c_sizeof(point) == 9_c_size_t * c_sizeof(scalar), &
         'population point layout is nine doubles')
    call check(c_sizeof(rates) == 5_c_size_t * c_sizeof(window), &
         'population rates layout is five beam windows')
  end subroutine verify_layout

  subroutine verify_population_point()
    type(fusion_pb_population_v1) :: out
    type(fusion_nuclear_mass_v1) :: proton, boron
    real(c_double), parameter :: energy = 148.0_c_double * joules_per_keV
    real(c_double) :: partition, expected_lab_energy
    integer(c_int) :: status, mass_status

    call fusion_nuclear_mass(FUSION_PROTON, proton, mass_status)
    call check(mass_status == PB11_STATUS_OK, &
         'canonical proton mass query succeeds')
    call fusion_nuclear_mass(FUSION_BORON11, boron, mass_status)
    call check(mass_status == PB11_STATUS_OK, &
         'canonical boron mass query succeeds')

    call fusion_pb_population(FUSION_ENDPOINT_S, FUSION_PB_LOW_TB, energy, &
         0.5_c_double, 1.0_c_double, out, status)
    call check(status == PB11_STATUS_OK, '148 keV population call returns OK')
    call check(out%total_cross_section_m2 >= 0.0_c_double .and. &
         out%narrow_fit_cross_section_m2 >= 0.0_c_double .and. &
         out%narrow_fit_fraction >= 0.0_c_double .and. &
         out%narrow_fit_fraction <= 1.0_c_double .and. &
         out%alpha0_peak_fraction >= 0.0_c_double .and. &
         out%alpha0_peak_fraction <= 1.0_c_double .and. &
         out%narrow_remainder_fraction >= 0.0_c_double .and. &
         out%other_remainder_fraction >= 0.0_c_double .and. &
         out%continuum_peak_fraction >= 0.0_c_double .and. &
         out%continuum_extrapolated_fraction >= 0.0_c_double .and. &
         ieee_is_finite(out%equivalent_proton_lab_energy_J), &
         '148 keV population fields are finite and nonnegative')
    partition = out%alpha0_peak_fraction + &
         out%narrow_remainder_fraction + out%other_remainder_fraction
    call check(close_relative(partition, 1.0_c_double, rel_tol), &
         'three population fractions sum to one')
    expected_lab_energy = energy * (proton%mass_kg + boron%mass_kg) / &
         boron%mass_kg
    call check(close_relative(out%equivalent_proton_lab_energy_J, &
         expected_lab_energy, 1.0e-14_c_double), &
         'equivalent proton lab energy uses canonical masses')
    call check(close_relative(out%narrow_fit_cross_section_m2, &
         out%narrow_fit_fraction * out%total_cross_section_m2, 1.0e-14_c_double), &
         'narrow fit cross section follows its fraction')
  end subroutine verify_population_point

  subroutine verify_invalid_population_clears()
    type(fusion_pb_population_v1) :: out
    integer(c_int) :: status

    out%total_cross_section_m2 = -1.0_c_double
    out%narrow_fit_cross_section_m2 = -1.0_c_double
    out%narrow_fit_fraction = -1.0_c_double
    out%alpha0_peak_fraction = -1.0_c_double
    out%narrow_remainder_fraction = -1.0_c_double
    out%other_remainder_fraction = -1.0_c_double
    out%continuum_peak_fraction = -1.0_c_double
    out%continuum_extrapolated_fraction = -1.0_c_double
    out%equivalent_proton_lab_energy_J = -1.0_c_double

    call fusion_pb_population(0_c_int, FUSION_PB_LOW_TB, &
         148.0_c_double * joules_per_keV, 0.5_c_double, 1.0_c_double, out, &
         status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'invalid population policy returns invalid-argument status')
    call check(population_is_zero(out), &
         'invalid population policy clears every point field')
  end subroutine verify_invalid_population_clears

  subroutine verify_thermal_population_sums()
    type(fusion_nuclear_mass_v1) :: proton, boron
    type(fusion_pb_population_rates_v1) :: out
    real(c_double), parameter :: temperature = 100.0_c_double * &
         joules_per_keV
    integer(c_int) :: status, mass_status

    call fusion_nuclear_mass(FUSION_PROTON, proton, mass_status)
    call fusion_nuclear_mass(FUSION_BORON11, boron, mass_status)
    call check(mass_status == PB11_STATUS_OK, &
         'canonical masses are available for thermal population rates')
    call fusion_pb_thermal_population_rates(FUSION_ENDPOINT_S, &
         FUSION_PB_LOW_TB, proton%mass_kg, boron%mass_kg, temperature, &
         temperature, 0.5_c_double, 1.0_c_double, out, status)
    call check(status == PB11_STATUS_OK, &
         'thermal population rates call returns OK')
    call check(window_is_finite(out%total) .and. &
         out%total%resolved_reactivity_m3_s > 0.0_c_double .and. &
         window_is_finite(out%alpha0_peak) .and. &
         window_is_finite(out%narrow_remainder) .and. &
         window_is_finite(out%other_remainder) .and. &
         window_is_finite(out%continuum_extrapolated), &
         'thermal population rate windows are finite')
    call check_rate_sums(out)
  end subroutine verify_thermal_population_sums

  subroutine check_rate_sums(out)
    type(fusion_pb_population_rates_v1), intent(in) :: out
    real(c_double) :: expected

    expected = out%alpha0_peak%resolved_reactivity_m3_s + &
         out%narrow_remainder%resolved_reactivity_m3_s + &
         out%other_remainder%resolved_reactivity_m3_s
    call check(close_relative(out%total%resolved_reactivity_m3_s, expected, &
         rel_tol), 'three population event rates sum to total')

    expected = out%alpha0_peak%projectile_energy_reactivity_J_m3_s + &
         out%narrow_remainder%projectile_energy_reactivity_J_m3_s + &
         out%other_remainder%projectile_energy_reactivity_J_m3_s
    call check(close_relative(out%total%projectile_energy_reactivity_J_m3_s, &
         expected, rel_tol), 'projectile energy moments sum to total')

    expected = out%alpha0_peak%target_energy_reactivity_J_m3_s + &
         out%narrow_remainder%target_energy_reactivity_J_m3_s + &
         out%other_remainder%target_energy_reactivity_J_m3_s
    call check(close_relative(out%total%target_energy_reactivity_J_m3_s, &
         expected, rel_tol), 'target energy moments sum to total')

    expected = out%alpha0_peak%relative_energy_reactivity_J_m3_s + &
         out%narrow_remainder%relative_energy_reactivity_J_m3_s + &
         out%other_remainder%relative_energy_reactivity_J_m3_s
    call check(close_relative(out%total%relative_energy_reactivity_J_m3_s, &
         expected, rel_tol), 'relative energy moments sum to total')

    expected = out%alpha0_peak%cm_energy_reactivity_J_m3_s + &
         out%narrow_remainder%cm_energy_reactivity_J_m3_s + &
         out%other_remainder%cm_energy_reactivity_J_m3_s
    call check(close_relative(out%total%cm_energy_reactivity_J_m3_s, expected, &
         rel_tol), 'CM energy moments sum to total')
  end subroutine check_rate_sums

  subroutine verify_beam_population_invalid_mass_clears()
    type(fusion_nuclear_mass_v1) :: proton, boron
    type(fusion_pb_population_rates_v1) :: out
    real(c_double), parameter :: target_temperature = 100.0_c_double * &
         joules_per_keV
    integer(c_int) :: status, mass_status

    call fusion_nuclear_mass(FUSION_PROTON, proton, mass_status)
    call fusion_nuclear_mass(FUSION_BORON11, boron, mass_status)
    call check(mass_status == PB11_STATUS_OK, &
         'canonical masses are available for beam population error test')
    call fill_rates_for_test(out)
    call fusion_pb_beam_population_rates(FUSION_ENDPOINT_S, FUSION_PB_LOW_TB, &
         proton%mass_kg * 1.001_c_double, boron%mass_kg, &
         148.0_c_double * joules_per_keV, target_temperature, 0.5_c_double, &
         1.0_c_double, out, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE, &
         'noncanonical beam mass returns out-of-range status')
    call check(rates_are_zero(out), &
         'invalid beam mass clears every population-rate window')
  end subroutine verify_beam_population_invalid_mass_clears

end program test_fusion_pb_population_fortran
