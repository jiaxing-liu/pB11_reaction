program test_fusion_two_component_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_kinetics_fortran, only : fusion_maxwellian_bath_v1, &
       PB11_STATUS_OK, PB11_STATUS_OUT_OF_RANGE, PB11_STATUS_INVALID_ARGUMENT
  use fusion_two_component_fortran
  implicit none

  real(c_double), parameter :: joules_per_keV = 1.602176634e-16_c_double
  integer :: failures

  failures = 0
  call verify_coulomb_transfer_rate()
  call verify_one_cell_no_bath()
  call verify_invalid_extent()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, &
          ' Fortran two-component binding test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran two-component binding tests passed'

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

  logical function ledger_is_zero(ledger)
    use fusion_kinetics_fortran, only : fusion_kinetic_ledger_v1
    type(fusion_two_component_ledger_v1), intent(in) :: ledger
    type(fusion_kinetic_ledger_v1) :: total

    total = ledger%total
    ledger_is_zero = total%initial_number_m3 == 0.0_c_double .and. &
         total%final_number_m3 == 0.0_c_double .and. &
         total%initial_energy_J_m3 == 0.0_c_double .and. &
         total%final_energy_J_m3 == 0.0_c_double .and. &
         total%born_number_m3 == 0.0_c_double .and. &
         total%born_energy_J_m3 == 0.0_c_double .and. &
         total%escaped_number_m3 == 0.0_c_double .and. &
         total%escaped_energy_J_m3 == 0.0_c_double .and. &
         total%thermalized_number_m3 == 0.0_c_double .and. &
         total%thermalized_energy_J_m3 == 0.0_c_double .and. &
         total%particle_balance_error_m3 == 0.0_c_double .and. &
         total%energy_balance_error_J_m3 == 0.0_c_double .and. &
         ledger%transferred_number_m3 == 0.0_c_double .and. &
         ledger%transferred_energy_J_m3 == 0.0_c_double
  end function ledger_is_zero

  subroutine verify_coulomb_transfer_rate()
    real(c_double), parameter :: elementary_charge = 1.602176634e-19_c_double
    real(c_double), parameter :: epsilon0 = 8.8541878188e-12_c_double
    real(c_double), parameter :: pi = 3.1415926535897932384626433832795_c_double
    real(c_double), parameter :: proton_mass = 1.67262192369e-27_c_double
    real(c_double), parameter :: deuteron_mass = 3.3435837724e-27_c_double
    real(c_double) :: energy, velocity, thermal_velocity, coulomb
    real(c_double) :: coefficient, expected, rate
    type(fusion_maxwellian_bath_v1) :: bath
    integer(c_int) :: status

    bath%density_m3 = 2.5e19_c_double
    bath%mass_kg = deuteron_mass
    bath%mean_charge_squared = 1.75_c_double
    bath%kT_J = 8.0_c_double * joules_per_keV
    bath%coulomb_log = 12.3_c_double
    energy = 15.0_c_double * joules_per_keV

    call fusion_coulomb_transfer_rate(energy, proton_mass, 1.0_c_double, &
         bath, rate, status)
    velocity = sqrt(2.0_c_double * energy / proton_mass)
    thermal_velocity = sqrt(2.0_c_double * bath%kT_J / bath%mass_kg)
    coulomb = elementary_charge**2 / (4.0_c_double * pi * epsilon0)
    coefficient = 4.0_c_double * pi * bath%density_m3 * &
         bath%mean_charge_squared * coulomb**2 * bath%coulomb_log / &
         (proton_mass * bath%mass_kg)
    expected = 4.0_c_double * coefficient / &
         (sqrt(pi) * thermal_velocity**3) * &
         exp(-(velocity / thermal_velocity)**2)
    call check(status == PB11_STATUS_OK, &
         'Coulomb transfer coefficient returns OK')
    call check(close_relative(rate, expected, 1.0e-12_c_double), &
         'Coulomb transfer coefficient matches independent formula')
    call check(rate > 0.0_c_double .and. ieee_is_finite(rate), &
         'Coulomb transfer coefficient is finite and positive')

    bath%kT_J = -1.0_c_double
    rate = -1.0_c_double
    call fusion_coulomb_transfer_rate(energy, proton_mass, 1.0_c_double, &
         bath, rate, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE, &
         'negative ion-bath temperature is rejected')
    call check(rate == 0.0_c_double .and. ieee_is_finite(rate), &
         'invalid transfer coefficient clears its output')
  end subroutine verify_coulomb_transfer_rate

  subroutine verify_one_cell_no_bath()
    real(c_double) :: edges(2), old_s(1), old_t(1)
    real(c_double) :: bath(0), diffusion(0), birth_s(1), birth_t(1)
    real(c_double) :: escape(1), transfer(1), trial_s(1), trial_t(1)
    real(c_double) :: trial_s_again(1), trial_t_again(1), heat(0)
    type(fusion_two_component_ledger_v1) :: ledger, ledger_again
    real(c_double) :: expected_s, expected_t, expected_born
    real(c_double) :: expected_escape, cell_energy, transfer_number
    real(c_double) :: initial_number, final_number
    integer(c_int) :: status, status_again
    real(c_double), parameter :: dt = 0.5_c_double

    edges = [1.0e-16_c_double, 3.0e-16_c_double]
    old_s = [4.0e20_c_double]
    old_t = [1.0e20_c_double]
    birth_s = [3.0e20_c_double]
    birth_t = [5.0e19_c_double]
    escape = [0.2_c_double]
    transfer = [0.4_c_double]
    trial_s = -1.0_c_double
    trial_t = -1.0_c_double
    heat = 0.0_c_double

    call fusion_two_component_trial(1_c_int, 0_c_int, dt, edges, old_s, &
         old_t, bath, diffusion, birth_s, birth_t, escape, transfer, &
         trial_s, trial_t, heat, ledger, status)
    call check(status == PB11_STATUS_OK, &
         'one-cell no-bath two-component trial returns OK')

    expected_s = (old_s(1) + dt * birth_s(1)) / &
         (1.0_c_double + dt * (escape(1) + transfer(1)))
    expected_t = (old_t(1) + dt * birth_t(1) + &
         dt * transfer(1) * expected_s) / &
         (1.0_c_double + dt * escape(1))
    call check(close_relative(trial_s(1), expected_s, 1.0e-14_c_double), &
         'one-cell S component matches analytic implicit solution')
    call check(close_relative(trial_t(1), expected_t, 1.0e-14_c_double), &
         'one-cell T component matches analytic implicit solution')

    cell_energy = 2.0e-16_c_double
    initial_number = old_s(1) + old_t(1)
    final_number = expected_s + expected_t
    expected_born = dt * (birth_s(1) + birth_t(1))
    expected_escape = dt * escape(1) * final_number
    transfer_number = dt * transfer(1) * expected_s
    call check(close_relative(ledger%total%initial_number_m3, &
         initial_number, 1.0e-14_c_double), &
         'nested ledger records both initial components')
    call check(close_relative(ledger%total%final_number_m3, final_number, &
         1.0e-14_c_double), 'nested ledger records both final components')
    call check(close_relative(ledger%total%born_number_m3, expected_born, &
         1.0e-14_c_double), 'nested ledger records external births once')
    call check(close_relative(ledger%total%escaped_number_m3, &
         expected_escape, 1.0e-14_c_double), &
         'nested ledger records physical escape once')
    call check(ledger%total%thermalized_number_m3 == 0.0_c_double .and. &
         ledger%total%thermalized_energy_J_m3 == 0.0_c_double, &
         'internal transfer is absent from the thermalized ledger fields')
    call check(close_relative(ledger%transferred_number_m3, transfer_number, &
         1.0e-14_c_double), 'internal transfer number is reported separately')
    call check(close_relative(ledger%transferred_energy_J_m3, &
         cell_energy * transfer_number, 1.0e-14_c_double), &
         'internal carried transfer energy is reported separately')
    call check(close_relative(ledger%total%initial_energy_J_m3, &
         cell_energy * initial_number, 1.0e-14_c_double), &
         'nested ledger records initial carried energy')
    call check(close_relative(ledger%total%born_energy_J_m3, &
         cell_energy * expected_born, 1.0e-14_c_double), &
         'nested ledger records birth carried energy')
    call check(close_relative(ledger%total%escaped_energy_J_m3, &
         cell_energy * expected_escape, 1.0e-14_c_double), &
         'nested ledger records physical escaped energy')

    ! An identical second trial must not consume either accepted input pool.
    call fusion_two_component_trial(1_c_int, 0_c_int, dt, edges, old_s, &
         old_t, bath, diffusion, birth_s, birth_t, escape, transfer, &
         trial_s_again, trial_t_again, heat, ledger_again, status_again)
    call check(status_again == PB11_STATUS_OK, 'repeated trial returns OK')
    call check(trial_s_again(1) == trial_s(1) .and. &
         trial_t_again(1) == trial_t(1), &
         'repeated trial is deterministic')
    call check(old_s(1) == 4.0e20_c_double .and. &
         old_t(1) == 1.0e20_c_double, &
         'repeated trial does not modify old component inputs')
    call check(ledger_again%total%thermalized_number_m3 == 0.0_c_double, &
         'repeated trial does not accumulate internal ash state')
  end subroutine verify_one_cell_no_bath

  subroutine verify_invalid_extent()
    real(c_double) :: edges(2), old_s(2), old_t(2), bath(0), diffusion(0)
    real(c_double) :: birth_s(2), birth_t(2), escape(2), transfer(2)
    real(c_double) :: trial_s(2), trial_t(2), heat(0)
    real(c_double) :: one_edge(2), one_old_s(1), one_old_t(1)
    real(c_double) :: one_birth_s(1), one_birth_t(1), one_escape(1)
    real(c_double) :: one_transfer(1), one_trial_s(1), one_trial_t(1)
    real(c_double) :: one_bath(0), one_diffusion(0), one_heat(1)
    type(fusion_two_component_ledger_v1) :: ledger
    integer(c_int) :: status

    edges = [1.0e-16_c_double, 2.0e-16_c_double]
    old_s = 1.0e20_c_double
    old_t = 2.0e20_c_double
    birth_s = 0.0_c_double
    birth_t = 0.0_c_double
    escape = 0.0_c_double
    transfer = 0.0_c_double
    trial_s = -1.0_c_double
    trial_t = -1.0_c_double
    call fusion_two_component_trial(2_c_int, 0_c_int, 0.1_c_double, edges, &
         old_s, old_t, bath, diffusion, birth_s, birth_t, escape, transfer, &
         trial_s, trial_t, heat, ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong edge extent is rejected before C_LOC')
    call check(all(trial_s == 0.0_c_double) .and. &
         all(trial_t == 0.0_c_double) .and. ledger_is_zero(ledger), &
         'wrong edge extent clears trial arrays and nested ledger')

    one_edge = [1.0e-16_c_double, 2.0e-16_c_double]
    one_old_s = 1.0e20_c_double
    one_old_t = 1.0e20_c_double
    one_birth_s = 0.0_c_double
    one_birth_t = 0.0_c_double
    one_escape = 0.0_c_double
    one_transfer = 0.0_c_double
    one_trial_s = -1.0_c_double
    one_trial_t = -1.0_c_double
    one_heat = -1.0_c_double
    call fusion_two_component_trial(1_c_int, 1_c_int, 0.1_c_double, one_edge, &
         one_old_s, one_old_t, one_bath, one_diffusion, one_birth_s, &
         one_birth_t, one_escape, one_transfer, one_trial_s, one_trial_t, &
         one_heat, ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'missing positive-bath extent is rejected before C_NULL_PTR misuse')
    call check(one_trial_s(1) == 0.0_c_double .and. &
         one_trial_t(1) == 0.0_c_double .and. one_heat(1) == 0.0_c_double &
         .and. ledger_is_zero(ledger), &
         'invalid bath extent clears all nonempty outputs')
  end subroutine verify_invalid_extent

end program test_fusion_two_component_fortran
