program test_fusion_kinetics_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, &
       ieee_value
  use fusion_kinetics_fortran
  implicit none

  real(c_double), parameter :: joules_per_keV = 1.602176634e-16_c_double
  real(c_double), parameter :: diffusion = 1.0e-30_c_double
  integer :: failures

  failures = 0
  call verify_coulomb_reference()
  call verify_one_cell()
  call verify_two_bath_exchange()
  call verify_noncontiguous_sections()
  call verify_input_errors()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran fusion-kinetics test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran fusion-kinetics tests passed'

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
    type(fusion_kinetic_ledger_v1), intent(in) :: ledger

    ledger_is_zero = ledger%initial_number_m3 == 0.0_c_double .and. &
         ledger%final_number_m3 == 0.0_c_double .and. &
         ledger%initial_energy_J_m3 == 0.0_c_double .and. &
         ledger%final_energy_J_m3 == 0.0_c_double .and. &
         ledger%born_number_m3 == 0.0_c_double .and. &
         ledger%born_energy_J_m3 == 0.0_c_double .and. &
         ledger%escaped_number_m3 == 0.0_c_double .and. &
         ledger%escaped_energy_J_m3 == 0.0_c_double .and. &
         ledger%thermalized_number_m3 == 0.0_c_double .and. &
         ledger%thermalized_energy_J_m3 == 0.0_c_double .and. &
         ledger%particle_balance_error_m3 == 0.0_c_double .and. &
         ledger%energy_balance_error_J_m3 == 0.0_c_double
  end function ledger_is_zero

  logical function coulomb_is_zero(out)
    type(fusion_coulomb_energy_v1), intent(in) :: out

    coulomb_is_zero = out%diffusion_J2_s == 0.0_c_double .and. &
         out%mean_energy_rate_J_s == 0.0_c_double
  end function coulomb_is_zero

  subroutine verify_coulomb_reference()
    real(c_double), parameter :: expected_diffusion = 7.67006381506318e-29_c_double
    real(c_double), parameter :: expected_drift = -5.33071110657554e-14_c_double
    type(fusion_maxwellian_bath_v1) :: bath
    type(fusion_coulomb_energy_v1) :: out
    integer(c_int) :: status

    ! Independent double-precision reference for the NRL/Trubnikov formula
    ! at 15 keV proton energy in an 8 keV deuterium bath.  The reference was
    ! evaluated separately from the wrapper and is intentionally hard-coded.
    bath%density_m3 = 2.5e19_c_double
    bath%mass_kg = 3.3435837724e-27_c_double
    bath%mean_charge_squared = 1.75_c_double
    bath%kT_J = 8.0_c_double * joules_per_keV
    bath%coulomb_log = 12.3_c_double

    call fusion_coulomb_energy(15.0_c_double * joules_per_keV, &
         1.67262192369e-27_c_double, -1.0_c_double, bath, out, status)
    call check(status == PB11_STATUS_OK, 'Coulomb reference call returns OK')
    call check(close_relative(out%diffusion_J2_s, expected_diffusion, &
         1.0e-11_c_double), 'Coulomb diffusion matches independent reference')
    call check(close_relative(out%mean_energy_rate_J_s, expected_drift, &
         1.0e-11_c_double), 'Coulomb drift matches independent reference')

    call fusion_coulomb_energy(-1.0_c_double, 1.0e-27_c_double, 1.0_c_double, &
         bath, out, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE, &
         'negative Coulomb energy is rejected')
    call check(coulomb_is_zero(out), 'Coulomb error clears output')
  end subroutine verify_coulomb_reference

  subroutine verify_one_cell()
    real(c_double) :: edges(2), old(1), birth(1), escape(1), trial(1)
    real(c_double) :: bath(0), no_diffusion(0), heat(0)
    type(fusion_kinetic_ledger_v1) :: ledger
    real(c_double) :: expected, expected_born, expected_escape
    real(c_double) :: expected_thermalized
    integer(c_int) :: status

    edges = [1.0e-16_c_double, 3.0e-16_c_double]
    old = [4.0e20_c_double]
    birth = [3.0e20_c_double]
    escape = [0.2_c_double]
    trial = -1.0_c_double
    ledger = fusion_kinetic_ledger_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double)

    call fusion_energy_fp_trial(1_c_int, 0_c_int, 0.5_c_double, edges, old, &
         bath, no_diffusion, birth, escape, 0.4_c_double, trial, heat, &
         ledger, status)
    call check(status == PB11_STATUS_OK, 'one-cell Fortran trial returns OK')
    expected = (old(1) + 0.5_c_double * birth(1)) / &
         (1.0_c_double + 0.5_c_double * (escape(1) + 0.4_c_double))
    call check(close_relative(trial(1), expected, 1.0e-14_c_double), &
         'one-cell birth/loss matches implicit analytic result')
    call check(old(1) == 4.0e20_c_double, &
         'one-cell trial leaves old input unchanged')
    call check(ledger%initial_number_m3 == old(1), &
         'one-cell ledger records initial number')
    call check(close_relative(ledger%final_number_m3, expected, &
         1.0e-14_c_double), 'one-cell ledger records final number')

    expected_born = 0.5_c_double * birth(1)
    expected_escape = 0.5_c_double * escape(1) * expected
    expected_thermalized = 0.5_c_double * 0.4_c_double * expected
    call check(close_relative(ledger%born_number_m3, expected_born, &
         1.0e-14_c_double), 'one-cell ledger records birth once')
    call check(close_relative(ledger%escaped_number_m3, expected_escape, &
         1.0e-14_c_double), 'one-cell ledger records escape once')
    call check(close_relative(ledger%thermalized_number_m3, &
         expected_thermalized, 1.0e-14_c_double), &
         'one-cell ledger records thermalization once')
    call check(abs(ledger%particle_balance_error_m3) < 1.0e10_c_double, &
         'one-cell particle ledger residual is small')
    call check(abs(ledger%energy_balance_error_J_m3) < 1.0e-5_c_double, &
         'one-cell energy ledger residual is small')
  end subroutine verify_one_cell

  subroutine verify_two_bath_exchange()
    real(c_double) :: edges(4), old(3), bath(2), diffusion_values(4)
    real(c_double) :: birth(3), escape(3), trial(3), heat(2)
    type(fusion_kinetic_ledger_v1) :: ledger
    real(c_double) :: expected_energy
    integer(c_int) :: status

    edges = [1.0e-16_c_double, 2.0e-16_c_double, 4.0e-16_c_double, &
         8.0e-16_c_double]
    old = [0.0_c_double, 1.0e20_c_double, 0.0_c_double]
    bath = [1.0e-16_c_double, 8.0e-16_c_double]
    ! C order is bath-major: [bath 0 faces, bath 1 faces].
    diffusion_values = [diffusion, diffusion, diffusion, diffusion]
    birth = 0.0_c_double
    escape = 0.0_c_double
    trial = -1.0_c_double
    heat = -1.0_c_double

    call fusion_energy_fp_trial(3_c_int, 2_c_int, 0.05_c_double, edges, old, &
         bath, diffusion_values, birth, escape, 0.0_c_double, trial, heat, &
         ledger, status)
    call check(status == PB11_STATUS_OK, 'two-bath Fortran trial returns OK')
    call check(all(ieee_is_finite(trial)) .and. all(trial >= 0.0_c_double), &
         'two-bath trial populations are finite and nonnegative')
    call check(heat(1) > 0.0_c_double, &
         'colder bath receives heat from the fast population')
    call check(heat(2) < 0.0_c_double, &
         'hotter bath supplies heat to the fast population')
    call check(close_relative(ledger%final_number_m3, &
         ledger%initial_number_m3, 1.0e-12_c_double), &
         'two-bath reflecting step conserves particles')
    expected_energy = ledger%initial_energy_J_m3 - heat(1) - heat(2)
    call check(close_relative(ledger%final_energy_J_m3, expected_energy, &
         1.0e-10_c_double), 'two-bath heat closes the energy ledger')
    call check(old(1) == 0.0_c_double .and. old(2) == 1.0e20_c_double .and. &
         old(3) == 0.0_c_double, 'two-bath trial leaves old input unchanged')
  end subroutine verify_two_bath_exchange

  subroutine verify_noncontiguous_sections()
    real(c_double) :: edges_store(5), old_store(4), birth_store(4)
    real(c_double) :: escape_store(4), trial_store(4)
    real(c_double) :: bath_store(4), diffusion_store(7), heat_store(4)
    real(c_double) :: edges_ref(3), old_ref(2), bath_ref(2), diffusion_ref(2)
    real(c_double) :: birth_ref(2), escape_ref(2), trial_ref(2), heat_ref(2)
    type(fusion_kinetic_ledger_v1) :: ledger, ledger_ref
    integer(c_int) :: status, status_ref

    ! The wrapper's CONTIGUOUS dummies must make safe compiler temporaries
    ! for these strided sections before taking C_LOC addresses.
    edges_store = [1.0e-16_c_double, -1.0_c_double, 2.0e-16_c_double, &
         -1.0_c_double, 4.0e-16_c_double]
    old_store = [0.0_c_double, -1.0_c_double, 1.0e20_c_double, &
         -1.0_c_double]
    birth_store = 0.0_c_double
    escape_store = 0.0_c_double
    trial_store = -1.0_c_double
    bath_store = [1.0e-16_c_double, -1.0_c_double, 8.0e-16_c_double, &
         -1.0_c_double]
    diffusion_store = -1.0_c_double
    diffusion_store(1) = diffusion
    diffusion_store(3) = diffusion
    diffusion_store(5) = diffusion
    diffusion_store(7) = diffusion
    heat_store = -1.0_c_double

    call fusion_energy_fp_trial(2_c_int, 2_c_int, 0.05_c_double, &
         edges_store(1:5:2), old_store(1:4:2), bath_store(1:4:2), &
         diffusion_store(1:3:2), birth_store(1:4:2), escape_store(1:4:2), &
         0.0_c_double, trial_store(1:4:2), heat_store(1:4:2), ledger, status)

    edges_ref = [1.0e-16_c_double, 2.0e-16_c_double, 4.0e-16_c_double]
    old_ref = [0.0_c_double, 1.0e20_c_double]
    bath_ref = [1.0e-16_c_double, 8.0e-16_c_double]
    diffusion_ref = [diffusion, diffusion]
    birth_ref = 0.0_c_double
    escape_ref = 0.0_c_double
    trial_ref = -1.0_c_double
    heat_ref = -1.0_c_double
    call fusion_energy_fp_trial(2_c_int, 2_c_int, 0.05_c_double, edges_ref, &
         old_ref, bath_ref, diffusion_ref, birth_ref, escape_ref, 0.0_c_double, &
         trial_ref, heat_ref, ledger_ref, status_ref)
    call check(status == PB11_STATUS_OK .and. status_ref == PB11_STATUS_OK, &
         'strided sections and contiguous reference return OK')
    call check(close_relative(trial_store(1), trial_ref(1), 1.0e-14_c_double) &
         .and. close_relative(trial_store(3), trial_ref(2), 1.0e-14_c_double), &
         'strided trial section matches contiguous reference')
    call check(close_relative(heat_store(1), heat_ref(1), 1.0e-14_c_double) &
         .and. close_relative(heat_store(3), heat_ref(2), 1.0e-14_c_double), &
         'strided bath-major heat section matches contiguous reference')
    call check(old_store(1) == 0.0_c_double .and. &
         old_store(3) == 1.0e20_c_double, &
         'strided old section remains unchanged')
  end subroutine verify_noncontiguous_sections

  subroutine verify_input_errors()
    real(c_double) :: edges(2), old(1), birth(1), escape(1), trial(1)
    real(c_double) :: bath(0), diffusion_values(0), heat(0)
    real(c_double) :: wrong_edges(2), wrong_old(1), wrong_birth(1)
    real(c_double) :: wrong_escape(1), wrong_trial(2), wrong_heat(0)
    real(c_double) :: nan_old(1), one_bath(0), one_diffusion(0), one_heat(1)
    type(fusion_kinetic_ledger_v1) :: ledger
    integer(c_int) :: status

    edges = [1.0e-16_c_double, 3.0e-16_c_double]
    old = [1.0e20_c_double]
    birth = 0.0_c_double
    escape = 0.0_c_double
    trial = -1.0_c_double
    call fusion_energy_fp_trial(1_c_int, 0_c_int, 0.0_c_double, edges, old, &
         bath, diffusion_values, birth, escape, 0.0_c_double, trial, heat, &
         ledger, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE, &
         'zero time step returns out-of-range status')
    call check(trial(1) == 0.0_c_double .and. ledger_is_zero(ledger), &
         'time-step error clears outputs')

    nan_old(1) = ieee_value(0.0_c_double, ieee_quiet_nan)
    trial = -1.0_c_double
    call fusion_energy_fp_trial(1_c_int, 0_c_int, 0.1_c_double, edges, nan_old, &
         bath, diffusion_values, birth, escape, 0.0_c_double, trial, heat, &
         ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'NaN population returns invalid-argument status')
    call check(trial(1) == 0.0_c_double .and. ledger_is_zero(ledger), &
         'NaN input clears outputs')

    ! Wrong extents are rejected before any C pointer is formed.
    wrong_edges = edges
    wrong_old = old
    wrong_birth = birth
    wrong_escape = escape
    wrong_trial = -1.0_c_double
    call fusion_energy_fp_trial(2_c_int, 0_c_int, 0.1_c_double, wrong_edges, &
         wrong_old, bath, diffusion_values, wrong_birth, wrong_escape, &
         0.0_c_double, wrong_trial, wrong_heat, ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong Fortran array extents are rejected')
    call check(all(wrong_trial == 0.0_c_double) .and. ledger_is_zero(ledger), &
         'array extent error clears outputs')

    one_heat = -1.0_c_double
    trial = -1.0_c_double
    call fusion_energy_fp_trial(1_c_int, 1_c_int, 0.1_c_double, edges, old, &
         one_bath, one_diffusion, birth, escape, 0.0_c_double, trial, &
         one_heat, ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'missing one-cell bath extent is rejected')
    call check(trial(1) == 0.0_c_double .and. one_heat(1) == 0.0_c_double .and. &
         ledger_is_zero(ledger), 'missing bath input clears outputs')
  end subroutine verify_input_errors

end program test_fusion_kinetics_fortran
