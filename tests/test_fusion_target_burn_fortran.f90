program test_fusion_target_burn_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_target_burn_fortran
  use fusion_beam_fortran, only : fusion_beam_window_v1, FUSION_PB11_3ALPHA
  implicit none

  real(c_double), parameter :: keV_J = 1.602176634e-16_c_double
  real(c_double), parameter :: proton_mass_kg = 1.67262192369e-27_c_double
  real(c_double), parameter :: amu_kg = 1.66053906660e-27_c_double
  real(c_double), parameter :: boron11_mass_kg = 11.0_c_double * amu_kg
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_one_bin_trial()
  call verify_negative_energy_rejection()
  call verify_physical_window()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran target-burn test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran target-burn tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function burn_is_zero(out)
    type(fusion_target_burn_v1), intent(in) :: out

    burn_is_zero = out%initial_fast_number_m3 == 0.0_c_double .and. &
         out%final_fast_number_m3 == 0.0_c_double .and. &
         out%initial_fast_energy_J_m3 == 0.0_c_double .and. &
         out%final_fast_energy_J_m3 == 0.0_c_double .and. &
         out%initial_target_number_m3 == 0.0_c_double .and. &
         out%final_target_number_m3 == 0.0_c_double .and. &
         out%initial_target_energy_J_m3 == 0.0_c_double .and. &
         out%final_target_energy_J_m3 == 0.0_c_double .and. &
         out%reactions_m3 == 0.0_c_double .and. &
         out%removed_fast_energy_J_m3 == 0.0_c_double .and. &
         out%removed_target_energy_J_m3 == 0.0_c_double .and. &
         out%fast_number_residual_m3 == 0.0_c_double .and. &
         out%target_number_residual_m3 == 0.0_c_double .and. &
         out%energy_residual_J_m3 == 0.0_c_double
  end function burn_is_zero

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

  logical function close_abs(actual, expected, tolerance)
    real(c_double), intent(in) :: actual, expected, tolerance

    close_abs = ieee_is_finite(actual) .and. ieee_is_finite(expected) .and. &
         abs(actual - expected) <= tolerance
  end function close_abs

  subroutine verify_layout()
    type(fusion_target_burn_v1) :: out
    real(c_double) :: d

    call check(c_sizeof(out) == 14 * c_sizeof(d), &
         'target-burn struct has fourteen c_double fields')
  end subroutine verify_layout

  subroutine verify_one_bin_trial()
    real(c_double), target :: edges(2), old_fast(1), reactivity(1)
    real(c_double), target :: target_moment(1), trial_fast(1), reaction_loss(1)
    real(c_double) :: old_before(1), expected, expected_loss
    real(c_double) :: trial_repeat(1), loss_repeat(1)
    type(fusion_target_burn_v1), target :: out, out_repeat
    integer(c_int) :: status, status_repeat

    edges = [0.0_c_double, 1.0_c_double]
    old_fast = 1.0_c_double
    old_before = old_fast
    reactivity = 1.0_c_double
    target_moment = 0.25_c_double
    expected = (sqrt(5.0_c_double) - 1.0_c_double) / 2.0_c_double
    expected_loss = 1.0_c_double - expected

    call fusion_target_burn_trial(1_c_int, 1.0_c_double, edges, old_fast, &
         1.0_c_double, 1.0_c_double, reactivity, target_moment, trial_fast, &
         reaction_loss, out, status)
    call check(status == PB11_STATUS_OK, 'one-bin trial returns OK')
    call check(close_abs(trial_fast(1), expected, 5.0e-14_c_double), &
         'one-bin Backward-Euler survivor solves the quadratic')
    call check(close_abs(reaction_loss(1), expected_loss, 5.0e-14_c_double), &
         'one-bin reaction loss is the consumed fast population')
    call check(old_fast(1) == old_before(1), 'trial leaves old input unchanged')
    call check(close_abs(out%initial_fast_number_m3, 1.0_c_double, 0.0_c_double) .and. &
         close_abs(out%final_fast_number_m3, expected, 5.0e-14_c_double) .and. &
         close_abs(out%initial_target_number_m3, 1.0_c_double, 0.0_c_double) .and. &
         close_abs(out%final_target_number_m3, expected, 5.0e-14_c_double), &
         'one-bin particle ledger matches the analytic solution')
    call check(close_abs(out%initial_fast_energy_J_m3, 0.5_c_double, 0.0_c_double) .and. &
         close_abs(out%final_fast_energy_J_m3, 0.5_c_double * expected, &
         5.0e-14_c_double), 'one-bin fast-energy ledger uses cell center')
    call check(close_abs(out%removed_target_energy_J_m3, &
         0.25_c_double * expected_loss, 5.0e-14_c_double) .and. &
         abs(out%energy_residual_J_m3) < 1.0e-13_c_double, &
         'one-bin target moment and energy ledger close')

    call fusion_target_burn_trial(1_c_int, 1.0_c_double, edges, old_fast, &
         1.0_c_double, 1.0_c_double, reactivity, target_moment, trial_repeat, &
         loss_repeat, out_repeat, status_repeat)
    call check(status_repeat == PB11_STATUS_OK .and. &
         trial_repeat(1) == trial_fast(1) .and. &
         loss_repeat(1) == reaction_loss(1) .and. &
         out_repeat%final_fast_number_m3 == out%final_fast_number_m3, &
         'repeated trial is deterministic and does not consume old input')
  end subroutine verify_one_bin_trial

  subroutine verify_negative_energy_rejection()
    real(c_double), target :: edges(2), old_fast(1), reactivity(1)
    real(c_double), target :: target_moment(1), trial_fast(1), reaction_loss(1)
    type(fusion_target_burn_v1), target :: out
    integer(c_int) :: status

    edges = [0.0_c_double, 1.0_c_double]
    old_fast = 1.0_c_double
    reactivity = 1.0_c_double
    target_moment = 0.0_c_double
    trial_fast = -1.0_c_double
    reaction_loss = -1.0_c_double
    out%initial_fast_number_m3 = -1.0_c_double
    call fusion_target_burn_trial(1_c_int, 1.0_c_double, edges, old_fast, &
         1.0_c_double, -1.0_c_double, reactivity, target_moment, trial_fast, &
         reaction_loss, out, status)
    call check(status /= PB11_STATUS_OK .and. trial_fast(1) == 0.0_c_double &
         .and. reaction_loss(1) == 0.0_c_double .and. burn_is_zero(out), &
         'negative target energy rejects and clears every output')
  end subroutine verify_negative_energy_rejection

  subroutine verify_physical_window()
    real(c_double), target :: edges(2), old_fast(1), trial_fast(1)
    real(c_double), target :: reaction_loss(1)
    type(fusion_beam_window_v1), target :: windows(1)
    type(fusion_target_burn_v1), target :: out
    integer(c_int) :: status

    edges = [9.0_c_double * keV_J, 11.0_c_double * keV_J]
    old_fast = 1.0e20_c_double
    trial_fast = -1.0_c_double
    reaction_loss = -1.0_c_double
    call fusion_beam_target_burn_window_trial(FUSION_PB11_3ALPHA, 1_c_int, &
         1.0_c_double, proton_mass_kg, boron11_mass_kg, edges, old_fast, &
         1.0e20_c_double, 0.0_c_double, trial_fast, reaction_loss, windows, &
         out, status)
    call check(status == PB11_STATUS_OK, &
         'cold pB window target-burn call returns OK')
    call check(windows(1)%resolved_reactivity_m3_s > 0.0_c_double .and. &
         windows(1)%domain_incomplete == 0_c_int, &
         'cold pB window has a nonzero resolved rate and complete domain')
    call check(ieee_is_finite(windows(1)%quadrature_error_m3_s) .and. &
         out%reactions_m3 >= 0.0_c_double .and. trial_fast(1) >= 0.0_c_double, &
         'physical wrapper returns finite diagnostics and trial arrays')

    windows(1)%resolved_reactivity_m3_s = -1.0_c_double
    out%reactions_m3 = -1.0_c_double
    trial_fast = -1.0_c_double
    reaction_loss = -1.0_c_double
    call fusion_beam_target_burn_window_trial(FUSION_PB11_3ALPHA, 1_c_int, &
         1.0_c_double, proton_mass_kg, boron11_mass_kg, edges, old_fast, &
         1.0e20_c_double, -1.0_c_double, trial_fast, reaction_loss, windows, &
         out, status)
    call check(status /= PB11_STATUS_OK .and. trial_fast(1) == 0.0_c_double &
         .and. reaction_loss(1) == 0.0_c_double .and. window_is_zero(windows(1)) &
         .and. burn_is_zero(out), &
         'negative window temperature rejects and clears outputs')
  end subroutine verify_physical_window

end program test_fusion_target_burn_fortran
