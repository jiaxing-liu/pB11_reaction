program test_fusion_thermal_burn_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use fusion_thermal_burn_fortran
  implicit none

  real(c_double), parameter :: relative_tolerance = 1.0e-10_c_double
  real(c_double), parameter :: absolute_tolerance = 1.0e-12_c_double
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_pb11_and_dd_stoichiometry()
  call verify_energy_rejection()
  call verify_bad_extents()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran thermal-burn test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran thermal-burn tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function near(actual, expected)
    real(c_double), intent(in) :: actual, expected
    real(c_double) :: scale

    scale = max(abs(actual), abs(expected), 1.0_c_double)
    near = abs(actual - expected) <= absolute_tolerance + &
         relative_tolerance * scale
  end function near

  logical function vector_near(actual, expected)
    real(c_double), intent(in) :: actual(:), expected(:)

    vector_near = size(actual) == size(expected) .and. &
         all(abs(actual - expected) <= absolute_tolerance + &
         relative_tolerance * max(abs(actual), abs(expected), 1.0_c_double))
  end function vector_near

  logical function ledger_near_zero(out)
    type(fusion_thermal_burn_v1), intent(in) :: out
    real(c_double) :: largest

    largest = max(maxval(abs(out%events_m3)), &
         maxval(abs(out%reactant_removed_m3)))
    largest = max(largest, maxval(abs(out%reactant_removed_energy_J_m3)))
    largest = max(largest, maxval(abs(out%fast_product_birth_m3)))
    largest = max(largest, abs(out%neutron_birth_m3))
    largest = max(largest, maxval(abs(out%number_residual_m3)))
    largest = max(largest, maxval(abs(out%energy_residual_J_m3)))
    ledger_near_zero = largest <= absolute_tolerance
  end function ledger_near_zero

  subroutine verify_layout()
    type(fusion_thermal_burn_v1) :: ledger
    real(c_double) :: dummy

    call check(c_sizeof(ledger) == 36_c_size_t * c_sizeof(dummy), &
         'thermal-burn ledger has 36 consecutive c_double values')
  end subroutine verify_layout

  subroutine verify_pb11_and_dd_stoichiometry()
    real(c_double) :: old_number(FUSION_SPECIES_COUNT)
    real(c_double) :: old_energy(FUSION_SPECIES_COUNT)
    real(c_double) :: old_number_before(FUSION_SPECIES_COUNT)
    real(c_double) :: old_energy_before(FUSION_SPECIES_COUNT)
    real(c_double) :: reactivity(FUSION_CHANNEL_COUNT)
    real(c_double) :: mean_a(FUSION_CHANNEL_COUNT)
    real(c_double) :: mean_b(FUSION_CHANNEL_COUNT)
    real(c_double) :: trial_number(FUSION_SPECIES_COUNT)
    real(c_double) :: trial_energy(FUSION_SPECIES_COUNT)
    real(c_double) :: pB_events, dd_events
    type(fusion_thermal_burn_v1) :: out
    integer(c_int) :: status
    integer :: p, d, t, he4, b, i

    p = FUSION_PROTON + 1
    d = FUSION_DEUTERON + 1
    t = FUSION_TRITON + 1
    he4 = FUSION_HELIUM4 + 1
    b = FUSION_BORON11 + 1

    old_number = 0.0_c_double
    old_energy = 0.0_c_double
    old_number(p) = 1.0e19_c_double
    old_number(d) = 2.0e19_c_double
    old_number(b) = 1.0e19_c_double
    old_energy(p) = 1.0e4_c_double
    old_energy(d) = 1.0e4_c_double
    old_energy(b) = 1.0e4_c_double
    old_number_before = old_number
    old_energy_before = old_energy

    reactivity = 0.0_c_double
    reactivity(FUSION_PB11_3ALPHA + 1) = 1.0e-22_c_double
    reactivity(FUSION_DD_TP + 1) = 1.0e-22_c_double
    mean_a = 0.0_c_double
    mean_b = 0.0_c_double
    mean_a(FUSION_PB11_3ALPHA + 1) = 2.0e-16_c_double
    mean_b(FUSION_PB11_3ALPHA + 1) = 3.0e-16_c_double
    mean_a(FUSION_DD_TP + 1) = 4.0e-16_c_double
    mean_b(FUSION_DD_TP + 1) = 4.0e-16_c_double
    trial_number = -1.0_c_double
    trial_energy = -1.0_c_double

    call fusion_thermal_burn_trial(1.0_c_double, old_number, old_energy, &
         reactivity, mean_a, mean_b, trial_number, trial_energy, out, status)
    call check(status == PB11_STATUS_OK, &
         'pB/DD reaction trial returns OK')

    pB_events = out%events_m3(FUSION_PB11_3ALPHA + 1)
    dd_events = out%events_m3(FUSION_DD_TP + 1)
    call check(pB_events > 0.0_c_double .and. dd_events > 0.0_c_double, &
         'pB and DD channels both produce positive event amounts')
    call check(near(out%reactant_removed_m3(p), pB_events), &
         'pB removes one proton per event')
    call check(near(out%reactant_removed_m3(b), pB_events), &
         'pB removes one boron per event')
    call check(near(out%reactant_removed_m3(d), 2.0_c_double * dd_events), &
         'DD removes two deuterons per event')
    call check(near(out%fast_product_birth_m3(p), dd_events), &
         'DD(T+p) births one proton per event')
    call check(near(out%fast_product_birth_m3(t), dd_events), &
         'DD(T+p) births one triton per event')
    call check(near(out%fast_product_birth_m3(he4), 3.0_c_double * pB_events), &
         'pB births three alpha particles per event')
    call check(near(out%neutron_birth_m3, 0.0_c_double), &
         'pB plus DD(T+p) produces no neutron')

    do i = 1, FUSION_SPECIES_COUNT
       call check(near(trial_number(i), old_number(i) - &
            out%reactant_removed_m3(i)), &
            'trial number equals old number minus reactant loss')
       call check(near(trial_energy(i), old_energy(i) - &
            out%reactant_removed_energy_J_m3(i)), &
            'trial energy equals old energy minus removed energy')
    end do
    call check(vector_near(old_number, old_number_before), &
         'old number state remains unchanged after trial')
    call check(vector_near(old_energy, old_energy_before), &
         'old energy state remains unchanged after trial')
  end subroutine verify_pb11_and_dd_stoichiometry

  subroutine verify_energy_rejection()
    real(c_double) :: old_number(FUSION_SPECIES_COUNT)
    real(c_double) :: old_energy(FUSION_SPECIES_COUNT)
    real(c_double) :: old_number_before(FUSION_SPECIES_COUNT)
    real(c_double) :: old_energy_before(FUSION_SPECIES_COUNT)
    real(c_double) :: reactivity(FUSION_CHANNEL_COUNT)
    real(c_double) :: mean_a(FUSION_CHANNEL_COUNT)
    real(c_double) :: mean_b(FUSION_CHANNEL_COUNT)
    real(c_double) :: trial_number(FUSION_SPECIES_COUNT)
    real(c_double) :: trial_energy(FUSION_SPECIES_COUNT)
    type(fusion_thermal_burn_v1) :: out
    integer(c_int) :: status
    integer :: p, b

    p = FUSION_PROTON + 1
    b = FUSION_BORON11 + 1
    old_number = 0.0_c_double
    old_energy = 0.0_c_double
    old_number(p) = 1.0e19_c_double
    old_number(b) = 1.0e19_c_double
    old_energy(p) = 1.0e-20_c_double
    old_energy(b) = 1.0e4_c_double
    old_number_before = old_number
    old_energy_before = old_energy
    reactivity = 0.0_c_double
    reactivity(FUSION_PB11_3ALPHA + 1) = 1.0e-22_c_double
    mean_a = 0.0_c_double
    mean_b = 0.0_c_double
    mean_a(FUSION_PB11_3ALPHA + 1) = 1.0e-15_c_double
    mean_b(FUSION_PB11_3ALPHA + 1) = 1.0e-15_c_double
    trial_number = -7.0_c_double
    trial_energy = -7.0_c_double
    out%events_m3 = -7.0_c_double
    out%reactant_removed_m3 = -7.0_c_double
    out%reactant_removed_energy_J_m3 = -7.0_c_double
    out%fast_product_birth_m3 = -7.0_c_double
    out%neutron_birth_m3 = -7.0_c_double
    out%number_residual_m3 = -7.0_c_double
    out%energy_residual_J_m3 = -7.0_c_double

    call fusion_thermal_burn_trial(1.0_c_double, old_number, old_energy, &
         reactivity, mean_a, mean_b, trial_number, trial_energy, out, status)
    call check(status == PB11_STATUS_NUMERICAL_FAILURE, &
         'negative remaining energy rejects the whole trial')
    call check(vector_near(trial_number, 0.0_c_double * trial_number), &
         'energy rejection clears trial number')
    call check(vector_near(trial_energy, 0.0_c_double * trial_energy), &
         'energy rejection clears trial energy')
    call check(ledger_near_zero(out), &
         'energy rejection clears every ledger field')
    call check(vector_near(old_number, old_number_before), &
         'energy rejection leaves old number state unchanged')
    call check(vector_near(old_energy, old_energy_before), &
         'energy rejection leaves old energy state unchanged')
  end subroutine verify_energy_rejection

  subroutine verify_bad_extents()
    real(c_double) :: bad_number(5)
    real(c_double) :: old_energy(FUSION_SPECIES_COUNT)
    real(c_double) :: reactivity(FUSION_CHANNEL_COUNT)
    real(c_double) :: mean_a(FUSION_CHANNEL_COUNT)
    real(c_double) :: mean_b(FUSION_CHANNEL_COUNT)
    real(c_double) :: trial_number(FUSION_SPECIES_COUNT)
    real(c_double) :: trial_energy(FUSION_SPECIES_COUNT)
    real(c_double) :: bad_trial_energy(5)
    type(fusion_thermal_burn_v1) :: out
    integer(c_int) :: status

    bad_number = 1.0_c_double
    old_energy = 1.0_c_double
    reactivity = 0.0_c_double
    mean_a = 0.0_c_double
    mean_b = 0.0_c_double
    trial_number = -3.0_c_double
    trial_energy = -3.0_c_double
    out%events_m3 = -3.0_c_double
    out%reactant_removed_m3 = -3.0_c_double
    out%reactant_removed_energy_J_m3 = -3.0_c_double
    out%fast_product_birth_m3 = -3.0_c_double
    out%neutron_birth_m3 = -3.0_c_double
    out%number_residual_m3 = -3.0_c_double
    out%energy_residual_J_m3 = -3.0_c_double

    call fusion_thermal_burn_trial(1.0_c_double, bad_number, old_energy, &
         reactivity, mean_a, mean_b, trial_number, trial_energy, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrapper rejects a wrong fixed species extent')
    call check(vector_near(trial_number, 0.0_c_double * trial_number) .and. &
         vector_near(trial_energy, 0.0_c_double * trial_energy) .and. &
         ledger_near_zero(out), &
         'species extent failure clears both trial arrays and ledger')

    bad_trial_energy = -4.0_c_double
    trial_number = -4.0_c_double
    out%events_m3 = -4.0_c_double
    out%reactant_removed_m3 = -4.0_c_double
    out%reactant_removed_energy_J_m3 = -4.0_c_double
    out%fast_product_birth_m3 = -4.0_c_double
    out%neutron_birth_m3 = -4.0_c_double
    out%number_residual_m3 = -4.0_c_double
    out%energy_residual_J_m3 = -4.0_c_double
    call fusion_thermal_burn_trial(1.0_c_double, old_energy, old_energy, &
         reactivity, mean_a, mean_b, trial_number, bad_trial_energy, out, &
         status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrapper rejects a wrong trial-energy extent')
    call check(vector_near(trial_number, 0.0_c_double * trial_number) .and. &
         all(abs(bad_trial_energy) <= absolute_tolerance) .and. &
         ledger_near_zero(out), &
         'trial-energy extent failure clears every output')
  end subroutine verify_bad_extents

end program test_fusion_thermal_burn_fortran
