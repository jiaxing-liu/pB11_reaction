program test_fusion_coupled_thermal_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_null_ptr, c_ptr, &
       c_size_t, c_sizeof
  use fusion_coupled_thermal_fortran
  use fusion_thermal_birth_fortran, only : &
       fusion_thermal_birth_options_v1, FUSION_ENDPOINT_S, FUSION_PB_LOW_TB
  use fusion_source_state_fortran, only : fusion_source_ledger_v1
  implicit none

  real(c_double), parameter :: kev_j = 1.602176634e-16_c_double
  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  real(c_double), parameter :: kT_j = 1.0_c_double * kev_j
  real(c_double), parameter :: u_kg = 1.66053906660e-27_c_double
  integer(c_int), parameter :: cells = 3_c_int
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_no_reaction_with_inert()
  call verify_no_inert_null_pointer()
  call verify_wrong_extents_clear_outputs()
  call verify_table_wrong_extent_clear_outputs()
  call verify_table_all_channels_closed_identity()
  call verify_table_dt_null_rejected()
  call verify_increment_mapping()
  call verify_increment_wrong_extent_clear()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, &
          ' Fortran coupled-thermal test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran coupled-thermal tests passed'

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
    near = abs(actual - expected) <= 2.0e-12_c_double * scale
  end function near

  subroutine fill_ledger(ledger, value)
    type(fusion_source_ledger_v1), intent(out) :: ledger
    real(c_double), intent(in) :: value

    ledger%events_m3 = value
    ledger%nuclear_born_number_m3 = value
    ledger%nuclear_born_energy_J_m3 = value
    ledger%external_born_number_m3 = value
    ledger%external_born_energy_J_m3 = value
    ledger%thermal_consumed_number_m3 = value
    ledger%thermal_consumed_energy_J_m3 = value
    ledger%fast_consumed_number_m3 = value
    ledger%fast_consumed_energy_J_m3 = value
    ledger%escaped_number_m3 = value
    ledger%escaped_energy_J_m3 = value
    ledger%handed_off_number_m3 = value
    ledger%handed_off_energy_J_m3 = value
    ledger%heat_to_bath_J_m3 = value
    ledger%neutron_number_m3 = value
    ledger%neutron_energy_J_m3 = value
  end subroutine fill_ledger

  subroutine fill_coupled(out, value)
    type(fusion_coupled_thermal_v1), intent(out) :: out
    real(c_double), intent(in) :: value

    call fill_ledger(out%ledger, value)
    out%inert_ion_heat_J_m3 = value
    out%electron_energy_J_m3 = value
    out%ion_energy_J_m3 = value
    out%particle_residual_m3 = value
    out%energy_residual_J_m3 = value
    out%max_source_rate_discrepancy = value
    out%max_source_debit_discrepancy = value
    out%handoff_L1 = value
    out%handoff_mean_error = value
    out%handoff_projected = int(value, c_int)
  end subroutine fill_coupled

  subroutine fill_increment(out, value)
    type(fusion_thermal_increment_v1), intent(out) :: out
    real(c_double), intent(in) :: value

    out%thermal_number_m3 = value
    out%electron_energy_J_m3 = value
    out%ion_energy_J_m3 = value
  end subroutine fill_increment

  logical function ledger_is_zero(ledger)
    type(fusion_source_ledger_v1), intent(in) :: ledger

    ledger_is_zero = all(ledger%events_m3 == 0.0_c_double) .and. &
         all(ledger%nuclear_born_number_m3 == 0.0_c_double) .and. &
         all(ledger%nuclear_born_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%external_born_number_m3 == 0.0_c_double) .and. &
         all(ledger%external_born_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%thermal_consumed_number_m3 == 0.0_c_double) .and. &
         all(ledger%thermal_consumed_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%fast_consumed_number_m3 == 0.0_c_double) .and. &
         all(ledger%fast_consumed_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%escaped_number_m3 == 0.0_c_double) .and. &
         all(ledger%escaped_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%handed_off_number_m3 == 0.0_c_double) .and. &
         all(ledger%handed_off_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%heat_to_bath_J_m3 == 0.0_c_double) .and. &
         ledger%neutron_number_m3 == 0.0_c_double .and. &
         ledger%neutron_energy_J_m3 == 0.0_c_double
  end function ledger_is_zero

  logical function coupled_is_zero(out)
    type(fusion_coupled_thermal_v1), intent(in) :: out

    coupled_is_zero = ledger_is_zero(out%ledger) .and. &
         all(out%inert_ion_heat_J_m3 == 0.0_c_double) .and. &
         out%electron_energy_J_m3 == 0.0_c_double .and. &
         out%ion_energy_J_m3 == 0.0_c_double .and. &
         all(out%particle_residual_m3 == 0.0_c_double) .and. &
         out%energy_residual_J_m3 == 0.0_c_double .and. &
         out%max_source_rate_discrepancy == 0.0_c_double .and. &
         out%max_source_debit_discrepancy == 0.0_c_double .and. &
         all(out%handoff_L1 == 0.0_c_double) .and. &
         all(out%handoff_mean_error == 0.0_c_double) .and. &
         all(out%handoff_projected == 0_c_int)
  end function coupled_is_zero

  subroutine make_options(options)
    type(fusion_coupled_thermal_options_v1), intent(out) :: options

    options%birth%relative_max_J = 100.0_c_double * kev_j
    options%birth%cm_max_kT = 40.0_c_double
    options%birth%ground_state_q_J = 91.84_c_double * kev_j
    options%birth%cutoff_J = 0.001_c_double * mev_j
    options%birth%l1_fraction = 0.76_c_double
    options%birth%relative_phase = 0.0_c_double
    options%birth%narrow_peak_fraction = 0.051_c_double
    options%birth%continuum_peak_scale = 1.0_c_double
    options%birth%continuation = FUSION_ENDPOINT_S
    options%birth%pb_low = FUSION_PB_LOW_TB
    options%birth%remainder_policy = 0_c_int
    options%birth%broad_mode = 13_c_int
    options%birth%fsci_policy = 0_c_int
    options%birth%relative_order = 8_c_int
    options%birth%cm_order = 8_c_int
    options%birth%nq = 4_c_int
    options%birth%ncos = 8_c_int
    options%max_source_rate_error = 1.0e-6_c_double
    options%max_source_debit_error = 1.0e-6_c_double
    options%handoff_max_L1 = 1.0e-3_c_double
    options%handoff_max_mean_error = 1.0e-3_c_double
    options%channels = 0_c_int
    options%handoff_enabled = 0_c_int
  end subroutine make_options

  subroutine make_null_tables(tables)
    type(c_ptr), intent(out) :: tables(:)

    tables = c_null_ptr
  end subroutine make_null_tables

  subroutine make_common_inputs(edges, thermal_number, charge_squared, &
       inert, logs, old_s, old_t, external_birth, escape, trial_thermal, &
       trial_s, trial_t, electron_energy, ion_energy, electron_density)
    real(c_double), intent(out), target :: edges(:)
    real(c_double), intent(out), target :: thermal_number(:)
    real(c_double), intent(out), target :: charge_squared(:)
    type(fusion_inert_ion_v1), intent(out), target :: inert(:)
    real(c_double), intent(out), target :: logs(:,:)
    real(c_double), intent(out), target :: old_s(:,:), old_t(:,:)
    real(c_double), intent(out), target :: external_birth(:,:), escape(:,:)
    real(c_double), intent(out), target :: trial_thermal(:)
    real(c_double), intent(out), target :: trial_s(:,:), trial_t(:,:)
    real(c_double), intent(out) :: electron_energy, ion_energy
    real(c_double), intent(out) :: electron_density

    integer :: i
    real(c_double) :: thermal_total, electron_kT

    do i = 1, size(edges)
       edges(i) = 20.0_c_double * mev_j * real(i - 1, c_double) / &
            real(size(edges) - 1, c_double)
    end do
    thermal_number = 0.0_c_double
    thermal_number(5) = 1.0e19_c_double
    charge_squared = 0.0_c_double
    charge_squared(5) = 4.0_c_double
    inert(1)%density_m3 = 2.0e19_c_double
    inert(1)%mass_kg = 12.0_c_double * u_kg
    inert(1)%mean_charge_squared = 16.0_c_double
    logs = 10.0_c_double
    old_s = 0.0_c_double
    old_t = 0.0_c_double
    external_birth = 0.0_c_double
    escape = 0.0_c_double
    trial_thermal = -7.0_c_double
    trial_s = -7.0_c_double
    trial_t = -7.0_c_double

    thermal_total = sum(thermal_number) + inert(1)%density_m3
    electron_density = 1.0e19_c_double
    electron_kT = kT_j
    electron_energy = 1.5_c_double * electron_kT * electron_density
    ion_energy = 1.5_c_double * kT_j * thermal_total
  end subroutine make_common_inputs

  subroutine verify_layout()
    type(fusion_coupled_thermal_options_v1) :: options
    type(fusion_coupled_thermal_v1) :: out
    type(fusion_thermal_increment_v1) :: increment
    type(fusion_inert_ion_v1) :: inert
    type(fusion_thermal_birth_options_v1) :: birth
    type(fusion_source_ledger_v1) :: ledger
    real(c_double) :: dummy

    call check(c_sizeof(birth) == 104_c_size_t, &
         'nested thermal-birth options retain their C ABI size')
    call check(c_sizeof(ledger) == 968_c_size_t, &
         'nested source ledger retains its 968-byte C ABI size')
    call check(c_sizeof(inert) == 24_c_size_t, &
         'inert-ion struct has three c_double fields')
    call check(c_sizeof(options) == 160_c_size_t .and. &
         c_sizeof(options) == 20_c_size_t * c_sizeof(dummy), &
         'coupled options have the 160-byte C ABI layout')
    call check(c_sizeof(out) == 1224_c_size_t .and. &
         c_sizeof(out) == 153_c_size_t * c_sizeof(dummy), &
         'coupled result has the 1224-byte C ABI layout')
    call check(c_sizeof(increment) == 64_c_size_t .and. &
         c_sizeof(increment) == 8_c_size_t * c_sizeof(dummy), &
         'thermal increment has the 64-byte C ABI layout')
  end subroutine verify_layout

  subroutine verify_no_reaction_with_inert()
    real(c_double), target :: edges(cells + 1)
    real(c_double), target :: thermal_number(6), charge_squared(6)
    type(fusion_inert_ion_v1), target :: inert(1)
    real(c_double), target :: logs(8,6)
    real(c_double), target :: old_s(cells,6), old_t(cells,6)
    real(c_double), target :: external_birth(cells,6), escape(cells,6)
    real(c_double), target :: trial_thermal(6), trial_s(cells,6), trial_t(cells,6)
    type(fusion_coupled_thermal_options_v1), target :: options
    type(fusion_coupled_thermal_v1), target :: out
    real(c_double) :: electron_energy, ion_energy, electron_density
    integer(c_int) :: status

    call make_options(options)
    call make_common_inputs(edges, thermal_number, charge_squared, inert, logs, &
         old_s, old_t, external_birth, escape, trial_thermal, trial_s, trial_t, &
         electron_energy, ion_energy, electron_density)
    call fill_coupled(out, -9.0_c_double)

    call fusion_coupled_thermal_trial(1.0e-6_c_double, options, cells, edges, &
         thermal_number, electron_energy, ion_energy, electron_density, &
         charge_squared, 1_c_int, inert, logs, old_s, old_t, external_birth, &
         escape, trial_thermal, trial_s, trial_t, out, status)

    call check(status == PB11_STATUS_OK, &
         'no-reaction inert-ion coupled trial returns OK')
    call check(all(trial_thermal == thermal_number), &
         'no-reaction trial preserves thermal network numbers')
    call check(all(trial_s == 0.0_c_double) .and. &
         all(trial_t == 0.0_c_double), &
         'no-reaction trial preserves empty fast populations')
    call check(near(out%electron_energy_J_m3, electron_energy) .and. &
         near(out%ion_energy_J_m3, ion_energy), &
         'no-reaction trial preserves electron and ion energy')
    call check(coupled_is_zero_except_energy(out), &
         'no-reaction trial has no ledger, heat, residual, or handoff amount')
  end subroutine verify_no_reaction_with_inert

  logical function coupled_is_zero_except_energy(out)
    type(fusion_coupled_thermal_v1), intent(in) :: out

    coupled_is_zero_except_energy = ledger_is_zero(out%ledger) .and. &
         all(out%inert_ion_heat_J_m3 == 0.0_c_double) .and. &
         all(out%particle_residual_m3 == 0.0_c_double) .and. &
         out%energy_residual_J_m3 == 0.0_c_double .and. &
         out%max_source_rate_discrepancy == 0.0_c_double .and. &
         out%max_source_debit_discrepancy == 0.0_c_double .and. &
         all(out%handoff_L1 == 0.0_c_double) .and. &
         all(out%handoff_mean_error == 0.0_c_double) .and. &
         all(out%handoff_projected == 0_c_int)
  end function coupled_is_zero_except_energy

  subroutine verify_no_inert_null_pointer()
    real(c_double), target :: edges(cells + 1), thermal_number(6)
    real(c_double), target :: charge_squared(6), logs(7,6)
    real(c_double), target :: old_s(cells,6), old_t(cells,6)
    real(c_double), target :: external_birth(cells,6), escape(cells,6)
    real(c_double), target :: trial_thermal(6), trial_s(cells,6), trial_t(cells,6)
    type(fusion_inert_ion_v1), target :: no_inert(0)
    type(fusion_coupled_thermal_options_v1), target :: options
    type(fusion_coupled_thermal_v1), target :: out
    real(c_double) :: electron_energy, ion_energy, electron_density
    integer(c_int) :: status
    integer :: i

    call make_options(options)
    do i = 1, size(edges)
       edges(i) = 20.0_c_double * mev_j * real(i - 1, c_double) / &
            real(size(edges) - 1, c_double)
    end do
    thermal_number = 0.0_c_double
    thermal_number(5) = 1.0e19_c_double
    charge_squared = 0.0_c_double
    charge_squared(5) = 4.0_c_double
    logs = 10.0_c_double
    old_s = 0.0_c_double
    old_t = 0.0_c_double
    external_birth = 0.0_c_double
    escape = 0.0_c_double
    trial_thermal = -7.0_c_double
    trial_s = -7.0_c_double
    trial_t = -7.0_c_double
    electron_density = 1.0e19_c_double
    electron_energy = 1.5_c_double * kT_j * electron_density
    ion_energy = 1.5_c_double * kT_j * sum(thermal_number)
    call fill_coupled(out, -9.0_c_double)

    call fusion_coupled_thermal_trial(1.0e-6_c_double, options, cells, edges, &
         thermal_number, electron_energy, ion_energy, electron_density, &
         charge_squared, 0_c_int, no_inert, logs, old_s, old_t, &
         external_birth, escape, trial_thermal, trial_s, trial_t, out, status)

    call check(status == PB11_STATUS_OK, &
         'zero inert-count trial accepts a zero-size inert array')
    call check(all(trial_thermal == thermal_number) .and. &
         all(trial_s == 0.0_c_double) .and. all(trial_t == 0.0_c_double), &
         'zero inert-count trial preserves the no-reaction state')
    call check(near(out%electron_energy_J_m3, electron_energy) .and. &
         near(out%ion_energy_J_m3, ion_energy), &
         'zero inert-count trial preserves both thermal energies')
  end subroutine verify_no_inert_null_pointer

  subroutine verify_wrong_extents_clear_outputs()
    real(c_double), target :: edges(cells + 1), bad_edges(cells)
    real(c_double), target :: thermal_number(6), charge_squared(6)
    type(fusion_inert_ion_v1), target :: inert(1)
    real(c_double), target :: logs(8,6), bad_logs(7,6)
    real(c_double), target :: old_s(cells,6), old_t(cells,6)
    real(c_double), target :: external_birth(cells,6), escape(cells,6)
    real(c_double), target :: trial_thermal(6), trial_s(cells,6), trial_t(cells,6)
    real(c_double), target :: bad_trial_t(cells-1,6)
    type(fusion_coupled_thermal_options_v1), target :: options
    type(fusion_coupled_thermal_v1), target :: out
    real(c_double) :: electron_energy, ion_energy, electron_density
    integer(c_int) :: status

    call make_options(options)
    call make_common_inputs(edges, thermal_number, charge_squared, inert, logs, &
         old_s, old_t, external_birth, escape, trial_thermal, trial_s, trial_t, &
         electron_energy, ion_energy, electron_density)
    bad_edges = 0.0_c_double
    bad_logs = 0.0_c_double
    bad_trial_t = -4.0_c_double

    call fill_coupled(out, -4.0_c_double)
    call fusion_coupled_thermal_trial(1.0e-6_c_double, options, cells, &
         bad_edges, thermal_number, electron_energy, ion_energy, &
         electron_density, charge_squared, 1_c_int, inert, logs, old_s, old_t, &
         external_birth, escape, trial_thermal, trial_s, trial_t, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong edge extent is rejected before C_LOC')
    call check(all(trial_thermal == 0.0_c_double) .and. &
         all(trial_s == 0.0_c_double) .and. all(trial_t == 0.0_c_double) .and. &
         coupled_is_zero(out), &
         'wrong edge extent clears all output arrays and result fields')

    trial_thermal = -5.0_c_double
    trial_s = -5.0_c_double
    trial_t = -5.0_c_double
    call fill_coupled(out, -5.0_c_double)
    call fusion_coupled_thermal_trial(1.0e-6_c_double, options, cells, edges, &
         thermal_number, electron_energy, ion_energy, electron_density, &
         charge_squared, 1_c_int, inert, bad_logs, old_s, old_t, &
         external_birth, escape, trial_thermal, trial_s, trial_t, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong Coulomb-log extent is rejected before C_LOC')
    call check(all(trial_thermal == 0.0_c_double) .and. &
         all(trial_s == 0.0_c_double) .and. all(trial_t == 0.0_c_double) .and. &
         coupled_is_zero(out), &
         'wrong Coulomb-log extent clears all outputs')

    trial_thermal = -6.0_c_double
    trial_s = -6.0_c_double
    bad_trial_t = -6.0_c_double
    call fill_coupled(out, -6.0_c_double)
    call fusion_coupled_thermal_trial(1.0e-6_c_double, options, cells, edges, &
         thermal_number, electron_energy, ion_energy, electron_density, &
         charge_squared, 1_c_int, inert, logs, old_s, old_t, external_birth, &
         escape, trial_thermal, trial_s, bad_trial_t, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong fast-trial extent is rejected before C_LOC')
    call check(all(trial_thermal == 0.0_c_double) .and. &
         all(trial_s == 0.0_c_double) .and. all(bad_trial_t == 0.0_c_double) .and. &
         coupled_is_zero(out), &
         'wrong fast-trial extent clears every actual output')
  end subroutine verify_wrong_extents_clear_outputs

  subroutine verify_table_wrong_extent_clear_outputs()
    real(c_double), target :: edges(cells + 1)
    real(c_double), target :: thermal_number(6), charge_squared(6)
    type(fusion_inert_ion_v1), target :: inert(1)
    real(c_double), target :: logs(8,6)
    real(c_double), target :: old_s(cells,6), old_t(cells,6)
    real(c_double), target :: external_birth(cells,6), escape(cells,6)
    real(c_double), target :: trial_thermal(6), trial_s(cells,6), trial_t(cells,6)
    type(c_ptr), target :: tables(4)
    type(fusion_coupled_thermal_options_v1), target :: options
    type(fusion_coupled_thermal_v1), target :: out
    real(c_double) :: electron_energy, ion_energy, electron_density
    integer(c_int) :: status

    call make_options(options)
    call make_common_inputs(edges, thermal_number, charge_squared, inert, logs, &
         old_s, old_t, external_birth, escape, trial_thermal, trial_s, trial_t, &
         electron_energy, ion_energy, electron_density)
    call make_null_tables(tables)
    call fill_coupled(out, -11.0_c_double)

    call fusion_coupled_thermal_table_trial(1.0e-6_c_double, options, tables, &
         cells, edges, thermal_number, electron_energy, ion_energy, &
         electron_density, charge_squared, 1_c_int, inert, logs, old_s, old_t, &
         external_birth, escape, trial_thermal, trial_s, trial_t, out, status)

    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'table trial rejects a handle array with the wrong extent')
    call check(all(trial_thermal == 0.0_c_double) .and. &
         all(trial_s == 0.0_c_double) .and. all(trial_t == 0.0_c_double) .and. &
         coupled_is_zero(out), &
         'wrong table extent clears all trial arrays and result fields')
  end subroutine verify_table_wrong_extent_clear_outputs

  subroutine verify_table_all_channels_closed_identity()
    real(c_double), target :: edges(cells + 1)
    real(c_double), target :: thermal_number(6), charge_squared(6)
    type(fusion_inert_ion_v1), target :: inert(1)
    real(c_double), target :: logs(8,6)
    real(c_double), target :: old_s(cells,6), old_t(cells,6)
    real(c_double), target :: external_birth(cells,6), escape(cells,6)
    real(c_double), target :: trial_thermal(6), trial_s(cells,6), trial_t(cells,6)
    type(c_ptr), target :: tables(5)
    type(fusion_coupled_thermal_options_v1), target :: options
    type(fusion_coupled_thermal_v1), target :: out
    real(c_double) :: electron_energy, ion_energy, electron_density
    integer(c_int) :: status

    call make_options(options)
    call make_common_inputs(edges, thermal_number, charge_squared, inert, logs, &
         old_s, old_t, external_birth, escape, trial_thermal, trial_s, trial_t, &
         electron_energy, ion_energy, electron_density)
    call make_null_tables(tables)
    call fill_coupled(out, -12.0_c_double)

    call fusion_coupled_thermal_table_trial(1.0e-6_c_double, options, tables, &
         cells, edges, thermal_number, electron_energy, ion_energy, &
         electron_density, charge_squared, 1_c_int, inert, logs, old_s, old_t, &
         external_birth, escape, trial_thermal, trial_s, trial_t, out, status)

    call check(status == PB11_STATUS_OK, &
         'all-disabled table trial with five null handles returns OK')
    call check(all(trial_thermal == thermal_number) .and. &
         all(trial_s == 0.0_c_double) .and. all(trial_t == 0.0_c_double), &
         'all-disabled table trial preserves the thermal/kinetic identity state')
    call check(near(out%electron_energy_J_m3, electron_energy) .and. &
         near(out%ion_energy_J_m3, ion_energy) .and. &
         coupled_is_zero_except_energy(out), &
         'all-disabled table trial preserves both energies and empty ledgers')
  end subroutine verify_table_all_channels_closed_identity

  subroutine verify_table_dt_null_rejected()
    real(c_double), target :: edges(cells + 1)
    real(c_double), target :: thermal_number(6), charge_squared(6)
    type(fusion_inert_ion_v1), target :: inert(1)
    real(c_double), target :: logs(8,6)
    real(c_double), target :: old_s(cells,6), old_t(cells,6)
    real(c_double), target :: external_birth(cells,6), escape(cells,6)
    real(c_double), target :: trial_thermal(6), trial_s(cells,6), trial_t(cells,6)
    type(c_ptr), target :: tables(5)
    type(fusion_coupled_thermal_options_v1), target :: options
    type(fusion_coupled_thermal_v1), target :: out
    real(c_double) :: electron_energy, ion_energy, electron_density
    integer(c_int) :: status

    call make_options(options)
    call make_common_inputs(edges, thermal_number, charge_squared, inert, logs, &
         old_s, old_t, external_birth, escape, trial_thermal, trial_s, trial_t, &
         electron_energy, ion_energy, electron_density)
    call make_null_tables(tables)
    options%channels = 0_c_int
    options%channels(4) = 1_c_int
    thermal_number = 0.0_c_double
    thermal_number(2) = 1.0e19_c_double
    thermal_number(3) = 1.0e19_c_double
    charge_squared = 0.0_c_double
    charge_squared(2) = 1.0_c_double
    charge_squared(3) = 1.0_c_double
    ion_energy = 1.5_c_double * kT_j * &
         (sum(thermal_number) + inert(1)%density_m3)
    call fill_coupled(out, -13.0_c_double)

    call fusion_coupled_thermal_table_trial(1.0e-6_c_double, options, tables, &
         cells, edges, thermal_number, electron_energy, ion_energy, &
         electron_density, charge_squared, 1_c_int, inert, logs, old_s, old_t, &
         external_birth, escape, trial_thermal, trial_s, trial_t, out, status)

    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'enabled DT channel with a null table is rejected without fallback')
    call check(all(trial_thermal == 0.0_c_double) .and. &
         all(trial_s == 0.0_c_double) .and. all(trial_t == 0.0_c_double) .and. &
         coupled_is_zero(out), &
         'missing enabled DT table clears every output')
  end subroutine verify_table_dt_null_rejected

  subroutine verify_increment_mapping()
    type(fusion_source_ledger_v1), target :: ledger
    real(c_double), target :: inert_heat(6)
    type(fusion_thermal_increment_v1), target :: out
    real(c_double) :: expected_numbers(6), expected_electron, expected_ion
    integer :: i, j
    integer(c_int) :: status

    call fill_ledger(ledger, 0.0_c_double)
    ledger%handed_off_number_m3 = 0.0_c_double
    ledger%handed_off_number_m3(2) = 2.5_c_double
    ledger%handed_off_energy_J_m3 = 0.0_c_double
    ledger%handed_off_energy_J_m3(2) = 3.5_c_double
    ledger%heat_to_bath_J_m3 = 0.0_c_double
    do i = 1, 6
       ledger%heat_to_bath_J_m3((i - 1) * 7 + 1) = &
            real(i, c_double)
       do j = 2, 7
          ledger%heat_to_bath_J_m3((i - 1) * 7 + j) = &
               real(i - j, c_double) / 4.0_c_double
       end do
    end do
    inert_heat = [0.5_c_double, -0.75_c_double, 1.25_c_double, &
         -1.5_c_double, 2.0_c_double, -2.25_c_double]
    expected_numbers = ledger%handed_off_number_m3 - &
         ledger%thermal_consumed_number_m3
    expected_electron = 0.0_c_double
    expected_ion = 0.0_c_double
    do i = 1, 6
       expected_electron = expected_electron + &
            ledger%heat_to_bath_J_m3((i - 1) * 7 + 1)
       expected_ion = expected_ion + inert_heat(i) + &
            ledger%handed_off_energy_J_m3(i) - &
            ledger%thermal_consumed_energy_J_m3(i) + &
            sum(ledger%heat_to_bath_J_m3((i - 1) * 7 + 2:(i - 1) * 7 + 7))
    end do
    call fill_increment(out, -99.0_c_double)

    call fusion_coupled_thermal_increment(ledger, inert_heat, out, status)

    call check(status == PB11_STATUS_OK, &
         'signed ledger maps to a thermal increment')
    call check(all(out%thermal_number_m3 == expected_numbers), &
         'increment particle change is handed-off minus thermal-consumed')
    call check(near(out%electron_energy_J_m3, expected_electron), &
         'increment electron energy sums the electron heat columns')
    call check(near(out%ion_energy_J_m3, expected_ion), &
         'increment ion energy counts ion heat, inert heat, and handoff once')
  end subroutine verify_increment_mapping

  subroutine verify_increment_wrong_extent_clear()
    type(fusion_source_ledger_v1), target :: ledger
    real(c_double), target :: bad_inert_heat(5)
    type(fusion_thermal_increment_v1), target :: out
    integer(c_int) :: status

    call fill_ledger(ledger, 0.0_c_double)
    bad_inert_heat = 1.0_c_double
    call fill_increment(out, -77.0_c_double)

    call fusion_coupled_thermal_increment(ledger, bad_inert_heat, out, status)

    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'increment rejects an inert heat array with the wrong extent')
    call check(all(out%thermal_number_m3 == 0.0_c_double) .and. &
         out%electron_energy_J_m3 == 0.0_c_double .and. &
         out%ion_energy_J_m3 == 0.0_c_double, &
         'wrong inert heat extent clears the complete increment')
  end subroutine verify_increment_wrong_extent_clear

end program test_fusion_coupled_thermal_fortran
