program test_fusion_coupled_sources_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_null_ptr, &
       c_ptr, c_size_t, c_sizeof
  use fusion_coupled_sources_fortran, only : &
       source_trial => fusion_coupled_sources_trial, &
       beam_entry_type => fusion_beam_table_entry_v1, &
       usage_type => fusion_beam_table_usage_v1, &
       fast_options_type => fusion_fast_target_options_v1, &
       coupled_options_type => fusion_coupled_thermal_options_v1, &
       coupled_output_type => fusion_coupled_thermal_v1, &
       diagnostics_type => fusion_handoff_diagnostics_v1, &
       inert_type => fusion_inert_ion_v1, &
       source_ok => PB11_STATUS_OK, source_invalid => PB11_STATUS_INVALID_ARGUMENT
  use fusion_coupled_fast_fortran, only : fast_trial => fusion_coupled_fast_trial
  use fusion_beam_birth_table_fortran, only : &
       fusion_beam_birth_options_v1, fusion_birth_table_control_v1, &
       fusion_beam_birth_table_create, fusion_beam_birth_table_destroy, &
       beam_table_ok => PB11_STATUS_OK, &
       beam_endpoint_s => FUSION_ENDPOINT_S, beam_pb_low_tb => FUSION_PB_LOW_TB, &
       beam_remainder_proxy => FUSION_PB_REMAINDER_ENTRANCE_PROXY
  implicit none

  integer(c_int), parameter :: cells = 32_c_int
  integer(c_int), parameter :: species = 6_c_int
  integer(c_int), parameter :: channels = 5_c_int
  integer(c_int), parameter :: baths = 7_c_int
  integer(c_int), parameter :: table_channel = 3_c_int
  integer(c_int), parameter :: table_slot = 0_c_int
  real(c_double), parameter :: ev_j = 1.602176634e-19_c_double
  real(c_double), parameter :: kev_j = 1.0e3_c_double * ev_j
  real(c_double), parameter :: mev_j = 1.0e6_c_double * ev_j
  real(c_double), parameter :: dt_s = 1.0e-4_c_double
  real(c_double), parameter :: target_energy = 100.0_c_double * kev_j
  real(c_double), parameter :: first_edge = 2.0e-6_c_double * ev_j
  real(c_double), parameter :: max_energy = 25.0_c_double * mev_j

  real(c_double), target :: edges(cells + 1_c_int), edges_bad(1)
  real(c_double), target :: thermal_number(species), thermal_charge_squared(species)
  real(c_double), target :: coulomb_logs(baths,species)
  real(c_double), target :: old_s(cells,species), old_t(cells,species)
  real(c_double), target :: external_birth(cells,species), escape(cells,species)
  real(c_double), target :: trial_number_fast(species), trial_s_fast(cells,species)
  real(c_double), target :: trial_t_fast(cells,species)
  real(c_double), target :: trial_number_source(species), trial_s_source(cells,species)
  real(c_double), target :: trial_t_source(cells,species)
  real(c_double), target :: trial_number_table(species), trial_s_table(cells,species)
  real(c_double), target :: trial_t_table(cells,species)
  real(c_double), target :: trial_number_bad(species), trial_s_bad(cells,species)
  real(c_double), target :: trial_t_bad(cells,species)
  real(c_double), target :: trial_number_shape(species), trial_s_shape(cells,species)
  real(c_double), target :: trial_t_shape(cells,species)
  real(c_double), target :: electron_energy, ion_energy, electron_density

  type(coupled_options_type), target :: options
  type(fast_options_type), target :: fast_options
  type(inert_type), target :: inert(0)
  type(coupled_output_type), target :: output_fast, output_source, output_table
  type(coupled_output_type), target :: output_bad, output_shape
  type(diagnostics_type), target :: diagnostics_source, diagnostics_table
  type(diagnostics_type), target :: diagnostics_bad, diagnostics_shape
  type(usage_type), target :: usage_source, usage_table, usage_bad, usage_shape
  type(c_ptr), target :: thermal_tables(channels)
  type(c_ptr), target :: table
  type(beam_entry_type), target :: entries(1), bad_entries(1)
  type(fusion_beam_birth_options_v1), target :: beam_source
  type(fusion_birth_table_control_v1), target :: beam_control
  integer(c_int) :: status_fast, status_source, status_table, status_bad
  integer(c_int) :: status_shape, status_table_create
  integer :: nearest, failures
  real(c_double) :: ratio, center, distance, best_distance

  failures = 0
  call check(c_sizeof(entries(1)) == 24_c_size_t, &
       'beam-table entry has the expected three-int-plus-pointer layout')
  call check(c_sizeof(usage_source) == 40_c_size_t, &
       'beam-table usage has the expected two-int-plus-four-double layout')

  call make_edges()
  edges_bad = [0.0_c_double]
  thermal_number = 0.0_c_double
  thermal_number(3) = 1.0e19_c_double
  thermal_charge_squared = [1.0_c_double, 1.0_c_double, 1.0_c_double, &
       4.0_c_double, 4.0_c_double, 25.0_c_double]
  coulomb_logs = 15.0_c_double
  old_s = 0.0_c_double
  old_t = 0.0_c_double
  external_birth = 0.0_c_double
  escape = 0.0_c_double
  nearest = nearest_cell()
  old_s(nearest, 2) = 1.0e15_c_double
  old_t(nearest, 2) = 2.0e15_c_double
  electron_density = 1.0e19_c_double
  electron_energy = 1.5_c_double * electron_density * 5.0e3_c_double * ev_j
  ion_energy = 1.5_c_double * thermal_number(3) * 10.0e3_c_double * ev_j
  call make_options()
  thermal_tables = c_null_ptr
  entries = beam_entry_type(0_c_int, 0_c_int, 0_c_int, c_null_ptr)
  bad_entries = beam_entry_type(0_c_int, 0_c_int, 0_c_int, c_null_ptr)
  table = c_null_ptr

  ! The no-table/no-entry source wrapper is the same direct-fast arithmetic.
  call fast_trial(dt_s, options, fast_options, cells, edges, thermal_number, &
       electron_energy, ion_energy, electron_density, thermal_charge_squared, &
       0_c_int, inert, coulomb_logs, old_s, old_t, external_birth, escape, &
       trial_number_fast, trial_s_fast, trial_t_fast, output_fast, status_fast)
  call source_trial(dt_s, options, fast_options, cells, edges, thermal_number, &
       electron_energy, ion_energy, electron_density, thermal_charge_squared, &
       0_c_int, inert, coulomb_logs, old_s, old_t, external_birth, escape, &
       trial_number_source, trial_s_source, trial_t_source, output_source, &
       diagnostics_source, usage_source, status_source)
  call check(status_fast == source_ok .and. status_source == source_ok, &
       'empty beam source wrapper and coupled-fast call succeed')
  call check(all(trial_number_fast == trial_number_source) .and. &
       all(trial_s_fast == trial_s_source) .and. all(trial_t_fast == trial_t_source) .and. &
       coupled_outputs_equal(output_fast, output_source), &
       'empty beam source wrapper preserves direct-fast outputs')
  call check(all(diagnostics_source%tested == 0_c_int) .and. &
       usage_source%direct_evaluations > 0_c_int .and. &
       usage_source%table_evaluations == 0_c_int, &
       'empty beam source reports direct evaluations and cleared diagnostics')

  call make_beam_source(beam_source)
  call make_beam_control(beam_control)
  call fusion_beam_birth_table_create(table_channel, table_slot, &
       0.5_c_double * (edges(nearest) + edges(nearest + 1)), 8.0_c_double * kev_j, &
       12.0_c_double * kev_j, beam_source, beam_control, edges, table, &
       status_table_create)
  call check(status_table_create == beam_table_ok, &
       'real fixed-energy beam table construction succeeds')
  if (status_table_create == beam_table_ok) then
     entries(1)%channel = table_channel
     entries(1)%projectile_slot = table_slot
     entries(1)%energy_cell = int(nearest - 1, c_int)
     entries(1)%table = table
     call source_trial(dt_s, options, fast_options, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, thermal_charge_squared, &
          0_c_int, inert, coulomb_logs, old_s, old_t, external_birth, escape, &
          trial_number_table, trial_s_table, trial_t_table, output_table, &
          diagnostics_table, usage_table, status_table, thermal_tables=thermal_tables, &
          beam_entries=entries, effective_charge=0_c_int)
     call check(status_table == source_ok .and. usage_table%table_evaluations > 0_c_int, &
          'created beam table dispatch succeeds through the source wrapper')

     bad_entries = entries
     bad_entries(1)%energy_cell = entries(1)%energy_cell + 1_c_int
     trial_number_bad = -7.0_c_double
     trial_s_bad = -7.0_c_double
     trial_t_bad = -7.0_c_double
     output_bad%electron_energy_J_m3 = -7.0_c_double
     diagnostics_bad%tested = 1_c_int
     usage_bad%table_evaluations = -7_c_int
     call source_trial(dt_s, options, fast_options, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, thermal_charge_squared, &
          0_c_int, inert, coulomb_logs, old_s, old_t, external_birth, escape, &
          trial_number_bad, trial_s_bad, trial_t_bad, output_bad, diagnostics_bad, &
          usage_bad, status_bad, beam_entries=bad_entries, effective_charge=0_c_int)
     call check(status_bad /= source_ok .and. source_outputs_zero(trial_number_bad, &
          trial_s_bad, trial_t_bad, output_bad, diagnostics_bad, usage_bad), &
          'mismatched beam entry rejects and clears every output')
     call fusion_beam_birth_table_destroy(table)
  end if

  trial_number_shape = -7.0_c_double
  trial_s_shape = -7.0_c_double
  trial_t_shape = -7.0_c_double
  output_shape%electron_energy_J_m3 = -7.0_c_double
  diagnostics_shape%tested = 1_c_int
  usage_shape%table_evaluations = -7_c_int
  call source_trial(dt_s, options, fast_options, cells, edges_bad, thermal_number, &
       electron_energy, ion_energy, electron_density, thermal_charge_squared, &
       0_c_int, inert, coulomb_logs, old_s, old_t, external_birth, escape, &
       trial_number_shape, trial_s_shape, trial_t_shape, output_shape, &
       diagnostics_shape, usage_shape, status_shape)
  call check(status_shape == source_invalid .and. source_outputs_zero( &
       trial_number_shape, trial_s_shape, trial_t_shape, output_shape, &
       diagnostics_shape, usage_shape), &
       'invalid edge extent rejects before pointer formation and clears outputs')

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' coupled-source Fortran test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'Fortran coupled-source binding smoke test passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  subroutine make_edges()
    integer :: j

    edges(1) = 0.0_c_double
    edges(2) = first_edge
    ratio = (max_energy / first_edge) ** &
         (1.0_c_double / real(cells - 1_c_int, c_double))
    do j = 3, cells + 1
       edges(j) = first_edge * ratio ** real(j - 2, c_double)
    end do
    edges(cells + 1) = max_energy
  end subroutine make_edges

  integer function nearest_cell()
    integer :: j

    best_distance = huge(1.0_c_double)
    nearest_cell = 1
    do j = 1, cells
       center = 0.5_c_double * (edges(j) + edges(j + 1))
       distance = abs(center - target_energy)
       if (distance < best_distance) then
          best_distance = distance
          nearest_cell = j
       end if
    end do
  end function nearest_cell

  subroutine make_options()
    options%birth%relative_max_J = 5.0e6_c_double * ev_j
    options%birth%cm_max_kT = 40.0_c_double
    options%birth%ground_state_q_J = 91.84e3_c_double * ev_j
    options%birth%cutoff_J = 2.0e3_c_double * ev_j
    options%birth%l1_fraction = 0.76_c_double
    options%birth%relative_phase = 0.0_c_double
    options%birth%narrow_peak_fraction = 0.051_c_double
    options%birth%continuum_peak_scale = 1.0_c_double
    options%birth%continuation = beam_endpoint_s
    options%birth%pb_low = beam_pb_low_tb
    options%birth%remainder_policy = beam_remainder_proxy
    options%birth%broad_mode = 13_c_int
    options%birth%fsci_policy = 0_c_int
    options%birth%relative_order = 8_c_int
    options%birth%cm_order = 8_c_int
    options%birth%nq = 8_c_int
    options%birth%ncos = 8_c_int
    options%max_source_rate_error = 2.0e-3_c_double
    options%max_source_debit_error = 2.0e-3_c_double
    options%handoff_max_L1 = 1.0e-3_c_double
    options%handoff_max_mean_error = 1.0e-3_c_double
    options%channels = 0_c_int
    options%handoff_enabled = 0_c_int
    fast_options%channels = 0_c_int
    fast_options%channels(table_channel + 1_c_int) = 1_c_int
    fast_options%angular_order = 16_c_int
    fast_options%angular_max_exponent = 40.0_c_double
  end subroutine make_options

  subroutine make_beam_source(source)
    type(fusion_beam_birth_options_v1), intent(out) :: source

    source%relative_max_J = options%birth%relative_max_J
    source%angular_max_exponent = fast_options%angular_max_exponent
    source%ground_state_q_J = options%birth%ground_state_q_J
    source%cutoff_J = options%birth%cutoff_J
    source%l1_fraction = options%birth%l1_fraction
    source%relative_phase = options%birth%relative_phase
    source%narrow_peak_fraction = options%birth%narrow_peak_fraction
    source%continuum_peak_scale = options%birth%continuum_peak_scale
    source%continuation = options%birth%continuation
    source%pb_low = options%birth%pb_low
    source%remainder_policy = options%birth%remainder_policy
    source%broad_mode = options%birth%broad_mode
    source%fsci_policy = options%birth%fsci_policy
    source%relative_order = options%birth%relative_order
    source%angular_order = fast_options%angular_order
    source%nq = options%birth%nq
    source%ncos = options%birth%ncos
  end subroutine make_beam_source

  subroutine make_beam_control(control)
    type(fusion_birth_table_control_v1), intent(out) :: control

    control%max_rate_error = 0.01_c_double
    control%max_debit_error = 0.01_c_double
    control%max_number_L1 = 0.03_c_double
    control%max_energy_L1 = 0.03_c_double
    control%max_direct_rate_discrepancy = 1.0e-5_c_double
    control%max_direct_debit_discrepancy = 1.0e-5_c_double
    control%max_knots = 256_c_int
    control%max_evaluations = 2048_c_int
    control%max_depth = 16_c_int
  end subroutine make_beam_control

  logical function coupled_outputs_equal(left, right)
    type(coupled_output_type), intent(in) :: left, right

    coupled_outputs_equal = all(left%ledger%events_m3 == right%ledger%events_m3) .and. &
         all(left%ledger%nuclear_born_number_m3 == right%ledger%nuclear_born_number_m3) .and. &
         all(left%ledger%nuclear_born_energy_J_m3 == right%ledger%nuclear_born_energy_J_m3) .and. &
         all(left%ledger%external_born_number_m3 == right%ledger%external_born_number_m3) .and. &
         all(left%ledger%external_born_energy_J_m3 == right%ledger%external_born_energy_J_m3) .and. &
         all(left%ledger%thermal_consumed_number_m3 == right%ledger%thermal_consumed_number_m3) .and. &
         all(left%ledger%thermal_consumed_energy_J_m3 == right%ledger%thermal_consumed_energy_J_m3) .and. &
         all(left%ledger%fast_consumed_number_m3 == right%ledger%fast_consumed_number_m3) .and. &
         all(left%ledger%fast_consumed_energy_J_m3 == right%ledger%fast_consumed_energy_J_m3) .and. &
         all(left%ledger%escaped_number_m3 == right%ledger%escaped_number_m3) .and. &
         all(left%ledger%escaped_energy_J_m3 == right%ledger%escaped_energy_J_m3) .and. &
         all(left%ledger%handed_off_number_m3 == right%ledger%handed_off_number_m3) .and. &
         all(left%ledger%handed_off_energy_J_m3 == right%ledger%handed_off_energy_J_m3) .and. &
         all(left%ledger%heat_to_bath_J_m3 == right%ledger%heat_to_bath_J_m3) .and. &
         left%ledger%neutron_number_m3 == right%ledger%neutron_number_m3 .and. &
         left%ledger%neutron_energy_J_m3 == right%ledger%neutron_energy_J_m3 .and. &
         all(left%inert_ion_heat_J_m3 == right%inert_ion_heat_J_m3) .and. &
         left%electron_energy_J_m3 == right%electron_energy_J_m3 .and. &
         left%ion_energy_J_m3 == right%ion_energy_J_m3 .and. &
         all(left%particle_residual_m3 == right%particle_residual_m3) .and. &
         left%energy_residual_J_m3 == right%energy_residual_J_m3 .and. &
         left%max_source_rate_discrepancy == right%max_source_rate_discrepancy .and. &
         left%max_source_debit_discrepancy == right%max_source_debit_discrepancy .and. &
         all(left%handoff_L1 == right%handoff_L1) .and. &
         all(left%handoff_mean_error == right%handoff_mean_error) .and. &
         all(left%handoff_projected == right%handoff_projected)
  end function coupled_outputs_equal

  logical function source_outputs_zero(numbers, s, t, out, diagnostics, usage)
    real(c_double), intent(in) :: numbers(:), s(:,:), t(:,:)
    type(coupled_output_type), intent(in) :: out
    type(diagnostics_type), intent(in) :: diagnostics
    type(usage_type), intent(in) :: usage

    source_outputs_zero = all(numbers == 0.0_c_double) .and. &
         all(s == 0.0_c_double) .and. all(t == 0.0_c_double) .and. &
         all(out%ledger%events_m3 == 0.0_c_double) .and. &
         all(out%ledger%nuclear_born_number_m3 == 0.0_c_double) .and. &
         all(out%ledger%nuclear_born_energy_J_m3 == 0.0_c_double) .and. &
         all(out%ledger%external_born_number_m3 == 0.0_c_double) .and. &
         all(out%ledger%external_born_energy_J_m3 == 0.0_c_double) .and. &
         all(out%ledger%thermal_consumed_number_m3 == 0.0_c_double) .and. &
         all(out%ledger%thermal_consumed_energy_J_m3 == 0.0_c_double) .and. &
         all(out%ledger%fast_consumed_number_m3 == 0.0_c_double) .and. &
         all(out%ledger%fast_consumed_energy_J_m3 == 0.0_c_double) .and. &
         all(out%ledger%escaped_number_m3 == 0.0_c_double) .and. &
         all(out%ledger%escaped_energy_J_m3 == 0.0_c_double) .and. &
         all(out%ledger%handed_off_number_m3 == 0.0_c_double) .and. &
         all(out%ledger%handed_off_energy_J_m3 == 0.0_c_double) .and. &
         all(out%ledger%heat_to_bath_J_m3 == 0.0_c_double) .and. &
         out%ledger%neutron_number_m3 == 0.0_c_double .and. &
         out%ledger%neutron_energy_J_m3 == 0.0_c_double .and. &
         all(out%inert_ion_heat_J_m3 == 0.0_c_double) .and. &
         out%electron_energy_J_m3 == 0.0_c_double .and. &
         out%ion_energy_J_m3 == 0.0_c_double .and. &
         all(out%particle_residual_m3 == 0.0_c_double) .and. &
         out%energy_residual_J_m3 == 0.0_c_double .and. &
         out%max_source_rate_discrepancy == 0.0_c_double .and. &
         out%max_source_debit_discrepancy == 0.0_c_double .and. &
         all(out%handoff_L1 == 0.0_c_double) .and. &
         all(out%handoff_mean_error == 0.0_c_double) .and. &
         all(out%handoff_projected == 0_c_int) .and. &
         all(diagnostics%candidate_number_m3 == 0.0_c_double) .and. &
         all(diagnostics%candidate_energy_J_m3 == 0.0_c_double) .and. &
         all(diagnostics%transferred_number_m3 == 0.0_c_double) .and. &
         all(diagnostics%transferred_energy_J_m3 == 0.0_c_double) .and. &
         all(diagnostics%target_kT_J == 0.0_c_double) .and. &
         all(diagnostics%tested == 0_c_int) .and. &
         usage%direct_evaluations == 0_c_int .and. &
         usage%table_evaluations == 0_c_int .and. &
         usage%max_validated_rate_error == 0.0_c_double .and. &
         usage%max_validated_debit_error == 0.0_c_double .and. &
         usage%max_validated_number_L1 == 0.0_c_double .and. &
         usage%max_validated_energy_L1 == 0.0_c_double
  end function source_outputs_zero

end program test_fusion_coupled_sources_fortran
