module fusion_coupled_sources_fortran
  !! ISO_C_BINDING wrapper for the additive coupled source trial.
  !!
  !! This wrapper reuses the coupled-fast ABI types and shape contract.  It
  !! only translates Fortran array descriptors to borrowed C pointers; table
  !! construction, ownership, interpolation, and all source physics remain in
  !! the C++ library.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_ptr, &
       c_null_ptr, c_loc
  use fusion_coupled_fast_fortran, only : &
       fusion_thermal_birth_options_v1, fusion_source_ledger_v1, &
       fusion_inert_ion_v1, fusion_fast_target_options_v1, &
       fusion_coupled_thermal_options_v1, fusion_coupled_thermal_v1, &
       fusion_handoff_diagnostics_v1, FUSION_COUPLED_THERMAL_SPECIES, &
       FUSION_COUPLED_THERMAL_CHANNELS, FUSION_COUPLED_THERMAL_MAX_CELLS, &
       FUSION_COUPLED_THERMAL_MAX_INERT, FUSION_COUPLED_THERMAL_BASE_BATHS, &
       FUSION_COUPLED_FAST_CHANNELS, PB11_STATUS_OK, &
       PB11_STATUS_NULL_OUTPUT, PB11_STATUS_INVALID_ARGUMENT, &
       PB11_STATUS_OUT_OF_RANGE, PB11_STATUS_NUMERICAL_FAILURE, &
       PB11_STATUS_EXCEPTION, PB11_STATUS_UNKNOWN_METHOD
  implicit none
  private

  ! Re-export the common ABI types and constants from the coupled-fast module.
  public :: fusion_thermal_birth_options_v1
  public :: fusion_source_ledger_v1
  public :: fusion_inert_ion_v1
  public :: fusion_fast_target_options_v1
  public :: fusion_coupled_thermal_options_v1
  public :: fusion_coupled_thermal_v1
  public :: fusion_handoff_diagnostics_v1
  public :: FUSION_COUPLED_THERMAL_SPECIES
  public :: FUSION_COUPLED_THERMAL_CHANNELS
  public :: FUSION_COUPLED_THERMAL_MAX_CELLS
  public :: FUSION_COUPLED_THERMAL_MAX_INERT
  public :: FUSION_COUPLED_THERMAL_BASE_BATHS
  public :: FUSION_COUPLED_FAST_CHANNELS
  public :: PB11_STATUS_OK
  public :: PB11_STATUS_NULL_OUTPUT
  public :: PB11_STATUS_INVALID_ARGUMENT
  public :: PB11_STATUS_OUT_OF_RANGE
  public :: PB11_STATUS_NUMERICAL_FAILURE
  public :: PB11_STATUS_EXCEPTION
  public :: PB11_STATUS_UNKNOWN_METHOD

  ! Exact C layout: three c_int values followed by one borrowed table handle.
  type, bind(C), public :: fusion_beam_table_entry_v1
     integer(c_int) :: channel
     integer(c_int) :: projectile_slot
     integer(c_int) :: energy_cell
     type(c_ptr) :: table
  end type fusion_beam_table_entry_v1

  ! Exact C layout: two c_int values followed by four c_double values.
  type, bind(C), public :: fusion_beam_table_usage_v1
     integer(c_int) :: direct_evaluations
     integer(c_int) :: table_evaluations
     real(c_double) :: max_validated_rate_error
     real(c_double) :: max_validated_debit_error
     real(c_double) :: max_validated_number_L1
     real(c_double) :: max_validated_energy_L1
  end type fusion_beam_table_usage_v1

  public :: fusion_coupled_sources_trial

  interface
     function c_fusion_c_coupled_sources_trial( &
          dt_s, options, fast_options, thermal_tables, beam_table_count, &
          beam_tables, effective_charge, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out, diagnostics, usage) &
          bind(C, name="fusion_c_coupled_sources_trial") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
       type(c_ptr), value :: thermal_tables
       integer(c_int), value :: beam_table_count
       type(c_ptr), value :: beam_tables
       integer(c_int), value :: effective_charge
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: thermal_number
       real(c_double), value :: electron_energy
       real(c_double), value :: ion_energy
       real(c_double), value :: electron_density
       type(c_ptr), value :: thermal_charge_squared
       integer(c_int), value :: inert_count
       type(c_ptr), value :: inert
       type(c_ptr), value :: coulomb_logs
       type(c_ptr), value :: old_s
       type(c_ptr), value :: old_t
       type(c_ptr), value :: external_birth
       type(c_ptr), value :: escape
       type(c_ptr), value :: trial_thermal_number
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       type(c_ptr), value :: out
       type(c_ptr), value :: diagnostics
       type(c_ptr), value :: usage
       integer(c_int) :: status
     end function c_fusion_c_coupled_sources_trial
  end interface

contains

  subroutine clear_source_ledger(ledger)
    type(fusion_source_ledger_v1), intent(out) :: ledger

    ledger%events_m3 = 0.0_c_double
    ledger%nuclear_born_number_m3 = 0.0_c_double
    ledger%nuclear_born_energy_J_m3 = 0.0_c_double
    ledger%external_born_number_m3 = 0.0_c_double
    ledger%external_born_energy_J_m3 = 0.0_c_double
    ledger%thermal_consumed_number_m3 = 0.0_c_double
    ledger%thermal_consumed_energy_J_m3 = 0.0_c_double
    ledger%fast_consumed_number_m3 = 0.0_c_double
    ledger%fast_consumed_energy_J_m3 = 0.0_c_double
    ledger%escaped_number_m3 = 0.0_c_double
    ledger%escaped_energy_J_m3 = 0.0_c_double
    ledger%handed_off_number_m3 = 0.0_c_double
    ledger%handed_off_energy_J_m3 = 0.0_c_double
    ledger%heat_to_bath_J_m3 = 0.0_c_double
    ledger%neutron_number_m3 = 0.0_c_double
    ledger%neutron_energy_J_m3 = 0.0_c_double
  end subroutine clear_source_ledger

  subroutine clear_source_output(out)
    type(fusion_coupled_thermal_v1), intent(out) :: out

    call clear_source_ledger(out%ledger)
    out%inert_ion_heat_J_m3 = 0.0_c_double
    out%electron_energy_J_m3 = 0.0_c_double
    out%ion_energy_J_m3 = 0.0_c_double
    out%particle_residual_m3 = 0.0_c_double
    out%energy_residual_J_m3 = 0.0_c_double
    out%max_source_rate_discrepancy = 0.0_c_double
    out%max_source_debit_discrepancy = 0.0_c_double
    out%handoff_L1 = 0.0_c_double
    out%handoff_mean_error = 0.0_c_double
    out%handoff_projected = 0_c_int
  end subroutine clear_source_output

  subroutine clear_source_diagnostics(diagnostics)
    type(fusion_handoff_diagnostics_v1), intent(out) :: diagnostics

    diagnostics%candidate_number_m3 = 0.0_c_double
    diagnostics%candidate_energy_J_m3 = 0.0_c_double
    diagnostics%transferred_number_m3 = 0.0_c_double
    diagnostics%transferred_energy_J_m3 = 0.0_c_double
    diagnostics%target_kT_J = 0.0_c_double
    diagnostics%tested = 0_c_int
  end subroutine clear_source_diagnostics

  subroutine clear_source_usage(usage)
    type(fusion_beam_table_usage_v1), intent(out) :: usage

    usage%direct_evaluations = 0_c_int
    usage%table_evaluations = 0_c_int
    usage%max_validated_rate_error = 0.0_c_double
    usage%max_validated_debit_error = 0.0_c_double
    usage%max_validated_number_L1 = 0.0_c_double
    usage%max_validated_energy_L1 = 0.0_c_double
  end subroutine clear_source_usage

  subroutine clear_source_outputs(trial_thermal_number_m3, trial_s_m3, &
       trial_t_m3, out, diagnostics, usage)
    real(c_double), intent(out) :: trial_thermal_number_m3(:)
    real(c_double), intent(out) :: trial_s_m3(:,:)
    real(c_double), intent(out) :: trial_t_m3(:,:)
    type(fusion_coupled_thermal_v1), intent(out) :: out
    type(fusion_handoff_diagnostics_v1), intent(out) :: diagnostics
    type(fusion_beam_table_usage_v1), intent(out) :: usage

    trial_thermal_number_m3 = 0.0_c_double
    trial_s_m3 = 0.0_c_double
    trial_t_m3 = 0.0_c_double
    call clear_source_output(out)
    call clear_source_diagnostics(diagnostics)
    call clear_source_usage(usage)
  end subroutine clear_source_outputs

  logical function valid_source_cells(cells)
    integer(c_int), intent(in) :: cells

    valid_source_cells = cells >= 1_c_int .and. &
         cells <= FUSION_COUPLED_THERMAL_MAX_CELLS
  end function valid_source_cells

  subroutine prepare_source_arguments(options, fast_options, cells, edges_J, &
       thermal_number_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, thermal_tables, &
       beam_entries, effective_charge, options_ptr, fast_options_ptr, &
       thermal_tables_ptr, beam_table_count, beam_tables_ptr, effective_charge_value, &
       edges_ptr, thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr, &
       coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, escape_ptr, &
       trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, status)
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in), target, contiguous :: thermal_charge_squared(:)
    integer(c_int), intent(in) :: inert_count
    type(fusion_inert_ion_v1), intent(in), target, contiguous :: inert(:)
    real(c_double), intent(in), target, contiguous :: coulomb_logs(:,:)
    real(c_double), intent(in), target, contiguous :: old_s_m3(:,:)
    real(c_double), intent(in), target, contiguous :: old_t_m3(:,:)
    real(c_double), intent(in), target, contiguous :: external_birth_m3_s(:,:)
    real(c_double), intent(in), target, contiguous :: escape_s_inv(:,:)
    real(c_double), intent(out), target, contiguous :: trial_thermal_number_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_s_m3(:,:)
    real(c_double), intent(out), target, contiguous :: trial_t_m3(:,:)
    type(c_ptr), intent(in), target, contiguous, optional :: thermal_tables(:)
    type(fusion_beam_table_entry_v1), intent(in), target, contiguous, optional :: &
         beam_entries(:)
    integer(c_int), intent(in), optional :: effective_charge
    type(c_ptr), intent(out) :: options_ptr, fast_options_ptr
    type(c_ptr), intent(out) :: thermal_tables_ptr, beam_tables_ptr
    integer(c_int), intent(out) :: beam_table_count, effective_charge_value
    type(c_ptr), intent(out) :: edges_ptr, thermal_number_ptr
    type(c_ptr), intent(out) :: thermal_charge_squared_ptr, inert_ptr
    type(c_ptr), intent(out) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr), intent(out) :: external_birth_ptr, escape_ptr
    type(c_ptr), intent(out) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: cells_size, inert_size, expected_baths
    integer(c_size_t) :: expected_edges, expected_logs, species_size
    integer(c_size_t) :: max_beam_entries, beam_size

    options_ptr = c_null_ptr
    fast_options_ptr = c_null_ptr
    thermal_tables_ptr = c_null_ptr
    beam_tables_ptr = c_null_ptr
    edges_ptr = c_null_ptr
    thermal_number_ptr = c_null_ptr
    thermal_charge_squared_ptr = c_null_ptr
    inert_ptr = c_null_ptr
    coulomb_logs_ptr = c_null_ptr
    old_s_ptr = c_null_ptr
    old_t_ptr = c_null_ptr
    external_birth_ptr = c_null_ptr
    escape_ptr = c_null_ptr
    trial_thermal_number_ptr = c_null_ptr
    trial_s_ptr = c_null_ptr
    trial_t_ptr = c_null_ptr
    beam_table_count = 0_c_int
    effective_charge_value = 0_c_int
    status = PB11_STATUS_INVALID_ARGUMENT

    if (present(effective_charge)) effective_charge_value = effective_charge
    if (effective_charge_value /= 0_c_int .and. &
         effective_charge_value /= 1_c_int) return

    ! Check signed values before converting to c_size_t.  No C_LOC is formed
    ! until every required extent, including optional arrays, has passed.
    if (.not. valid_source_cells(cells)) return
    if (inert_count < 0_c_int .or. &
         inert_count > FUSION_COUPLED_THERMAL_MAX_INERT) return
    cells_size = int(cells, kind=c_size_t)
    inert_size = int(inert_count, kind=c_size_t)
    species_size = int(FUSION_COUPLED_THERMAL_SPECIES, kind=c_size_t)
    expected_edges = cells_size + 1_c_size_t
    expected_baths = int(FUSION_COUPLED_THERMAL_BASE_BATHS, kind=c_size_t) + &
         inert_size
    expected_logs = expected_baths * species_size
    max_beam_entries = 10_c_size_t * cells_size

    if (size(edges_J, kind=c_size_t) /= expected_edges) return
    if (size(thermal_number_m3, kind=c_size_t) /= species_size) return
    if (size(thermal_charge_squared, kind=c_size_t) /= species_size) return
    if (size(trial_thermal_number_m3, kind=c_size_t) /= species_size) return
    if (size(inert, kind=c_size_t) /= inert_size) return
    if (size(coulomb_logs, 1, kind=c_size_t) /= expected_baths) return
    if (size(coulomb_logs, 2, kind=c_size_t) /= species_size) return
    if (size(coulomb_logs, kind=c_size_t) /= expected_logs) return
    if (size(old_s_m3, 1, kind=c_size_t) /= cells_size) return
    if (size(old_s_m3, 2, kind=c_size_t) /= species_size) return
    if (size(old_t_m3, 1, kind=c_size_t) /= cells_size) return
    if (size(old_t_m3, 2, kind=c_size_t) /= species_size) return
    if (size(external_birth_m3_s, 1, kind=c_size_t) /= cells_size) return
    if (size(external_birth_m3_s, 2, kind=c_size_t) /= species_size) return
    if (size(escape_s_inv, 1, kind=c_size_t) /= cells_size) return
    if (size(escape_s_inv, 2, kind=c_size_t) /= species_size) return
    if (size(trial_s_m3, 1, kind=c_size_t) /= cells_size) return
    if (size(trial_s_m3, 2, kind=c_size_t) /= species_size) return
    if (size(trial_t_m3, 1, kind=c_size_t) /= cells_size) return
    if (size(trial_t_m3, 2, kind=c_size_t) /= species_size) return

    if (present(thermal_tables)) then
       if (size(thermal_tables, kind=c_size_t) /= &
            int(FUSION_COUPLED_THERMAL_CHANNELS, kind=c_size_t)) return
    end if
    beam_size = 0_c_size_t
    if (present(beam_entries)) beam_size = size(beam_entries, kind=c_size_t)
    if (beam_size > max_beam_entries) return
    beam_table_count = int(beam_size, kind=c_int)

    options_ptr = c_loc(options)
    fast_options_ptr = c_loc(fast_options)
    edges_ptr = c_loc(edges_J(1))
    thermal_number_ptr = c_loc(thermal_number_m3(1))
    thermal_charge_squared_ptr = c_loc(thermal_charge_squared(1))
    coulomb_logs_ptr = c_loc(coulomb_logs(1,1))
    old_s_ptr = c_loc(old_s_m3(1,1))
    old_t_ptr = c_loc(old_t_m3(1,1))
    external_birth_ptr = c_loc(external_birth_m3_s(1,1))
    escape_ptr = c_loc(escape_s_inv(1,1))
    trial_thermal_number_ptr = c_loc(trial_thermal_number_m3(1))
    trial_s_ptr = c_loc(trial_s_m3(1,1))
    trial_t_ptr = c_loc(trial_t_m3(1,1))
    if (present(thermal_tables)) thermal_tables_ptr = c_loc(thermal_tables(1))
    if (beam_size > 0_c_size_t) beam_tables_ptr = c_loc(beam_entries(1))
    if (beam_size == 0_c_size_t) beam_tables_ptr = c_null_ptr
    status = PB11_STATUS_OK
  end subroutine prepare_source_arguments

  subroutine fusion_coupled_sources_trial(dt_s, options, fast_options, cells, &
       edges_J, thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, diagnostics, &
       usage, status, thermal_tables, beam_entries, effective_charge)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in) :: electron_energy_J_m3
    real(c_double), intent(in) :: ion_energy_J_m3
    real(c_double), intent(in) :: electron_density_m3
    real(c_double), intent(in), target, contiguous :: thermal_charge_squared(:)
    integer(c_int), intent(in) :: inert_count
    type(fusion_inert_ion_v1), intent(in), target, contiguous :: inert(:)
    real(c_double), intent(in), target, contiguous :: coulomb_logs(:,:)
    real(c_double), intent(in), target, contiguous :: old_s_m3(:,:)
    real(c_double), intent(in), target, contiguous :: old_t_m3(:,:)
    real(c_double), intent(in), target, contiguous :: external_birth_m3_s(:,:)
    real(c_double), intent(in), target, contiguous :: escape_s_inv(:,:)
    real(c_double), intent(out), target, contiguous :: &
         trial_thermal_number_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_s_m3(:,:)
    real(c_double), intent(out), target, contiguous :: trial_t_m3(:,:)
    type(fusion_coupled_thermal_v1), intent(out), target :: out
    type(fusion_handoff_diagnostics_v1), intent(out), target :: diagnostics
    type(fusion_beam_table_usage_v1), intent(out), target :: usage
    integer(c_int), intent(out) :: status
    type(c_ptr), intent(in), target, contiguous, optional :: thermal_tables(:)
    type(fusion_beam_table_entry_v1), intent(in), target, contiguous, optional :: &
         beam_entries(:)
    integer(c_int), intent(in), optional :: effective_charge

    type(c_ptr) :: options_ptr, fast_options_ptr, thermal_tables_ptr
    type(c_ptr) :: beam_tables_ptr, edges_ptr, thermal_number_ptr
    type(c_ptr) :: thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr
    type(c_ptr) :: old_s_ptr, old_t_ptr, external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    integer(c_int) :: beam_table_count, effective_charge_value

    call clear_source_outputs(trial_thermal_number_m3, trial_s_m3, trial_t_m3, &
         out, diagnostics, usage)
    call prepare_source_arguments(options, fast_options, cells, edges_J, &
         thermal_number_m3, thermal_charge_squared, inert_count, inert, &
         coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
         trial_thermal_number_m3, trial_s_m3, trial_t_m3, thermal_tables, &
         beam_entries, effective_charge, options_ptr, fast_options_ptr, &
         thermal_tables_ptr, beam_table_count, beam_tables_ptr, &
         effective_charge_value, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, status)
    if (status /= PB11_STATUS_OK) return

    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_fusion_c_coupled_sources_trial(dt_s, options_ptr, &
         fast_options_ptr, thermal_tables_ptr, beam_table_count, &
         beam_tables_ptr, effective_charge_value, cells, edges_ptr, &
         thermal_number_ptr, electron_energy_J_m3, ion_energy_J_m3, &
         electron_density_m3, thermal_charge_squared_ptr, inert_count, &
         inert_ptr, coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         c_loc(out), c_loc(diagnostics), c_loc(usage))
    if (status /= PB11_STATUS_OK) then
       call clear_source_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out, diagnostics, usage)
    end if
  end subroutine fusion_coupled_sources_trial

end module fusion_coupled_sources_fortran
