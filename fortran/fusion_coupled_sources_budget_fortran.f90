module fusion_coupled_sources_budget_fortran
  !! Explicit highest packet entry; source policy and physical checks belong to C.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_ptr, &
       c_null_ptr, c_loc, c_int64_t, c_int8_t, c_f_pointer, c_sizeof
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
  use fusion_coupled_sources_fortran, only : fusion_beam_table_entry_v1, fusion_beam_table_usage_v1
  use fusion_coupled_sources_floor_fortran, only : fusion_coupled_floor_limits_v1, &
       fusion_birth_packets_v1, fusion_birth_floor_ledger_v1, &
       FUSION_BEAM_TABLE_STRICT, FUSION_BEAM_TABLE_DIRECT_OUTSIDE
  use fusion_birth_table_budget_fortran, only : fusion_birth_table_budget_info_v1, &
       fusion_thermal_birth_budget_policy_v1, fusion_thermal_birth_budget_report_v1
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


  public :: fusion_beam_table_entry_v1, fusion_beam_table_usage_v1
  public :: fusion_coupled_floor_limits_v1, fusion_birth_floor_ledger_v1, fusion_birth_packets_v1
  public :: FUSION_BEAM_TABLE_STRICT, FUSION_BEAM_TABLE_DIRECT_OUTSIDE
  public :: fusion_thermal_birth_budget_policy_v1, fusion_thermal_birth_budget_report_v1
  public :: fusion_birth_table_budget_info_v1, fusion_coupled_sources_packets_budget_trial
  type, bind(C), public :: fusion_coupled_source_accuracy_v1
    integer(c_int) :: version, mode
    real(c_double) :: number_L1_estimate, energy_L1_estimate
    type(fusion_thermal_birth_budget_report_v1) :: direct
    type(fusion_birth_table_budget_info_v1) :: table
  end type fusion_coupled_source_accuracy_v1
  public :: c_fusion_c_coupled_sources_covered_budget_trial
  interface
     function c_fusion_c_coupled_sources_covered_budget_trial( &
          dt_s, options, fast_options, thermal_tables, beam_table_count, &
          beam_tables, effective_charge, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, new_thermal_number, new_s, new_t, &
          out, diagnostics, usage, floor_limits, floor_ledger, table_domain_policy, outside_count, &
          pb_policy, accuracy) &
          bind(C, name="fusion_c_coupled_sources_covered_budget_trial") result(status)
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
       type(c_ptr), value :: new_thermal_number
       type(c_ptr), value :: new_s
       type(c_ptr), value :: new_t
       type(c_ptr), value :: out
       type(c_ptr), value :: diagnostics
       type(c_ptr), value :: usage
       type(c_ptr), value :: floor_limits
       type(c_ptr), value :: floor_ledger
       integer(c_int), value :: table_domain_policy
       type(c_ptr), value :: outside_count
       type(c_ptr), value :: pb_policy, accuracy
       integer(c_int) :: status
     end function c_fusion_c_coupled_sources_covered_budget_trial
     function c_budget_trial( &
          dt_s, options, fast_options, thermal_tables, beam_table_count, &
          beam_tables, effective_charge, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, new_thermal_number, new_s, new_t, &
          out, diagnostics, usage, floor_limits, floor_ledger, table_domain_policy, outside_count, &
          mapped_packets, packets, pb_policy, accuracy) &
          bind(C, name="fusion_c_coupled_sources_packets_budget_trial") result(status)
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
       type(c_ptr), value :: new_thermal_number
       type(c_ptr), value :: new_s
       type(c_ptr), value :: new_t
       type(c_ptr), value :: out
       type(c_ptr), value :: diagnostics
       type(c_ptr), value :: usage
       type(c_ptr), value :: floor_limits
       type(c_ptr), value :: floor_ledger
       integer(c_int), value :: table_domain_policy
       type(c_ptr), value :: outside_count
       type(c_ptr), value :: mapped_packets, packets, pb_policy, accuracy
       integer(c_int) :: status
     end function c_budget_trial
  end interface
contains
  subroutine clear_bytes(address, count)
    type(c_ptr), intent(in) :: address
    integer(c_size_t), intent(in) :: count
    integer(c_int8_t), pointer :: bytes(:)
    call c_f_pointer(address, bytes, [count])
    bytes = 0_c_int8_t
  end subroutine clear_bytes

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

  subroutine fusion_coupled_sources_packets_budget_trial(dt_s, options, fast_options, cells, &
       edges_J, thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, diagnostics, &
       usage, floor_limits, floor_ledger, table_domain_policy, outside_count, &
       mapped_packets, packets, pb_policy, accuracy, status, &
       thermal_tables, beam_entries, effective_charge)
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
    type(fusion_coupled_floor_limits_v1), intent(in), target :: floor_limits
    type(fusion_birth_floor_ledger_v1), intent(out), target :: floor_ledger
    integer(c_int), intent(in) :: table_domain_policy
    integer(c_int64_t), intent(out), target :: outside_count
    ! C 2*5*7*cells includes neutron slot 7; thermal populations remain six.
    ! C [source][channel][species][cell] maps to (cell,species,channel,source).
    real(c_double), intent(out), target, contiguous :: mapped_packets(:,:,:,:)
    type(fusion_birth_packets_v1), intent(out), target :: packets
    type(fusion_thermal_birth_budget_policy_v1), intent(in), target :: pb_policy
    type(fusion_coupled_source_accuracy_v1), intent(out), target :: accuracy
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

    call clear_outputs()
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. valid_source_cells(cells)) return
    if (size(mapped_packets,1,kind=c_size_t) /= int(cells,c_size_t)) return
    if (size(mapped_packets,2,kind=c_size_t) /= 7_c_size_t) return
    if (size(mapped_packets,3,kind=c_size_t) /= 5_c_size_t) return
    if (size(mapped_packets,4,kind=c_size_t) /= 2_c_size_t) return
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
    if (status /= PB11_STATUS_OK) then
       call clear_outputs()
       return
    end if

    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_budget_trial(dt_s, options_ptr, &
         fast_options_ptr, thermal_tables_ptr, beam_table_count, &
         beam_tables_ptr, effective_charge_value, cells, edges_ptr, &
         thermal_number_ptr, electron_energy_J_m3, ion_energy_J_m3, &
         electron_density_m3, thermal_charge_squared_ptr, inert_count, &
         inert_ptr, coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         c_loc(out), c_loc(diagnostics), c_loc(usage), c_loc(floor_limits), &
         c_loc(floor_ledger), table_domain_policy, c_loc(outside_count), &
         c_loc(mapped_packets(1,1,1,1)), c_loc(packets), c_loc(pb_policy), c_loc(accuracy))
    if (status /= PB11_STATUS_OK) then
       call clear_outputs()
    end if
  contains
    subroutine clear_outputs()
      trial_thermal_number_m3 = 0.0_c_double
      trial_s_m3 = 0.0_c_double
      trial_t_m3 = 0.0_c_double
      mapped_packets = 0.0_c_double
      outside_count = 0_c_int64_t
      call clear_bytes(c_loc(out), c_sizeof(out))
      call clear_bytes(c_loc(diagnostics), c_sizeof(diagnostics))
      call clear_bytes(c_loc(usage), c_sizeof(usage))
      call clear_bytes(c_loc(floor_ledger), c_sizeof(floor_ledger))
      call clear_bytes(c_loc(packets), c_sizeof(packets))
      call clear_bytes(c_loc(accuracy), c_sizeof(accuracy))
    end subroutine clear_outputs
  end subroutine fusion_coupled_sources_packets_budget_trial

end module fusion_coupled_sources_budget_fortran
