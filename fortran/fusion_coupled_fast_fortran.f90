module fusion_coupled_fast_fortran
  !! ISO_C_BINDING wrapper for the bounded coupled fast-target trial.
  !!
  !! The thermal options, inert-ion record, source ledger, and coupled output
  !! are imported from fusion_coupled_thermal_fortran.  This module owns only
  !! the additional fast-target options and the three fast entry points.
  !! Fortran arrays use the same (cells,6) shape as the thermal wrapper; the
  !! first dimension is contiguous and therefore presents C's species-major
  !! [species][cell] storage.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_ptr, &
       c_null_ptr, c_loc
  use fusion_thermal_birth_fortran, only : fusion_thermal_birth_options_v1
  use fusion_source_state_fortran, only : fusion_source_ledger_v1
  use fusion_coupled_thermal_fortran, only : &
       fusion_inert_ion_v1, fusion_coupled_thermal_options_v1, &
       fusion_coupled_thermal_v1, &
       fusion_handoff_diagnostics_v1, &
       FUSION_COUPLED_THERMAL_SPECIES, FUSION_COUPLED_THERMAL_CHANNELS, &
       FUSION_COUPLED_THERMAL_MAX_CELLS, FUSION_COUPLED_THERMAL_MAX_INERT, &
       FUSION_COUPLED_THERMAL_BASE_BATHS, PB11_STATUS_OK, &
       PB11_STATUS_NULL_OUTPUT, PB11_STATUS_INVALID_ARGUMENT, &
       PB11_STATUS_OUT_OF_RANGE, PB11_STATUS_NUMERICAL_FAILURE, &
       PB11_STATUS_EXCEPTION, PB11_STATUS_UNKNOWN_METHOD
  implicit none
  private

  ! Re-export the common ABI types and status/shape constants.  Their layouts
  ! remain defined in the existing thermal/source-state bindings.
  public :: fusion_thermal_birth_options_v1
  public :: fusion_source_ledger_v1
  public :: fusion_inert_ion_v1
  public :: fusion_coupled_thermal_options_v1
  public :: fusion_coupled_thermal_v1
  public :: fusion_handoff_diagnostics_v1
  public :: FUSION_COUPLED_THERMAL_SPECIES
  public :: FUSION_COUPLED_THERMAL_CHANNELS
  public :: FUSION_COUPLED_THERMAL_MAX_CELLS
  public :: FUSION_COUPLED_THERMAL_MAX_INERT
  public :: FUSION_COUPLED_THERMAL_BASE_BATHS
  public :: PB11_STATUS_OK
  public :: PB11_STATUS_NULL_OUTPUT
  public :: PB11_STATUS_INVALID_ARGUMENT
  public :: PB11_STATUS_OUT_OF_RANGE
  public :: PB11_STATUS_NUMERICAL_FAILURE
  public :: PB11_STATUS_EXCEPTION
  public :: PB11_STATUS_UNKNOWN_METHOD

  integer(c_int), parameter, public :: FUSION_COUPLED_FAST_CHANNELS = &
       FUSION_COUPLED_THERMAL_CHANNELS

  ! Exact C layout: five c_int channel selectors, one c_int angular order,
  ! then one c_double angular exponent.
  type, bind(C), public :: fusion_fast_target_options_v1
     integer(c_int) :: channels(5)
     integer(c_int) :: angular_order
     real(c_double) :: angular_max_exponent
  end type fusion_fast_target_options_v1

  public :: fusion_coupled_fast_trial
  public :: fusion_coupled_fast_table_trial
  public :: fusion_coupled_fast_table_trial_effective_charge
  public :: fusion_coupled_fast_trial_diagnosed
  public :: fusion_coupled_fast_table_trial_diagnosed
  public :: fusion_coupled_fast_table_trial_effective_charge_diagnosed

  interface
     function c_fusion_c_coupled_fast_trial( &
          dt_s, options, fast_options, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out) bind(C, name="fusion_c_coupled_fast_trial") &
          result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
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
       integer(c_int) :: status
     end function c_fusion_c_coupled_fast_trial

     function c_fusion_c_coupled_fast_table_trial( &
          dt_s, options, fast_options, tables, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out) bind(C, name="fusion_c_coupled_fast_table_trial") &
          result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
       type(c_ptr), value :: tables
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
       integer(c_int) :: status
     end function c_fusion_c_coupled_fast_table_trial

     function c_fusion_c_coupled_fast_table_trial_effective_charge( &
          dt_s, options, fast_options, tables, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out) bind(C, &
          name="fusion_c_coupled_fast_table_trial_effective_charge") &
          result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
       type(c_ptr), value :: tables
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
       integer(c_int) :: status
     end function c_fusion_c_coupled_fast_table_trial_effective_charge

     function c_fusion_c_coupled_fast_trial_diagnosed( &
          dt_s, options, fast_options, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out, diagnostics) bind(C, &
          name="fusion_c_coupled_fast_trial_diagnosed") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
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
       integer(c_int) :: status
     end function c_fusion_c_coupled_fast_trial_diagnosed

     function c_fusion_c_coupled_fast_table_trial_diagnosed( &
          dt_s, options, fast_options, tables, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out, diagnostics) bind(C, &
          name="fusion_c_coupled_fast_table_trial_diagnosed") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
       type(c_ptr), value :: tables
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
       integer(c_int) :: status
     end function c_fusion_c_coupled_fast_table_trial_diagnosed

     function c_fusion_c_coupled_fast_table_trial_effective_charge_diagnosed( &
          dt_s, options, fast_options, tables, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out, diagnostics) bind(C, &
          name="fusion_c_coupled_fast_table_trial_effective_charge_diagnosed") &
          result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
       type(c_ptr), value :: tables
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
       integer(c_int) :: status
     end function c_fusion_c_coupled_fast_table_trial_effective_charge_diagnosed
  end interface

contains

  subroutine clear_fast_ledger(ledger)
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
  end subroutine clear_fast_ledger

  subroutine clear_fast_output(out)
    type(fusion_coupled_thermal_v1), intent(out) :: out

    call clear_fast_ledger(out%ledger)
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
  end subroutine clear_fast_output

  subroutine clear_fast_diagnostics(diagnostics)
    type(fusion_handoff_diagnostics_v1), intent(out) :: diagnostics

    diagnostics%candidate_number_m3 = 0.0_c_double
    diagnostics%candidate_energy_J_m3 = 0.0_c_double
    diagnostics%transferred_number_m3 = 0.0_c_double
    diagnostics%transferred_energy_J_m3 = 0.0_c_double
    diagnostics%target_kT_J = 0.0_c_double
    diagnostics%tested = 0_c_int
  end subroutine clear_fast_diagnostics

  subroutine clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
       trial_t_m3, out)
    real(c_double), intent(out) :: trial_thermal_number_m3(:)
    real(c_double), intent(out) :: trial_s_m3(:,:)
    real(c_double), intent(out) :: trial_t_m3(:,:)
    type(fusion_coupled_thermal_v1), intent(out) :: out

    call clear_fast_output(out)
    trial_thermal_number_m3 = 0.0_c_double
    trial_s_m3 = 0.0_c_double
    trial_t_m3 = 0.0_c_double
  end subroutine clear_fast_outputs

  logical function valid_fast_cells(cells)
    integer(c_int), intent(in) :: cells

    valid_fast_cells = cells >= 1_c_int .and. &
         cells <= FUSION_COUPLED_THERMAL_MAX_CELLS
  end function valid_fast_cells

  subroutine prepare_fast_arguments(options, fast_options, cells, edges_J, &
       thermal_number_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
       fast_options_ptr, edges_ptr, thermal_number_ptr, &
       thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
       old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
       trial_s_ptr, trial_t_ptr, tables_ptr, status, tables)
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in), target, contiguous :: &
         thermal_charge_squared(:)
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
    type(c_ptr), intent(out) :: options_ptr, fast_options_ptr, edges_ptr
    type(c_ptr), intent(out) :: thermal_number_ptr
    type(c_ptr), intent(out) :: thermal_charge_squared_ptr, inert_ptr
    type(c_ptr), intent(out) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr), intent(out) :: external_birth_ptr, escape_ptr
    type(c_ptr), intent(out) :: trial_thermal_number_ptr, trial_s_ptr
    type(c_ptr), intent(out) :: trial_t_ptr, tables_ptr
    integer(c_int), intent(out) :: status
    type(c_ptr), intent(in), target, contiguous, optional :: tables(:)

    integer(c_size_t) :: cells_size, inert_size, expected_baths
    integer(c_size_t) :: expected_edges, expected_logs, species_size
    integer(c_size_t) :: channel_size

    options_ptr = c_null_ptr
    fast_options_ptr = c_null_ptr
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
    tables_ptr = c_null_ptr
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Check signed c_int values before converting them to c_size_t.  Every
    ! subsequent extent comparison stays in c_size_t, avoiding narrowing.
    if (.not. valid_fast_cells(cells)) return
    if (inert_count < 0_c_int .or. &
         inert_count > FUSION_COUPLED_THERMAL_MAX_INERT) return
    cells_size = int(cells, kind=c_size_t)
    inert_size = int(inert_count, kind=c_size_t)
    species_size = int(FUSION_COUPLED_THERMAL_SPECIES, kind=c_size_t)
    channel_size = int(FUSION_COUPLED_THERMAL_CHANNELS, kind=c_size_t)
    expected_edges = cells_size + 1_c_size_t
    expected_baths = int(FUSION_COUPLED_THERMAL_BASE_BATHS, kind=c_size_t) + &
         inert_size
    expected_logs = expected_baths * species_size

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

    if (present(tables)) then
       if (size(tables, kind=c_size_t) /= channel_size) return
    end if

    ! No C_LOC is formed until all extents above have passed.  In particular,
    ! an empty inert array is never indexed when inert_count is zero.
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
    if (present(tables)) tables_ptr = c_loc(tables(1))
    status = PB11_STATUS_OK
  end subroutine prepare_fast_arguments

  subroutine fusion_coupled_fast_trial(dt_s, options, fast_options, cells, &
       edges_J, thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in) :: electron_energy_J_m3
    real(c_double), intent(in) :: ion_energy_J_m3
    real(c_double), intent(in) :: electron_density_m3
    real(c_double), intent(in), target, contiguous :: &
         thermal_charge_squared(:)
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
    integer(c_int), intent(out) :: status

    type(c_ptr) :: options_ptr, fast_options_ptr, edges_ptr
    type(c_ptr) :: thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr
    type(c_ptr) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr) :: external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr, tables_ptr

    call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
         trial_t_m3, out)
    call prepare_fast_arguments(options, fast_options, cells, edges_J, &
         thermal_number_m3, thermal_charge_squared, inert_count, inert, &
         coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
         trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
         fast_options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status)
    if (status /= PB11_STATUS_OK) return

    out_ptr = c_loc(out)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_fusion_c_coupled_fast_trial(dt_s, options_ptr, &
         fast_options_ptr, cells, edges_ptr, thermal_number_ptr, &
         electron_energy_J_m3, ion_energy_J_m3, electron_density_m3, &
         thermal_charge_squared_ptr, inert_count, inert_ptr, &
         coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         out_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
    end if
  end subroutine fusion_coupled_fast_trial

  subroutine fusion_coupled_fast_table_trial(dt_s, options, fast_options, &
       tables, cells, edges_J, thermal_number_m3, electron_energy_J_m3, &
       ion_energy_J_m3, electron_density_m3, thermal_charge_squared, &
       inert_count, inert, coulomb_logs, old_s_m3, old_t_m3, &
       external_birth_m3_s, escape_s_inv, trial_thermal_number_m3, &
       trial_s_m3, trial_t_m3, out, status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    type(c_ptr), intent(in), target, contiguous :: tables(:)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in) :: electron_energy_J_m3
    real(c_double), intent(in) :: ion_energy_J_m3
    real(c_double), intent(in) :: electron_density_m3
    real(c_double), intent(in), target, contiguous :: &
         thermal_charge_squared(:)
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
    integer(c_int), intent(out) :: status

    type(c_ptr) :: options_ptr, fast_options_ptr, edges_ptr
    type(c_ptr) :: thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr
    type(c_ptr) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr) :: external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr, tables_ptr

    call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
         trial_t_m3, out)
    call prepare_fast_arguments(options, fast_options, cells, edges_J, &
         thermal_number_m3, thermal_charge_squared, inert_count, inert, &
         coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
         trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
         fast_options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status, tables=tables)
    if (status /= PB11_STATUS_OK) return

    out_ptr = c_loc(out)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_fusion_c_coupled_fast_table_trial(dt_s, options_ptr, &
         fast_options_ptr, tables_ptr, cells, edges_ptr, thermal_number_ptr, &
         electron_energy_J_m3, ion_energy_J_m3, electron_density_m3, &
         thermal_charge_squared_ptr, inert_count, inert_ptr, &
         coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         out_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
    end if
  end subroutine fusion_coupled_fast_table_trial

  subroutine fusion_coupled_fast_table_trial_effective_charge( &
       dt_s, options, fast_options, tables, cells, edges_J, &
       thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    type(c_ptr), intent(in), target, contiguous :: tables(:)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in) :: electron_energy_J_m3
    real(c_double), intent(in) :: ion_energy_J_m3
    real(c_double), intent(in) :: electron_density_m3
    real(c_double), intent(in), target, contiguous :: &
         thermal_charge_squared(:)
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
    integer(c_int), intent(out) :: status

    type(c_ptr) :: options_ptr, fast_options_ptr, edges_ptr
    type(c_ptr) :: thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr
    type(c_ptr) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr) :: external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr, tables_ptr

    call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
         trial_t_m3, out)
    call prepare_fast_arguments(options, fast_options, cells, edges_J, &
         thermal_number_m3, thermal_charge_squared, inert_count, inert, &
         coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
         trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
         fast_options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status, tables=tables)
    if (status /= PB11_STATUS_OK) return

    out_ptr = c_loc(out)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_fusion_c_coupled_fast_table_trial_effective_charge( &
         dt_s, options_ptr, fast_options_ptr, tables_ptr, cells, edges_ptr, &
         thermal_number_ptr, electron_energy_J_m3, ion_energy_J_m3, &
         electron_density_m3, thermal_charge_squared_ptr, inert_count, &
         inert_ptr, coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         out_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
    end if
  end subroutine fusion_coupled_fast_table_trial_effective_charge

  subroutine fusion_coupled_fast_trial_diagnosed( &
       dt_s, options, fast_options, cells, edges_J, thermal_number_m3, &
       electron_energy_J_m3, ion_energy_J_m3, electron_density_m3, &
       thermal_charge_squared, inert_count, inert, coulomb_logs, old_s_m3, &
       old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, diagnostics, &
       status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in) :: electron_energy_J_m3
    real(c_double), intent(in) :: ion_energy_J_m3
    real(c_double), intent(in) :: electron_density_m3
    real(c_double), intent(in), target, contiguous :: &
         thermal_charge_squared(:)
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
    integer(c_int), intent(out) :: status

    type(c_ptr) :: options_ptr, fast_options_ptr, edges_ptr
    type(c_ptr) :: thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr
    type(c_ptr) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr) :: external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr, diagnostics_ptr, tables_ptr

    call clear_fast_diagnostics(diagnostics)
    call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
         trial_t_m3, out)
    call prepare_fast_arguments(options, fast_options, cells, edges_J, &
         thermal_number_m3, thermal_charge_squared, inert_count, inert, &
         coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
         trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
         fast_options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status)
    if (status /= PB11_STATUS_OK) return

    out_ptr = c_loc(out)
    diagnostics_ptr = c_loc(diagnostics)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_fusion_c_coupled_fast_trial_diagnosed(dt_s, options_ptr, &
         fast_options_ptr, cells, edges_ptr, thermal_number_ptr, &
         electron_energy_J_m3, ion_energy_J_m3, electron_density_m3, &
         thermal_charge_squared_ptr, inert_count, inert_ptr, &
         coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         out_ptr, diagnostics_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
       call clear_fast_diagnostics(diagnostics)
    end if
  end subroutine fusion_coupled_fast_trial_diagnosed

  subroutine fusion_coupled_fast_table_trial_diagnosed( &
       dt_s, options, fast_options, tables, cells, edges_J, &
       thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, diagnostics, &
       status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    type(c_ptr), intent(in), target, contiguous :: tables(:)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in) :: electron_energy_J_m3
    real(c_double), intent(in) :: ion_energy_J_m3
    real(c_double), intent(in) :: electron_density_m3
    real(c_double), intent(in), target, contiguous :: &
         thermal_charge_squared(:)
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
    integer(c_int), intent(out) :: status

    type(c_ptr) :: options_ptr, fast_options_ptr, edges_ptr
    type(c_ptr) :: thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr
    type(c_ptr) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr) :: external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr, diagnostics_ptr, tables_ptr

    call clear_fast_diagnostics(diagnostics)
    call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
         trial_t_m3, out)
    call prepare_fast_arguments(options, fast_options, cells, edges_J, &
         thermal_number_m3, thermal_charge_squared, inert_count, inert, &
         coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
         trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
         fast_options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status, tables=tables)
    if (status /= PB11_STATUS_OK) return

    out_ptr = c_loc(out)
    diagnostics_ptr = c_loc(diagnostics)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_fusion_c_coupled_fast_table_trial_diagnosed( &
         dt_s, options_ptr, fast_options_ptr, tables_ptr, cells, edges_ptr, &
         thermal_number_ptr, electron_energy_J_m3, ion_energy_J_m3, &
         electron_density_m3, thermal_charge_squared_ptr, inert_count, &
         inert_ptr, coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         out_ptr, diagnostics_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
       call clear_fast_diagnostics(diagnostics)
    end if
  end subroutine fusion_coupled_fast_table_trial_diagnosed

  subroutine fusion_coupled_fast_table_trial_effective_charge_diagnosed( &
       dt_s, options, fast_options, tables, cells, edges_J, &
       thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, diagnostics, &
       status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
    type(fusion_fast_target_options_v1), intent(in), target :: fast_options
    type(c_ptr), intent(in), target, contiguous :: tables(:)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: thermal_number_m3(:)
    real(c_double), intent(in) :: electron_energy_J_m3
    real(c_double), intent(in) :: ion_energy_J_m3
    real(c_double), intent(in) :: electron_density_m3
    real(c_double), intent(in), target, contiguous :: &
         thermal_charge_squared(:)
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
    integer(c_int), intent(out) :: status

    type(c_ptr) :: options_ptr, fast_options_ptr, edges_ptr
    type(c_ptr) :: thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr
    type(c_ptr) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr) :: external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr, diagnostics_ptr, tables_ptr

    call clear_fast_diagnostics(diagnostics)
    call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
         trial_t_m3, out)
    call prepare_fast_arguments(options, fast_options, cells, edges_J, &
         thermal_number_m3, thermal_charge_squared, inert_count, inert, &
         coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
         trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
         fast_options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status, tables=tables)
    if (status /= PB11_STATUS_OK) return

    out_ptr = c_loc(out)
    diagnostics_ptr = c_loc(diagnostics)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))
    status = c_fusion_c_coupled_fast_table_trial_effective_charge_diagnosed( &
         dt_s, options_ptr, fast_options_ptr, tables_ptr, cells, edges_ptr, &
         thermal_number_ptr, electron_energy_J_m3, ion_energy_J_m3, &
         electron_density_m3, thermal_charge_squared_ptr, inert_count, &
         inert_ptr, coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
         escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
         out_ptr, diagnostics_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_fast_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
       call clear_fast_diagnostics(diagnostics)
    end if
  end subroutine fusion_coupled_fast_table_trial_effective_charge_diagnosed

end module fusion_coupled_fast_fortran
