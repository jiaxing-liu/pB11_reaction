module fusion_coupled_thermal_fortran
  !! ISO_C_BINDING wrapper for the stateless coupled thermal/kinetic trial.
  !!
  !! The nested C layouts are imported from the existing Fortran bindings so
  !! that there is one source of truth for the thermal-birth options and the
  !! source ledger.  This layer only validates Fortran extents, translates
  !! the species-major C buffers to Fortran shapes, and clears every output
  !! on a rejected call.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_ptr, &
       c_null_ptr, c_loc
  use fusion_thermal_birth_fortran, only : fusion_thermal_birth_options_v1
  use fusion_source_state_fortran, only : fusion_source_ledger_v1, &
       PB11_STATUS_OK, PB11_STATUS_NULL_OUTPUT, PB11_STATUS_INVALID_ARGUMENT, &
       PB11_STATUS_OUT_OF_RANGE, PB11_STATUS_NUMERICAL_FAILURE, &
       PB11_STATUS_EXCEPTION, PB11_STATUS_UNKNOWN_METHOD
  implicit none
  private

  integer(c_int), parameter, public :: FUSION_COUPLED_THERMAL_SPECIES = &
       6_c_int
  integer(c_int), parameter, public :: FUSION_COUPLED_THERMAL_CHANNELS = &
       5_c_int
  integer(c_int), parameter, public :: FUSION_COUPLED_THERMAL_MAX_CELLS = &
       100000_c_int
  integer(c_int), parameter, public :: FUSION_COUPLED_THERMAL_MAX_INERT = &
       32_c_int
  integer(c_int), parameter, public :: FUSION_COUPLED_THERMAL_BASE_BATHS = &
       7_c_int

  ! Keep the status names available from this binding, as in the other
  ! Fortran wrappers.  Their definitions come from source_state_fortran.
  public :: PB11_STATUS_OK
  public :: PB11_STATUS_NULL_OUTPUT
  public :: PB11_STATUS_INVALID_ARGUMENT
  public :: PB11_STATUS_OUT_OF_RANGE
  public :: PB11_STATUS_NUMERICAL_FAILURE
  public :: PB11_STATUS_EXCEPTION
  public :: PB11_STATUS_UNKNOWN_METHOD

  ! Exact C layout: three consecutive c_double fields.
  type, bind(C), public :: fusion_inert_ion_v1
     real(c_double) :: density_m3
     real(c_double) :: mass_kg
     real(c_double) :: mean_charge_squared
  end type fusion_inert_ion_v1

  ! Exact C layout: an imported 104-byte birth-options object, four doubles,
  ! and six c_int values.  BIND(C) supplies the final ABI padding.
  type, bind(C), public :: fusion_coupled_thermal_options_v1
     type(fusion_thermal_birth_options_v1) :: birth
     real(c_double) :: max_source_rate_error
     real(c_double) :: max_source_debit_error
     real(c_double) :: handoff_max_L1
     real(c_double) :: handoff_max_mean_error
     integer(c_int) :: channels(5)
     integer(c_int) :: handoff_enabled
  end type fusion_coupled_thermal_options_v1

  ! Exact C layout: imported source ledger followed by all coupled thermal,
  ! residual, and handoff diagnostics.
  type, bind(C), public :: fusion_coupled_thermal_v1
     type(fusion_source_ledger_v1) :: ledger
     real(c_double) :: inert_ion_heat_J_m3(6)
     real(c_double) :: electron_energy_J_m3
     real(c_double) :: ion_energy_J_m3
     real(c_double) :: particle_residual_m3(6)
     real(c_double) :: energy_residual_J_m3
     real(c_double) :: max_source_rate_discrepancy
     real(c_double) :: max_source_debit_discrepancy
     real(c_double) :: handoff_L1(6)
     real(c_double) :: handoff_mean_error(6)
     integer(c_int) :: handoff_projected(6)
  end type fusion_coupled_thermal_v1

  ! Exact C layout: eight consecutive c_double values.
  type, bind(C), public :: fusion_thermal_increment_v1
     real(c_double) :: thermal_number_m3(6)
     real(c_double) :: electron_energy_J_m3
     real(c_double) :: ion_energy_J_m3
  end type fusion_thermal_increment_v1

  public :: fusion_coupled_thermal_trial
  public :: fusion_coupled_thermal_table_trial
  public :: fusion_coupled_thermal_increment

  interface
     function c_fusion_c_coupled_thermal_trial(dt_s, options, cells, edges, &
          thermal_number, electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, trial_thermal_number, trial_s, &
          trial_t, out) bind(C, name="fusion_c_coupled_thermal_trial") &
          result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
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
     end function c_fusion_c_coupled_thermal_trial

     function c_fusion_c_coupled_thermal_table_trial(dt_s, options, tables, &
          cells, edges, thermal_number, electron_energy, ion_energy, &
          electron_density, thermal_charge_squared, inert_count, inert, &
          coulomb_logs, old_s, old_t, external_birth, escape, &
          trial_thermal_number, trial_s, trial_t, out) bind(C, &
          name="fusion_c_coupled_thermal_table_trial") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
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
     end function c_fusion_c_coupled_thermal_table_trial

     function c_fusion_c_coupled_thermal_increment(ledger, inert_heat, out) &
          bind(C, name="fusion_c_coupled_thermal_increment") result(status)
       import :: c_int, c_ptr
       type(c_ptr), value :: ledger
       type(c_ptr), value :: inert_heat
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_c_coupled_thermal_increment
  end interface

contains

  subroutine clear_ledger(ledger)
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
  end subroutine clear_ledger

  subroutine clear_coupled(out)
    type(fusion_coupled_thermal_v1), intent(out) :: out

    call clear_ledger(out%ledger)
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
  end subroutine clear_coupled

  subroutine clear_increment(out)
    type(fusion_thermal_increment_v1), intent(out) :: out

    out%thermal_number_m3 = 0.0_c_double
    out%electron_energy_J_m3 = 0.0_c_double
    out%ion_energy_J_m3 = 0.0_c_double
  end subroutine clear_increment

  logical function valid_cells(cells)
    integer(c_int), intent(in) :: cells

    valid_cells = cells >= 1_c_int .and. &
         cells <= FUSION_COUPLED_THERMAL_MAX_CELLS
  end function valid_cells

  subroutine clear_trial_outputs(trial_thermal_number_m3, trial_s_m3, &
       trial_t_m3, out)
    real(c_double), intent(out) :: trial_thermal_number_m3(:)
    real(c_double), intent(out) :: trial_s_m3(:,:)
    real(c_double), intent(out) :: trial_t_m3(:,:)
    type(fusion_coupled_thermal_v1), intent(out) :: out

    call clear_coupled(out)
    trial_thermal_number_m3 = 0.0_c_double
    trial_s_m3 = 0.0_c_double
    trial_t_m3 = 0.0_c_double
  end subroutine clear_trial_outputs

  subroutine prepare_trial_arguments(options, cells, edges_J, &
       thermal_number_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, options_ptr, &
       edges_ptr, thermal_number_ptr, thermal_charge_squared_ptr, inert_ptr, &
       coulomb_logs_ptr, old_s_ptr, old_t_ptr, external_birth_ptr, &
       escape_ptr, trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, &
       tables_ptr, status, tables)
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
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
    type(c_ptr), intent(out) :: options_ptr, edges_ptr, thermal_number_ptr
    type(c_ptr), intent(out) :: thermal_charge_squared_ptr, inert_ptr
    type(c_ptr), intent(out) :: coulomb_logs_ptr, old_s_ptr, old_t_ptr
    type(c_ptr), intent(out) :: external_birth_ptr, escape_ptr
    type(c_ptr), intent(out) :: trial_thermal_number_ptr, trial_s_ptr
    type(c_ptr), intent(out) :: trial_t_ptr, tables_ptr
    integer(c_int), intent(out) :: status
    type(c_ptr), intent(in), target, contiguous, optional :: tables(:)

    integer(c_int) :: expected_baths
    integer(c_size_t) :: expected_edges, expected_logs

    options_ptr = c_null_ptr
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

    ! Validate every extent, including the optional table handle array, before
    ! taking any C_LOC.  A zero-size inert/table actual is safe on rejection.
    if (.not. valid_cells(cells)) return
    expected_edges = int(cells, c_size_t) + 1_c_size_t
    if (size(edges_J, kind=c_size_t) /= expected_edges) return
    if (size(thermal_number_m3, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (size(thermal_charge_squared, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (size(trial_thermal_number_m3, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return

    if (inert_count < 0_c_int .or. &
         inert_count > FUSION_COUPLED_THERMAL_MAX_INERT) return
    if (size(inert, kind=c_size_t) /= int(inert_count, c_size_t)) return

    expected_baths = FUSION_COUPLED_THERMAL_BASE_BATHS + inert_count
    if (size(coulomb_logs, 1, kind=c_size_t) /= &
         int(expected_baths, c_size_t)) return
    if (size(coulomb_logs, 2, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    expected_logs = int(expected_baths, c_size_t) * &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)
    if (size(coulomb_logs, kind=c_size_t) /= expected_logs) return

    if (size(old_s_m3, 1, kind=c_size_t) /= int(cells, c_size_t)) return
    if (size(old_s_m3, 2, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (size(old_t_m3, 1, kind=c_size_t) /= int(cells, c_size_t)) return
    if (size(old_t_m3, 2, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (size(external_birth_m3_s, 1, kind=c_size_t) /= &
         int(cells, c_size_t)) return
    if (size(external_birth_m3_s, 2, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (size(escape_s_inv, 1, kind=c_size_t) /= int(cells, c_size_t)) return
    if (size(escape_s_inv, 2, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (size(trial_s_m3, 1, kind=c_size_t) /= int(cells, c_size_t)) return
    if (size(trial_s_m3, 2, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (size(trial_t_m3, 1, kind=c_size_t) /= int(cells, c_size_t)) return
    if (size(trial_t_m3, 2, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return
    if (present(tables)) then
       if (size(tables, kind=c_size_t) /= &
            int(FUSION_COUPLED_THERMAL_CHANNELS, c_size_t)) return
    end if

    options_ptr = c_loc(options)
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
    status = PB11_STATUS_OK
    if (present(tables)) tables_ptr = c_loc(tables(1))
  end subroutine prepare_trial_arguments

  subroutine fusion_coupled_thermal_trial(dt_s, options, cells, edges_J, &
       thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
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

    type(c_ptr) :: options_ptr, edges_ptr, thermal_number_ptr
    type(c_ptr) :: thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr
    type(c_ptr) :: old_s_ptr, old_t_ptr, external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr, tables_ptr

    call clear_trial_outputs(trial_thermal_number_m3, trial_s_m3, trial_t_m3, out)
    call prepare_trial_arguments(options, cells, edges_J, thermal_number_m3, &
         thermal_charge_squared, inert_count, inert, coulomb_logs, old_s_m3, &
         old_t_m3, external_birth_m3_s, escape_s_inv, trial_thermal_number_m3, &
         trial_s_m3, trial_t_m3, options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status)
    if (status /= PB11_STATUS_OK) return
    out_ptr = c_loc(out)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))

    status = c_fusion_c_coupled_thermal_trial(dt_s, options_ptr, cells, &
         edges_ptr, thermal_number_ptr, electron_energy_J_m3, &
         ion_energy_J_m3, electron_density_m3, thermal_charge_squared_ptr, &
         inert_count, inert_ptr, coulomb_logs_ptr, old_s_ptr, old_t_ptr, &
         external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_trial_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
    end if
  end subroutine fusion_coupled_thermal_trial

  subroutine fusion_coupled_thermal_table_trial(dt_s, options, tables, cells, &
       edges_J, thermal_number_m3, electron_energy_J_m3, ion_energy_J_m3, &
       electron_density_m3, thermal_charge_squared, inert_count, inert, &
       coulomb_logs, old_s_m3, old_t_m3, external_birth_m3_s, escape_s_inv, &
       trial_thermal_number_m3, trial_s_m3, trial_t_m3, out, status)
    real(c_double), intent(in) :: dt_s
    type(fusion_coupled_thermal_options_v1), intent(in), target :: options
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

    type(c_ptr) :: options_ptr, tables_ptr, edges_ptr, thermal_number_ptr
    type(c_ptr) :: thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr
    type(c_ptr) :: old_s_ptr, old_t_ptr, external_birth_ptr, escape_ptr
    type(c_ptr) :: trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr
    type(c_ptr) :: out_ptr

    call clear_trial_outputs(trial_thermal_number_m3, trial_s_m3, trial_t_m3, out)
    call prepare_trial_arguments(options, cells, edges_J, thermal_number_m3, &
         thermal_charge_squared, inert_count, inert, coulomb_logs, old_s_m3, &
         old_t_m3, external_birth_m3_s, escape_s_inv, trial_thermal_number_m3, &
         trial_s_m3, trial_t_m3, options_ptr, edges_ptr, thermal_number_ptr, &
         thermal_charge_squared_ptr, inert_ptr, coulomb_logs_ptr, old_s_ptr, &
         old_t_ptr, external_birth_ptr, escape_ptr, trial_thermal_number_ptr, &
         trial_s_ptr, trial_t_ptr, tables_ptr, status, &
         tables=tables)
    if (status /= PB11_STATUS_OK) return
    out_ptr = c_loc(out)
    inert_ptr = c_null_ptr
    if (inert_count > 0_c_int) inert_ptr = c_loc(inert(1))

    status = c_fusion_c_coupled_thermal_table_trial(dt_s, options_ptr, &
         tables_ptr, cells, edges_ptr, thermal_number_ptr, &
         electron_energy_J_m3, ion_energy_J_m3, electron_density_m3, &
         thermal_charge_squared_ptr, inert_count, inert_ptr, coulomb_logs_ptr, &
         old_s_ptr, old_t_ptr, external_birth_ptr, escape_ptr, &
         trial_thermal_number_ptr, trial_s_ptr, trial_t_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       call clear_trial_outputs(trial_thermal_number_m3, trial_s_m3, &
            trial_t_m3, out)
    end if
  end subroutine fusion_coupled_thermal_table_trial

  subroutine fusion_coupled_thermal_increment(ledger, inert_heat_J_m3, out, &
       status)
    type(fusion_source_ledger_v1), intent(in), target :: ledger
    real(c_double), intent(in), target, contiguous :: inert_heat_J_m3(:)
    type(fusion_thermal_increment_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_increment(out)
    status = PB11_STATUS_INVALID_ARGUMENT
    if (size(inert_heat_J_m3, kind=c_size_t) /= &
         int(FUSION_COUPLED_THERMAL_SPECIES, c_size_t)) return

    status = c_fusion_c_coupled_thermal_increment(c_loc(ledger), &
         c_loc(inert_heat_J_m3(1)), c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_increment(out)
  end subroutine fusion_coupled_thermal_increment

end module fusion_coupled_thermal_fortran
