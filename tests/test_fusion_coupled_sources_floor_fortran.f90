program test_fusion_coupled_sources_floor_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_null_ptr, &
       c_loc, c_ptr, c_size_t, c_sizeof
  use fusion_coupled_sources_floor_fortran, only : &
       floor_call => c_fusion_c_coupled_sources_floor_trial, &
       floor_limits_type => fusion_coupled_floor_limits_v1, &
       floor_ledger_type => fusion_birth_floor_ledger_v1, &
       options_type => fusion_coupled_thermal_options_v1, &
       fast_options_type => fusion_fast_target_options_v1, &
       coupled_output_type => fusion_coupled_thermal_v1, &
       diagnostics_type => fusion_handoff_diagnostics_v1, &
       usage_type => fusion_beam_table_usage_v1, &
       status_ok => PB11_STATUS_OK, status_bad => PB11_STATUS_OUT_OF_RANGE
  implicit none

  integer(c_int), parameter :: cells = 8_c_int
  integer(c_int), parameter :: species = 6_c_int
  integer(c_int), parameter :: baths = 7_c_int
  real(c_double), parameter :: ev_j = 1.602176634e-19_c_double
  real(c_double), parameter :: mev_j = 1.0e6_c_double * ev_j
  real(c_double), parameter :: dt_s = 1.0e-5_c_double

  real(c_double), target :: edges(cells + 1_c_int)
  real(c_double), target :: thermal_number(species)
  real(c_double), target :: thermal_charge_squared(species)
  real(c_double), target :: coulomb_logs(baths,species)
  real(c_double), target :: old_s(cells,species), old_t(cells,species)
  real(c_double), target :: external_birth(cells,species), escape(cells,species)
  real(c_double), target :: new_thermal_number(species)
  real(c_double), target :: new_s(cells,species), new_t(cells,species)
  real(c_double), target :: new_thermal_bad(species)
  real(c_double), target :: new_s_bad(cells,species), new_t_bad(cells,species)
  real(c_double), target :: electron_energy, ion_energy, electron_density

  type(options_type), target :: options
  type(fast_options_type), target :: fast_options
  type(coupled_output_type), target :: output, output_bad
  type(diagnostics_type), target :: diagnostics, diagnostics_bad
  type(usage_type), target :: usage, usage_bad
  type(floor_limits_type), target :: limits
  type(floor_ledger_type), target :: ledger, ledger_bad
  integer(c_int) :: status, status_bad_call
  integer :: failures

  failures = 0
  call check(c_sizeof(options) == 160_c_size_t, &
       'coupled options retain the existing 160-byte C layout')
  call check(c_sizeof(output) == 1224_c_size_t, &
       'coupled output retains the existing 1224-byte C layout')
  call check(c_sizeof(diagnostics) == 264_c_size_t, &
       'handoff diagnostics retain the existing 264-byte C layout')
  call check(c_sizeof(limits) == 16_c_size_t, &
       'floor limits has the two-double C layout')
  call check(c_sizeof(ledger) == 296_c_size_t, &
       'floor ledger has the expected six arrays plus residual')

  call make_inputs()
  limits%max_center_over_ion_kT = 1.0_c_double
  limits%max_ion_energy_fraction = 1.0_c_double
  call clear_sentinels()

  ! All optional borrowed arrays are absent.  Every channel is disabled, so
  ! this is the direct no-source path and the floor ledger must remain empty.
  status = floor_call(dt_s, c_loc(options), c_loc(fast_options), c_null_ptr, &
       0_c_int, c_null_ptr, 0_c_int, cells, c_loc(edges), &
       c_loc(thermal_number), electron_energy, ion_energy, electron_density, &
       c_loc(thermal_charge_squared), 0_c_int, c_null_ptr, c_loc(coulomb_logs), &
       c_loc(old_s), c_loc(old_t), c_loc(external_birth), c_loc(escape), &
       c_loc(new_thermal_number), c_loc(new_s), c_loc(new_t), c_loc(output), &
       c_loc(diagnostics), c_loc(usage), c_loc(limits), c_loc(ledger))
  call check(status == status_ok, 'null optional pointers and inactive channels succeed')
  call check(all(ledger%born_number_m3 == 0.0_c_double) .and. &
       all(ledger%born_energy_J_m3 == 0.0_c_double) .and. &
       all(ledger%mapped_number_m3 == 0.0_c_double) .and. &
       all(ledger%ion_energy_correction_J_m3 == 0.0_c_double), &
       'inactive channels leave the floor ledger empty')
  call check(all(output%ledger%events_m3 == 0.0_c_double) .and. &
       all(new_s == 0.0_c_double) .and. all(new_t == 0.0_c_double), &
       'inactive channels produce no source or fast population')

  ! A rejected limits object must clear every public output.  This also checks
  ! that the final two arguments are passed as pointers, with no hidden value
  ! copy or Fortran descriptor.
  limits%max_center_over_ion_kT = -1.0_c_double
  new_thermal_bad = -7.0_c_double
  new_s_bad = -7.0_c_double
  new_t_bad = -7.0_c_double
  output_bad%electron_energy_J_m3 = -7.0_c_double
  diagnostics_bad%tested = 1_c_int
  usage_bad%table_evaluations = -7_c_int
  ledger_bad%born_number_m3 = -7.0_c_double
  status_bad_call = floor_call(dt_s, c_loc(options), c_loc(fast_options), &
       c_null_ptr, 0_c_int, c_null_ptr, 0_c_int, cells, c_loc(edges), &
       c_loc(thermal_number), electron_energy, ion_energy, electron_density, &
       c_loc(thermal_charge_squared), 0_c_int, c_null_ptr, c_loc(coulomb_logs), &
       c_loc(old_s), c_loc(old_t), c_loc(external_birth), c_loc(escape), &
       c_loc(new_thermal_bad), c_loc(new_s_bad), c_loc(new_t_bad), &
       c_loc(output_bad), c_loc(diagnostics_bad), c_loc(usage_bad), &
       c_loc(limits), c_loc(ledger_bad))
  call check(status_bad_call == status_bad, 'invalid floor limit is rejected')
  call check(all(new_thermal_bad == 0.0_c_double) .and. &
       all(new_s_bad == 0.0_c_double) .and. all(new_t_bad == 0.0_c_double) .and. &
       all(ledger_bad%born_number_m3 == 0.0_c_double) .and. &
       output_bad%electron_energy_J_m3 == 0.0_c_double .and. &
       all(diagnostics_bad%tested == 0_c_int) .and. &
       usage_bad%table_evaluations == 0_c_int, &
       'rejected call clears arrays, structs, diagnostics, usage, and ledger')

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' coupled-floor Fortran check(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'Fortran coupled-floor direct ABI smoke test passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  subroutine make_inputs()
    integer :: j

    do j = 1, cells + 1
       edges(j) = 25.0_c_double * mev_j * real(j - 1, c_double) / &
            real(cells, c_double)
    end do
    thermal_number = 0.0_c_double
    thermal_number(3) = 1.0e19_c_double
    thermal_charge_squared = [1.0_c_double, 1.0_c_double, 1.0_c_double, &
         4.0_c_double, 4.0_c_double, 25.0_c_double]
    coulomb_logs = 15.0_c_double
    old_s = 0.0_c_double
    old_t = 0.0_c_double
    external_birth = 0.0_c_double
    escape = 0.0_c_double
    electron_density = 1.0e19_c_double
    electron_energy = 1.5_c_double * electron_density * 5.0e3_c_double * ev_j
    ion_energy = 1.5_c_double * thermal_number(3) * 10.0e3_c_double * ev_j

    options%birth%relative_max_J = 5.0e6_c_double * ev_j
    options%birth%cm_max_kT = 40.0_c_double
    options%birth%ground_state_q_J = 91.84e3_c_double * ev_j
    options%birth%cutoff_J = 2.0e3_c_double * ev_j
    options%birth%l1_fraction = 0.76_c_double
    options%birth%relative_phase = 0.0_c_double
    options%birth%narrow_peak_fraction = 0.051_c_double
    options%birth%continuum_peak_scale = 1.0_c_double
    options%birth%continuation = 1_c_int
    options%birth%pb_low = 0_c_int
    options%birth%remainder_policy = 0_c_int
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
    fast_options%angular_order = 16_c_int
    fast_options%angular_max_exponent = 40.0_c_double
  end subroutine make_inputs

  subroutine clear_sentinels()
    new_thermal_number = -7.0_c_double
    new_s = -7.0_c_double
    new_t = -7.0_c_double
    output%electron_energy_J_m3 = -7.0_c_double
    diagnostics%tested = 1_c_int
    usage%table_evaluations = -7_c_int
    ledger%born_number_m3 = -7.0_c_double
  end subroutine clear_sentinels

end program test_fusion_coupled_sources_floor_fortran
