program test_fusion_coupled_fast_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_ptr
  use fusion_coupled_fast_fortran, only : &
       fast_options_type => fusion_fast_target_options_v1, &
       coupled_options_type => fusion_coupled_thermal_options_v1, &
       coupled_output_type => fusion_coupled_thermal_v1, &
       inert_type => fusion_inert_ion_v1, &
       fast_trial => fusion_coupled_fast_trial, &
       fast_table_trial => fusion_coupled_fast_table_trial, &
       fast_effective_trial => fusion_coupled_fast_table_trial_effective_charge, &
       fast_ok => PB11_STATUS_OK, fast_invalid => PB11_STATUS_INVALID_ARGUMENT
  use fusion_coupled_thermal_fortran, only : &
       thermal_trial => fusion_coupled_thermal_trial
  implicit none

  integer(c_int), parameter :: cells = 1_c_int
  integer(c_int), parameter :: species = 6_c_int
  integer(c_int), parameter :: channels = 5_c_int
  integer(c_int), parameter :: baths = 7_c_int
  real(c_double), parameter :: ev_j = 1.602176634e-19_c_double
  real(c_double), parameter :: dt_s = 1.0e-9_c_double
  real(c_double) :: edges(2), edges_bad(1)
  real(c_double) :: thermal_number(species), thermal_charge_squared(species)
  real(c_double) :: coulomb_logs(baths,species)
  real(c_double) :: old_s(cells,species), old_t(cells,species)
  real(c_double) :: external_birth(cells,species), escape(cells,species)
  real(c_double) :: trial_number_a(species), trial_s_a(cells,species)
  real(c_double) :: trial_t_a(cells,species), trial_number_b(species)
  real(c_double) :: trial_s_b(cells,species), trial_t_b(cells,species)
  real(c_double) :: trial_number_bad(species), trial_s_bad(cells,species)
  real(c_double) :: trial_t_bad(cells,species)
  real(c_double) :: electron_energy, ion_energy, electron_density
  type(coupled_options_type) :: options
  type(fast_options_type) :: fast_options
  type(inert_type) :: inert(0)
  type(coupled_output_type) :: output_a, output_b, output_bad
  type(c_ptr) :: tables(channels)
  integer(c_int) :: status_a, status_b, status_bad
  integer :: i

  edges = [0.0_c_double, 1.0e-12_c_double]
  edges_bad = [0.0_c_double]
  thermal_number = 0.0_c_double
  thermal_number(2) = 5.0e19_c_double
  thermal_charge_squared = [1.0_c_double, 1.0_c_double, 1.0_c_double, &
       4.0_c_double, 4.0_c_double, 25.0_c_double]
  coulomb_logs = 15.0_c_double
  old_s = 0.0_c_double
  old_t = 0.0_c_double
  external_birth = 0.0_c_double
  escape = 0.0_c_double
  electron_density = 1.0e20_c_double
  electron_energy = 1.5_c_double * electron_density * 10.0e3_c_double * ev_j
  ion_energy = 1.5_c_double * thermal_number(2) * 10.0e3_c_double * ev_j

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
  fast_options%angular_order = 8_c_int
  fast_options%angular_max_exponent = 40.0_c_double
  tables = c_null_ptr

  ! A wrong edge extent must be rejected and all caller-visible outputs clear
  ! before any C_LOC or C++ entry point is reached.
  trial_number_bad = -7.0_c_double
  trial_s_bad = -7.0_c_double
  trial_t_bad = -7.0_c_double
  output_bad%electron_energy_J_m3 = -7.0_c_double
  call fast_trial(dt_s, options, fast_options, cells, edges_bad, &
       thermal_number, electron_energy, ion_energy, electron_density, &
       thermal_charge_squared, 0_c_int, inert, coulomb_logs, old_s, old_t, &
       external_birth, escape, trial_number_bad, trial_s_bad, trial_t_bad, &
       output_bad, status_bad)
  if (status_bad /= fast_invalid .or. any(trial_number_bad /= 0.0_c_double) &
       .or. any(trial_s_bad /= 0.0_c_double) .or. &
       any(trial_t_bad /= 0.0_c_double) .or. &
       output_bad%electron_energy_J_m3 /= 0.0_c_double) error stop 1

  call thermal_trial(dt_s, options, cells, edges, thermal_number, &
       electron_energy, ion_energy, electron_density, thermal_charge_squared, &
       0_c_int, inert, coulomb_logs, old_s, old_t, external_birth, escape, &
       trial_number_a, trial_s_a, trial_t_a, output_a, status_a)
  call fast_trial(dt_s, options, fast_options, cells, edges, thermal_number, &
       electron_energy, ion_energy, electron_density, thermal_charge_squared, &
       0_c_int, inert, coulomb_logs, old_s, old_t, external_birth, escape, &
       trial_number_b, trial_s_b, trial_t_b, output_b, status_b)
  if (status_a /= fast_ok .or. status_b /= fast_ok) error stop 2
  if (any(trial_number_a /= trial_number_b) .or. &
       any(trial_s_a /= trial_s_b) .or. any(trial_t_a /= trial_t_b) .or. &
       output_a%electron_energy_J_m3 /= output_b%electron_energy_J_m3 .or. &
       output_a%ion_energy_J_m3 /= output_b%ion_energy_J_m3) error stop 3

  call fast_table_trial(dt_s, options, fast_options, tables, cells, edges, &
       thermal_number, electron_energy, ion_energy, electron_density, &
       thermal_charge_squared, 0_c_int, inert, coulomb_logs, old_s, old_t, &
       external_birth, escape, trial_number_b, trial_s_b, trial_t_b, output_b, &
       status_b)
  if (status_b /= fast_ok) error stop 4
  call fast_effective_trial(dt_s, options, fast_options, tables, cells, edges, &
       thermal_number, electron_energy, ion_energy, electron_density, &
       thermal_charge_squared, 0_c_int, inert, coulomb_logs, old_s, old_t, &
       external_birth, escape, trial_number_b, trial_s_b, trial_t_b, output_b, &
       status_b)
  if (status_b /= fast_ok) error stop 5

  do i = 1, species
     if (trial_number_b(i) /= trial_number_a(i)) error stop 6
  end do
  call run_nonzero_fast_dt()
  print *, 'fusion_coupled_fast_fortran smoke: PASS'

contains

  subroutine run_nonzero_fast_dt()
    ! A bounded integration fixture with one enabled fast-target channel.
    ! Keep the first edge at zero and the first positive edge at 2 micro-eV;
    ! the remaining 256 intervals are geometric through 25 MeV.
    integer, parameter :: n = 257
    real(c_double), parameter :: target_energy = 100.0e3_c_double * ev_j
    real(c_double), parameter :: max_energy = 25.0e6_c_double * ev_j
    real(c_double), parameter :: first_edge = 2.0e-6_c_double * ev_j
    real(c_double), target :: local_edges(n + 1)
    real(c_double), target :: local_thermal_number(species)
    real(c_double), target :: local_charge_squared(species)
    real(c_double), target :: local_logs(baths,species)
    real(c_double), target :: local_old_s(n,species), local_old_t(n,species)
    real(c_double), target :: local_external(n,species), local_escape(n,species)
    real(c_double), target :: local_trial_number(species)
    real(c_double), target :: local_trial_s(n,species), local_trial_t(n,species)
    type(inert_type), target :: local_inert(0)
    type(coupled_options_type), target :: local_options
    type(fast_options_type), target :: local_fast_options
    type(coupled_output_type), target :: local_output
    real(c_double) :: ratio, center, distance, best_distance
    real(c_double) :: electron_density_local, electron_energy_local
    real(c_double) :: ion_energy_local, events, fast_d, thermal_t
    integer :: i, nearest
    integer(c_int) :: local_status

    local_edges(1) = 0.0_c_double
    local_edges(2) = first_edge
    ratio = (max_energy / first_edge) ** &
         (1.0_c_double / real(n - 1, c_double))
    do i = 3, n + 1
       local_edges(i) = first_edge * ratio ** real(i - 2, c_double)
    end do
    local_edges(n + 1) = max_energy

    nearest = 1
    best_distance = huge(1.0_c_double)
    do i = 1, n
       center = 0.5_c_double * (local_edges(i) + local_edges(i + 1))
       distance = abs(center - target_energy)
       if (distance < best_distance) then
          best_distance = distance
          nearest = i
       end if
    end do

    local_thermal_number = 0.0_c_double
    local_thermal_number(3) = 1.0e19_c_double
    local_charge_squared = [1.0_c_double, 1.0_c_double, 1.0_c_double, &
         4.0_c_double, 4.0_c_double, 25.0_c_double]
    local_logs = 15.0_c_double
    local_old_s = 0.0_c_double
    local_old_t = 0.0_c_double
    local_old_s(nearest, 2) = 1.0e15_c_double
    local_old_t(nearest, 2) = 2.0e15_c_double
    local_external = 0.0_c_double
    local_escape = 0.0_c_double

    local_options%birth%relative_max_J = 5.0e6_c_double * ev_j
    local_options%birth%cm_max_kT = 40.0_c_double
    local_options%birth%ground_state_q_J = 91.84e3_c_double * ev_j
    local_options%birth%cutoff_J = 2.0e3_c_double * ev_j
    local_options%birth%l1_fraction = 0.76_c_double
    local_options%birth%relative_phase = 0.0_c_double
    local_options%birth%narrow_peak_fraction = 0.051_c_double
    local_options%birth%continuum_peak_scale = 1.0_c_double
    local_options%birth%continuation = 1_c_int
    local_options%birth%pb_low = 0_c_int
    local_options%birth%remainder_policy = 0_c_int
    local_options%birth%broad_mode = 13_c_int
    local_options%birth%fsci_policy = 0_c_int
    local_options%birth%relative_order = 8_c_int
    local_options%birth%cm_order = 8_c_int
    local_options%birth%nq = 8_c_int
    local_options%birth%ncos = 8_c_int
    local_options%max_source_rate_error = 2.0e-3_c_double
    local_options%max_source_debit_error = 2.0e-3_c_double
    local_options%handoff_max_L1 = 1.0e-3_c_double
    local_options%handoff_max_mean_error = 1.0e-3_c_double
    local_options%channels = 0_c_int
    local_options%handoff_enabled = 0_c_int
    local_fast_options%channels = 0_c_int
    local_fast_options%channels(4) = 1_c_int
    local_fast_options%angular_order = 16_c_int
    local_fast_options%angular_max_exponent = 40.0_c_double

    electron_density_local = 1.0e19_c_double
    electron_energy_local = 1.5_c_double * electron_density_local * &
         5.0e3_c_double * ev_j
    ion_energy_local = 1.5_c_double * local_thermal_number(3) * &
         10.0e3_c_double * ev_j

    call fast_trial(1.0e-4_c_double, local_options, local_fast_options, &
         int(n, c_int), local_edges, local_thermal_number, &
         electron_energy_local, ion_energy_local, electron_density_local, &
         local_charge_squared, 0_c_int, local_inert, local_logs, local_old_s, &
         local_old_t, local_external, local_escape, local_trial_number, &
         local_trial_s, local_trial_t, local_output, local_status)
    if (local_status /= fast_ok) then
       write(*, '(A,I0)') 'nonzero fast DT wrapper status=', local_status
       error stop 7
    end if

    events = local_output%ledger%events_m3(4)
    fast_d = local_output%ledger%fast_consumed_number_m3(2)
    thermal_t = local_output%ledger%thermal_consumed_number_m3(3)
    if (events <= 0.0_c_double .or. fast_d <= 0.0_c_double .or. &
         thermal_t <= 0.0_c_double .or. &
         local_output%ledger%neutron_number_m3 <= 0.0_c_double .or. &
         local_output%ledger%nuclear_born_number_m3(5) <= 0.0_c_double) then
       write(*, '(A,5(1X,ES14.6))') 'nonzero fast DT amounts=', events, &
            fast_d, thermal_t, local_output%ledger%neutron_number_m3, &
            local_output%ledger%nuclear_born_number_m3(5)
       error stop 8
    end if
    if (abs(fast_d - events) > 1.0e-8_c_double * &
         max(abs(fast_d), abs(events)) .or. &
         abs(thermal_t - events) > 1.0e-8_c_double * &
         max(abs(thermal_t), abs(events)) .or. &
         abs(local_output%ledger%neutron_number_m3 - events) > &
         1.0e-8_c_double * max(abs(local_output%ledger%neutron_number_m3), &
         abs(events)) .or. &
       abs(local_output%ledger%nuclear_born_number_m3(5) - events) > &
       1.0e-8_c_double * &
         max(abs(local_output%ledger%nuclear_born_number_m3(5)), abs(events))) then
       write(*, '(A,4(1X,ES14.6))') 'nonzero fast DT mismatch=', events, &
            fast_d, thermal_t, local_output%ledger%neutron_number_m3
       error stop 9
    end if
  end subroutine run_nonzero_fast_dt
end program test_fusion_coupled_fast_fortran
