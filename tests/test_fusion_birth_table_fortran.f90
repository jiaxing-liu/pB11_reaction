program test_fusion_birth_table_fortran
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_null_ptr, c_ptr, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_thermal_birth_fortran, only : &
       fusion_thermal_birth_options_v1, FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, &
       FUSION_PB_LOW_TB, FUSION_PB_REMAINDER_ENTRANCE_PROXY
  use fusion_birth_table_fortran, only : PB11_STATUS_OK, &
       PB11_STATUS_INVALID_ARGUMENT, PB11_STATUS_OUT_OF_RANGE, &
       PB11_STATUS_NUMERICAL_FAILURE, &
       fusion_birth_table_control_v1, fusion_birth_table_info_v1, &
       fusion_birth_coefficients_v1, fusion_birth_table_create, &
       fusion_birth_table_destroy, fusion_birth_table_info, &
       fusion_birth_table_evaluate
  implicit none

  real(c_double), parameter :: kev_j = 1.602176634e-16_c_double
  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  real(c_double), parameter :: lower_kT_j = 20.0_c_double * kev_j
  real(c_double), parameter :: upper_kT_j = 20.2_c_double * kev_j
  real(c_double), parameter :: intermediate_kT_j = 20.1_c_double * kev_j
  integer(c_int), parameter :: cells = 100_c_int
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_table_lifecycle()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran birth-table test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran birth-table tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function close_scaled(actual, expected, relative, absolute)
    real(c_double), intent(in) :: actual, expected, relative, absolute
    real(c_double) :: scale

    if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
       close_scaled = .false.
       return
    end if
    scale = max(abs(actual), abs(expected))
    close_scaled = abs(actual - expected) <= absolute + relative * scale
  end function close_scaled

  subroutine make_source(source)
    type(fusion_thermal_birth_options_v1), intent(out) :: source

    source%relative_max_J = 2.5_c_double * mev_j
    source%cm_max_kT = 40.0_c_double
    source%ground_state_q_J = 91.84_c_double * kev_j
    ! Use the source API's canonical one-keV endpoint spelling.
    source%cutoff_J = 0.001_c_double * mev_j
    source%l1_fraction = 0.76_c_double
    source%relative_phase = 0.0_c_double
    source%narrow_peak_fraction = 0.051_c_double
    source%continuum_peak_scale = 1.0_c_double
    source%continuation = FUSION_ENDPOINT_S
    source%pb_low = FUSION_PB_LOW_TB
    source%remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY
    source%broad_mode = 13_c_int
    source%fsci_policy = 0_c_int
    source%relative_order = 8_c_int
    source%cm_order = 8_c_int
    source%nq = 4_c_int
    source%ncos = 4_c_int
  end subroutine make_source

  subroutine make_control(control)
    type(fusion_birth_table_control_v1), intent(out) :: control

    control%max_rate_error = 0.03_c_double
    control%max_debit_error = 0.03_c_double
    control%max_number_L1 = 0.03_c_double
    control%max_energy_L1 = 0.03_c_double
    control%max_direct_rate_discrepancy = 0.001_c_double
    control%max_direct_debit_discrepancy = 0.001_c_double
    control%max_knots = 16_c_int
    control%max_evaluations = 128_c_int
    control%max_depth = 8_c_int
  end subroutine make_control

  subroutine make_edges(edges)
    real(c_double), intent(out) :: edges(:)
    integer :: i

    do i = 1, size(edges)
       edges(i) = 25.0_c_double * mev_j * real(i - 1, c_double) / &
            real(size(edges) - 1, c_double)
    end do
  end subroutine make_edges

  logical function coefficients_are_zero(out)
    type(fusion_birth_coefficients_v1), intent(in) :: out

    coefficients_are_zero = out%reactivity_m3_s == 0.0_c_double .and. &
         all(out%reactant_energy_moment_J_m3_s == 0.0_c_double) .and. &
         all(out%below_number_m3_s == 0.0_c_double) .and. &
         all(out%below_energy_J_m3_s == 0.0_c_double) .and. &
         all(out%above_number_m3_s == 0.0_c_double) .and. &
         all(out%above_energy_J_m3_s == 0.0_c_double)
  end function coefficients_are_zero

  logical function coefficients_are_finite(out)
    type(fusion_birth_coefficients_v1), intent(in) :: out

    coefficients_are_finite = ieee_is_finite(out%reactivity_m3_s) .and. &
         all(ieee_is_finite(out%reactant_energy_moment_J_m3_s)) .and. &
         all(ieee_is_finite(out%below_number_m3_s)) .and. &
         all(ieee_is_finite(out%below_energy_J_m3_s)) .and. &
         all(ieee_is_finite(out%above_number_m3_s)) .and. &
         all(ieee_is_finite(out%above_energy_J_m3_s))
  end function coefficients_are_finite

  subroutine fill_coefficients(out, value)
    type(fusion_birth_coefficients_v1), intent(out) :: out
    real(c_double), intent(in) :: value

    out%reactivity_m3_s = value
    out%reactant_energy_moment_J_m3_s = value
    out%below_number_m3_s = value
    out%below_energy_J_m3_s = value
    out%above_number_m3_s = value
    out%above_energy_J_m3_s = value
  end subroutine fill_coefficients

  subroutine fill_info(info, value)
    type(fusion_birth_table_info_v1), intent(out) :: info
    real(c_double), intent(in) :: value

    info%lower_kT_J = value
    info%upper_kT_J = value
    info%max_validated_rate_error = value
    info%max_validated_debit_error = value
    info%max_validated_number_L1 = value
    info%max_validated_energy_L1 = value
    info%max_sampled_direct_rate_discrepancy = value
    info%max_sampled_direct_debit_discrepancy = value
    info%channel = int(value, c_int)
    info%cells = int(value, c_int)
    info%knots = int(value, c_int)
    info%direct_evaluations = int(value, c_int)
    info%source%relative_max_J = value
    info%source%cm_max_kT = value
    info%source%ground_state_q_J = value
    info%source%cutoff_J = value
    info%source%l1_fraction = value
    info%source%relative_phase = value
    info%source%narrow_peak_fraction = value
    info%source%continuum_peak_scale = value
    info%source%continuation = int(value, c_int)
    info%source%pb_low = int(value, c_int)
    info%source%remainder_policy = int(value, c_int)
    info%source%broad_mode = int(value, c_int)
    info%source%fsci_policy = int(value, c_int)
    info%source%relative_order = int(value, c_int)
    info%source%cm_order = int(value, c_int)
    info%source%nq = int(value, c_int)
    info%source%ncos = int(value, c_int)
    info%control%max_rate_error = value
    info%control%max_debit_error = value
    info%control%max_number_L1 = value
    info%control%max_energy_L1 = value
    info%control%max_direct_rate_discrepancy = value
    info%control%max_direct_debit_discrepancy = value
    info%control%max_knots = int(value, c_int)
    info%control%max_evaluations = int(value, c_int)
    info%control%max_depth = int(value, c_int)
  end subroutine fill_info

  logical function info_is_zero(info)
    type(fusion_birth_table_info_v1), intent(in) :: info

    info_is_zero = info%lower_kT_J == 0.0_c_double .and. &
         info%upper_kT_J == 0.0_c_double .and. &
         info%max_validated_rate_error == 0.0_c_double .and. &
         info%max_validated_debit_error == 0.0_c_double .and. &
         info%max_validated_number_L1 == 0.0_c_double .and. &
         info%max_validated_energy_L1 == 0.0_c_double .and. &
         info%max_sampled_direct_rate_discrepancy == 0.0_c_double .and. &
         info%max_sampled_direct_debit_discrepancy == 0.0_c_double .and. &
         info%channel == 0_c_int .and. info%cells == 0_c_int .and. &
         info%knots == 0_c_int .and. info%direct_evaluations == 0_c_int .and. &
         info%source%relative_max_J == 0.0_c_double .and. &
         info%source%cm_max_kT == 0.0_c_double .and. &
         info%source%ground_state_q_J == 0.0_c_double .and. &
         info%source%cutoff_J == 0.0_c_double .and. &
         info%source%l1_fraction == 0.0_c_double .and. &
         info%source%relative_phase == 0.0_c_double .and. &
         info%source%narrow_peak_fraction == 0.0_c_double .and. &
         info%source%continuum_peak_scale == 0.0_c_double .and. &
         info%source%continuation == 0_c_int .and. &
         info%source%pb_low == 0_c_int .and. &
         info%source%remainder_policy == 0_c_int .and. &
         info%source%broad_mode == 0_c_int .and. &
         info%source%fsci_policy == 0_c_int .and. &
         info%source%relative_order == 0_c_int .and. &
         info%source%cm_order == 0_c_int .and. &
         info%source%nq == 0_c_int .and. info%source%ncos == 0_c_int .and. &
         info%control%max_rate_error == 0.0_c_double .and. &
         info%control%max_debit_error == 0.0_c_double .and. &
         info%control%max_number_L1 == 0.0_c_double .and. &
         info%control%max_energy_L1 == 0.0_c_double .and. &
         info%control%max_direct_rate_discrepancy == 0.0_c_double .and. &
         info%control%max_direct_debit_discrepancy == 0.0_c_double .and. &
         info%control%max_knots == 0_c_int .and. &
         info%control%max_evaluations == 0_c_int .and. &
         info%control%max_depth == 0_c_int
  end function info_is_zero

  subroutine verify_layout()
    type(fusion_thermal_birth_options_v1) :: source
    type(fusion_birth_table_control_v1) :: control
    type(fusion_birth_table_info_v1) :: info
    type(fusion_birth_coefficients_v1) :: coefficients

    call check(c_sizeof(source) == 104_c_size_t, &
         'nested source options retain the 104-byte ABI layout')
    call check(c_sizeof(control) == 64_c_size_t, &
         'table control has the padded 64-byte C layout')
    call check(c_sizeof(info) == 248_c_size_t, &
         'table info preserves nested source/control order and padding')
    call check(c_sizeof(coefficients) == 31_c_size_t * c_sizeof(0.0_c_double), &
         'birth coefficients contain exactly 31 c_double values')
  end subroutine verify_layout

  subroutine verify_table_lifecycle()
    real(c_double), target :: edges(cells + 1), birth(cells, 7)
    real(c_double), target :: bad_birth(cells - 1, 7), bad_species(cells, 6)
    type(fusion_thermal_birth_options_v1), target :: source
    type(fusion_birth_table_control_v1), target :: control
    type(fusion_birth_table_info_v1), target :: info
    type(fusion_birth_coefficients_v1), target :: coefficients
    type(c_ptr) :: table
    type(c_ptr) :: rejected_table
    type(fusion_thermal_birth_options_v1), target :: narrow_source
    real(c_double) :: event_rate, alpha_number, neutron_number
    integer(c_int) :: status

    call make_source(source)
    call make_control(control)
    call make_edges(edges)
    table = c_null_ptr

    ! A 100 keV source window truncates the DT reacting-energy tail at these
    ! temperatures.  The direct discrepancy gate must reject that source;
    ! refining the temperature table cannot repair a source-window error.
    narrow_source = source
    narrow_source%relative_max_J = 100.0_c_double * kev_j
    rejected_table = c_null_ptr
    call fusion_birth_table_create(FUSION_DT_ALPHAN, lower_kT_j, upper_kT_j, &
         narrow_source, control, edges, rejected_table, status)
    call check(status == PB11_STATUS_NUMERICAL_FAILURE .and. &
         .not. c_associated(rejected_table), &
         '100 keV source window is rejected by the direct-integral gate')

    call fusion_birth_table_create(FUSION_DT_ALPHAN, lower_kT_j, upper_kT_j, &
         source, control, edges, table, status)
    if (status /= PB11_STATUS_OK) write(*, '(A, I0)') &
         'birth-table create status: ', status
    call check(status == PB11_STATUS_OK .and. c_associated(table), &
         'DT narrow-range table construction succeeds')
    if (.not. c_associated(table)) return

    call fill_info(info, -1.0_c_double)
    call fusion_birth_table_info(table, info, status)
    call check(status == PB11_STATUS_OK, 'table info query succeeds')
    call check(info%channel == FUSION_DT_ALPHAN .and. info%cells == cells .and. &
         info%lower_kT_J == lower_kT_j .and. info%upper_kT_J == upper_kT_j, &
         'table info reports channel, extent, and temperature domain')
    call check(info%knots >= 2_c_int .and. info%knots <= control%max_knots .and. &
         info%direct_evaluations >= 5_c_int .and. &
         info%direct_evaluations <= control%max_evaluations, &
         'table info reports bounded construction work')
    call check(info%source%relative_order == source%relative_order .and. &
         info%source%cm_order == source%cm_order .and. &
         info%source%nq == source%nq .and. info%source%ncos == source%ncos .and. &
         info%source%cutoff_J == source%cutoff_J, &
         'table info preserves the nested source controls')
    call check(info%control%max_knots == control%max_knots .and. &
         info%control%max_evaluations == control%max_evaluations .and. &
         info%control%max_depth == control%max_depth .and. &
         info%control%max_rate_error == control%max_rate_error, &
         'table info preserves the nested table gates')
    call check(ieee_is_finite(info%max_validated_rate_error) .and. &
         ieee_is_finite(info%max_validated_debit_error) .and. &
         ieee_is_finite(info%max_validated_number_L1) .and. &
         ieee_is_finite(info%max_validated_energy_L1) .and. &
         ieee_is_finite(info%max_sampled_direct_rate_discrepancy) .and. &
         ieee_is_finite(info%max_sampled_direct_debit_discrepancy), &
         'table validation diagnostics are finite')

    birth = -1.0_c_double
    call fill_coefficients(coefficients, -1.0_c_double)
    call fusion_birth_table_evaluate(table, intermediate_kT_j, birth, &
         coefficients, status)
    call check(status == PB11_STATUS_OK, &
         'table intermediate-temperature evaluation succeeds')
    call check(coefficients_are_finite(coefficients) .and. &
         coefficients%reactivity_m3_s > 0.0_c_double, &
         'intermediate coefficients are finite and have a positive rate')
    call check(all(ieee_is_finite(birth)) .and. all(birth >= 0.0_c_double), &
         'intermediate species-major birth table is finite and nonnegative')

    event_rate = coefficients%reactivity_m3_s
    alpha_number = sum(birth(:, 5)) + coefficients%below_number_m3_s(5) + &
         coefficients%above_number_m3_s(5)
    neutron_number = sum(birth(:, 7)) + coefficients%below_number_m3_s(7) + &
         coefficients%above_number_m3_s(7)
    call check(close_scaled(alpha_number, event_rate, 2.0e-9_c_double, &
         1.0e-300_c_double), &
         'table alpha number closes across mapped and spill coefficients')
    call check(close_scaled(neutron_number, event_rate, 2.0e-9_c_double, &
         1.0e-300_c_double), &
         'table neutron number closes across mapped and spill coefficients')

    birth = 3.0_c_double
    call fill_coefficients(coefficients, 3.0_c_double)
    call fusion_birth_table_evaluate(table, upper_kT_j + 0.1_c_double * kev_j, &
         birth, coefficients, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         all(birth == 0.0_c_double) .and. coefficients_are_zero(coefficients), &
         'out-of-domain table evaluation clears every output')

    bad_birth = 4.0_c_double
    call fill_coefficients(coefficients, 4.0_c_double)
    call fusion_birth_table_evaluate(table, intermediate_kT_j, bad_birth, &
         coefficients, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(bad_birth == 0.0_c_double) .and. coefficients_are_zero(coefficients), &
         'wrong cell extent is rejected from queried table info and cleared')

    bad_species = 5.0_c_double
    call fill_coefficients(coefficients, 5.0_c_double)
    call fusion_birth_table_evaluate(table, intermediate_kT_j, bad_species, &
         coefficients, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(bad_species == 0.0_c_double) .and. coefficients_are_zero(coefficients), &
         'wrong species extent is rejected before C_LOC and cleared')

    call fusion_birth_table_destroy(table)
    call check(.not. c_associated(table), 'destroy clears the caller handle')
    call fill_info(info, -1.0_c_double)
    call fusion_birth_table_info(table, info, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. info_is_zero(info), &
         'query through a destroyed handle clears the info output')
  end subroutine verify_table_lifecycle

end program test_fusion_birth_table_fortran
