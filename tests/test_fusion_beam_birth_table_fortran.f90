program test_fusion_beam_birth_table_fortran
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_int64_t, c_null_ptr, c_ptr, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_beam_birth_table_fortran, only : &
       PB11_STATUS_OK, PB11_STATUS_INVALID_ARGUMENT, &
       PB11_STATUS_OUT_OF_RANGE, FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, &
       FUSION_PB_LOW_TB, FUSION_PB_REMAINDER_ENTRANCE_PROXY, &
       FUSION_BEAM_BIRTH_TABLE_SPECIES, fusion_beam_birth_options_v1, &
       fusion_birth_table_control_v1, fusion_beam_birth_table_info_v1, &
       fusion_birth_coefficients_v1, fusion_beam_birth_table_create, &
       fusion_beam_birth_table_destroy, fusion_beam_birth_table_info, &
       fusion_beam_birth_table_evaluate
  implicit none

  real(c_double), parameter :: kev_j = 1.602176634e-16_c_double
  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  real(c_double), parameter :: projectile_energy_j = 100.0_c_double * kev_j
  real(c_double), parameter :: lower_kT_j = 0.08_c_double * kev_j
  real(c_double), parameter :: upper_kT_j = 0.12_c_double * kev_j
  integer(c_int), parameter :: cells = 16_c_int
  real(c_double), target :: edges(cells + 1_c_int)
  real(c_double), target :: birth(cells, FUSION_BEAM_BIRTH_TABLE_SPECIES)
  type(fusion_beam_birth_options_v1), target :: source
  type(fusion_birth_table_control_v1), target :: control
  type(fusion_beam_birth_table_info_v1), target :: info
  type(fusion_birth_coefficients_v1), target :: coefficients
  type(c_ptr) :: table
  integer(c_int) :: status
  integer :: failures

  failures = 0
  call verify_layout()
  call make_source(source)
  call make_control(control)
  call make_edges(edges)
  table = c_null_ptr

  call fusion_beam_birth_table_create(FUSION_DT_ALPHAN, 0_c_int, &
       projectile_energy_j, lower_kT_j, upper_kT_j, source, control, edges, &
       table, status)
  call check(status == PB11_STATUS_OK .and. c_associated(table), &
       'fixed-energy DT table construction succeeds')
  if (.not. c_associated(table)) then
     write(*, '(A, I0)') 'create status: ', status
  else
     call verify_info()
     call verify_endpoints()
     call verify_out_of_range()
     call fusion_beam_birth_table_destroy(table)
     call check(.not. c_associated(table), 'destroy clears the caller handle')
     call fill_info(info, -1.0_c_double)
     call fusion_beam_birth_table_info(table, info, status)
     call check(status == PB11_STATUS_INVALID_ARGUMENT .and. info_is_zero(info), &
          'query through a destroyed handle clears info output')
  end if

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran beam-birth-table test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'Fortran beam-birth-table smoke test passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  subroutine verify_layout()
    type(fusion_beam_birth_options_v1) :: local_source
    type(fusion_birth_table_control_v1) :: local_control
    type(fusion_beam_birth_table_info_v1) :: local_info
    type(fusion_birth_coefficients_v1) :: local_coefficients

    call check(c_sizeof(local_source) == 104_c_size_t, &
         'beam source options retain the 104-byte C ABI layout')
    call check(c_sizeof(local_control) == 64_c_size_t, &
         'table control retains the 64-byte C ABI layout')
    call check(c_sizeof(local_info) == 280_c_size_t, &
         'beam table info matches the nested C ABI layout')
    call check(c_sizeof(local_coefficients) == &
         31_c_size_t * c_sizeof(0.0_c_double), &
         'birth coefficients are reused as 31 c_double values')
  end subroutine verify_layout

  subroutine make_source(o)
    type(fusion_beam_birth_options_v1), intent(out) :: o

    o%relative_max_J = 2.5_c_double * mev_j
    o%angular_max_exponent = 40.0_c_double
    o%ground_state_q_J = 91.84_c_double * kev_j
    ! Keep the one-keV source cutoff as the exact requested expression.
    o%cutoff_J = 0.001_c_double * mev_j
    o%l1_fraction = 0.76_c_double
    o%relative_phase = 0.0_c_double
    o%narrow_peak_fraction = 0.051_c_double
    o%continuum_peak_scale = 1.0_c_double
    o%continuation = FUSION_ENDPOINT_S
    o%pb_low = FUSION_PB_LOW_TB
    o%remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY
    o%broad_mode = 13_c_int
    o%fsci_policy = 0_c_int
    o%relative_order = 16_c_int
    o%angular_order = 8_c_int
    o%nq = 4_c_int
    o%ncos = 4_c_int
  end subroutine make_source

  subroutine make_control(o)
    type(fusion_birth_table_control_v1), intent(out) :: o

    o%max_rate_error = 0.01_c_double
    o%max_debit_error = 0.01_c_double
    o%max_number_L1 = 0.03_c_double
    o%max_energy_L1 = 0.03_c_double
    o%max_direct_rate_discrepancy = 1.0e-5_c_double
    o%max_direct_debit_discrepancy = 1.0e-5_c_double
    o%max_knots = 256_c_int
    o%max_evaluations = 2048_c_int
    o%max_depth = 16_c_int
  end subroutine make_control

  subroutine make_edges(e)
    real(c_double), intent(out) :: e(:)
    integer :: i

    do i = 1, size(e)
       e(i) = 25.0_c_double * mev_j * real(i - 1, c_double) / &
            real(size(e) - 1, c_double)
    end do
  end subroutine make_edges

  subroutine verify_info()
    call fill_info(info, -1.0_c_double)
    call fusion_beam_birth_table_info(table, info, status)
    call check(status == PB11_STATUS_OK, 'table info query succeeds')
    call check(info%projectile_energy_J == projectile_energy_j .and. &
         info%lower_kT_J == lower_kT_j .and. info%upper_kT_J == upper_kT_j .and. &
         info%channel == FUSION_DT_ALPHAN .and. info%projectile_slot == 0_c_int .and. &
         info%cells == cells, 'info reports fixed energy, domain, and extent')
    call check(info%knots >= 2_c_int .and. info%knots <= control%max_knots .and. &
         info%direct_evaluations >= 5_c_int .and. &
         info%direct_evaluations <= control%max_evaluations, &
         'info reports bounded table construction work')
    call check(info%spectral_entries_evaluated > 0_c_int64_t .and. &
         info%stored_spectral_entries > 0_c_int64_t, &
         'info reports nonzero int64 spectral counters')
    call check(info%source%angular_order == 8_c_int .and. &
         info%source%relative_order == 16_c_int .and. &
         info%source%cutoff_J == 0.001_c_double * mev_j .and. &
         info%control%max_rate_error == 0.01_c_double .and. &
         info%control%max_debit_error == 0.01_c_double .and. &
         info%control%max_number_L1 == 0.03_c_double .and. &
         info%control%max_energy_L1 == 0.03_c_double .and. &
         info%control%max_direct_rate_discrepancy == 1.0e-5_c_double .and. &
         info%control%max_direct_debit_discrepancy == 1.0e-5_c_double .and. &
         info%control%max_knots == 256_c_int .and. &
         info%control%max_evaluations == 2048_c_int .and. &
         info%control%max_depth == 16_c_int, &
         'info preserves source and construction controls')
    call check(ieee_is_finite(info%max_validated_rate_error) .and. &
         ieee_is_finite(info%max_validated_debit_error) .and. &
         ieee_is_finite(info%max_validated_number_L1) .and. &
         ieee_is_finite(info%max_validated_energy_L1) .and. &
         ieee_is_finite(info%max_sampled_direct_rate_discrepancy) .and. &
         ieee_is_finite(info%max_sampled_direct_debit_discrepancy), &
         'info validation diagnostics are finite')
  end subroutine verify_info

  subroutine verify_endpoints()
    real(c_double), parameter :: endpoint_values(2) = [lower_kT_j, upper_kT_j]
    integer :: i

    do i = 1, size(endpoint_values)
       birth = -1.0_c_double
       call fill_coefficients(coefficients, -1.0_c_double)
       call fusion_beam_birth_table_evaluate(table, endpoint_values(i), birth, &
            coefficients, status)
       call check(status == PB11_STATUS_OK, 'table endpoint evaluation succeeds')
       call check(ieee_is_finite(coefficients%reactivity_m3_s) .and. &
            coefficients%reactivity_m3_s > 0.0_c_double .and. &
            all(ieee_is_finite(birth)) .and. all(birth >= 0.0_c_double), &
            'endpoint output is finite, positive-rate, and nonnegative')
    end do
  end subroutine verify_endpoints

  subroutine verify_out_of_range()
    birth = 3.0_c_double
    call fill_coefficients(coefficients, 3.0_c_double)
    call fusion_beam_birth_table_evaluate(table, upper_kT_j + 0.001_c_double * &
         kev_j, birth, coefficients, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         all(birth == 0.0_c_double) .and. coefficients_are_zero(coefficients), &
         'out-of-domain evaluation clears dense birth and coefficients')
  end subroutine verify_out_of_range

  subroutine fill_coefficients(o, value)
    type(fusion_birth_coefficients_v1), intent(out) :: o
    real(c_double), intent(in) :: value

    o%reactivity_m3_s = value
    o%reactant_energy_moment_J_m3_s = value
    o%below_number_m3_s = value
    o%below_energy_J_m3_s = value
    o%above_number_m3_s = value
    o%above_energy_J_m3_s = value
  end subroutine fill_coefficients

  logical function coefficients_are_zero(o)
    type(fusion_birth_coefficients_v1), intent(in) :: o

    coefficients_are_zero = o%reactivity_m3_s == 0.0_c_double .and. &
         all(o%reactant_energy_moment_J_m3_s == 0.0_c_double) .and. &
         all(o%below_number_m3_s == 0.0_c_double) .and. &
         all(o%below_energy_J_m3_s == 0.0_c_double) .and. &
         all(o%above_number_m3_s == 0.0_c_double) .and. &
         all(o%above_energy_J_m3_s == 0.0_c_double)
  end function coefficients_are_zero

  subroutine fill_info(o, value)
    type(fusion_beam_birth_table_info_v1), intent(out) :: o
    real(c_double), intent(in) :: value

    o%projectile_energy_J = value
    o%lower_kT_J = value
    o%upper_kT_J = value
    o%max_validated_rate_error = value
    o%max_validated_debit_error = value
    o%max_validated_number_L1 = value
    o%max_validated_energy_L1 = value
    o%max_sampled_direct_rate_discrepancy = value
    o%max_sampled_direct_debit_discrepancy = value
    o%channel = int(value, c_int)
    o%projectile_slot = int(value, c_int)
    o%cells = int(value, c_int)
    o%knots = int(value, c_int)
    o%direct_evaluations = int(value, c_int)
    o%spectral_entries_evaluated = int(value, c_int64_t)
    o%stored_spectral_entries = int(value, c_int64_t)
    o%source%relative_max_J = value
    o%source%angular_max_exponent = value
    o%source%ground_state_q_J = value
    o%source%cutoff_J = value
    o%source%l1_fraction = value
    o%source%relative_phase = value
    o%source%narrow_peak_fraction = value
    o%source%continuum_peak_scale = value
    o%source%continuation = int(value, c_int)
    o%source%pb_low = int(value, c_int)
    o%source%remainder_policy = int(value, c_int)
    o%source%broad_mode = int(value, c_int)
    o%source%fsci_policy = int(value, c_int)
    o%source%relative_order = int(value, c_int)
    o%source%angular_order = int(value, c_int)
    o%source%nq = int(value, c_int)
    o%source%ncos = int(value, c_int)
    o%control%max_rate_error = value
    o%control%max_debit_error = value
    o%control%max_number_L1 = value
    o%control%max_energy_L1 = value
    o%control%max_direct_rate_discrepancy = value
    o%control%max_direct_debit_discrepancy = value
    o%control%max_knots = int(value, c_int)
    o%control%max_evaluations = int(value, c_int)
    o%control%max_depth = int(value, c_int)
  end subroutine fill_info

  logical function info_is_zero(o)
    type(fusion_beam_birth_table_info_v1), intent(in) :: o

    info_is_zero = o%projectile_energy_J == 0.0_c_double .and. &
         o%lower_kT_J == 0.0_c_double .and. o%upper_kT_J == 0.0_c_double .and. &
         o%max_validated_rate_error == 0.0_c_double .and. &
         o%max_validated_debit_error == 0.0_c_double .and. &
         o%max_validated_number_L1 == 0.0_c_double .and. &
         o%max_validated_energy_L1 == 0.0_c_double .and. &
         o%max_sampled_direct_rate_discrepancy == 0.0_c_double .and. &
         o%max_sampled_direct_debit_discrepancy == 0.0_c_double .and. &
         o%channel == 0_c_int .and. o%projectile_slot == 0_c_int .and. &
         o%cells == 0_c_int .and. o%knots == 0_c_int .and. &
         o%direct_evaluations == 0_c_int .and. &
         o%spectral_entries_evaluated == 0_c_int64_t .and. &
         o%stored_spectral_entries == 0_c_int64_t .and. &
         o%source%relative_max_J == 0.0_c_double .and. &
         o%source%angular_max_exponent == 0.0_c_double .and. &
         o%source%ground_state_q_J == 0.0_c_double .and. &
         o%source%cutoff_J == 0.0_c_double .and. &
         o%source%l1_fraction == 0.0_c_double .and. &
         o%source%relative_phase == 0.0_c_double .and. &
         o%source%narrow_peak_fraction == 0.0_c_double .and. &
         o%source%continuum_peak_scale == 0.0_c_double .and. &
         o%source%continuation == 0_c_int .and. o%source%pb_low == 0_c_int .and. &
         o%source%remainder_policy == 0_c_int .and. &
         o%source%broad_mode == 0_c_int .and. o%source%fsci_policy == 0_c_int .and. &
         o%source%relative_order == 0_c_int .and. &
         o%source%angular_order == 0_c_int .and. o%source%nq == 0_c_int .and. &
         o%source%ncos == 0_c_int .and. &
         o%control%max_rate_error == 0.0_c_double .and. &
         o%control%max_debit_error == 0.0_c_double .and. &
         o%control%max_number_L1 == 0.0_c_double .and. &
         o%control%max_energy_L1 == 0.0_c_double .and. &
         o%control%max_direct_rate_discrepancy == 0.0_c_double .and. &
         o%control%max_direct_debit_discrepancy == 0.0_c_double .and. &
         o%control%max_knots == 0_c_int .and. &
         o%control%max_evaluations == 0_c_int .and. &
         o%control%max_depth == 0_c_int
  end function info_is_zero

end program test_fusion_beam_birth_table_fortran
