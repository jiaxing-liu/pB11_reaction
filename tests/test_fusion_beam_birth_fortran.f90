program beam_birth_binding_test
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fusion_beam_birth_fortran
  implicit none

  real(c_double), parameter :: kev = 1.602176634e-16_c_double
  real(c_double), parameter :: mev = 1.602176634e-13_c_double
  integer(c_int), parameter :: cells = 8_c_int
  real(c_double), target :: edges(cells + 1)
  real(c_double), target :: birth_slot0(cells, FUSION_BEAM_BIRTH_SPECIES)
  real(c_double), target :: birth_slot1(cells, FUSION_BEAM_BIRTH_SPECIES)
  real(c_double), allocatable, target :: zero_birth(:,:), zero_edges(:)
  type(fusion_beam_birth_options_v1), target :: options
  type(fusion_beam_birth_options_v1), target :: base_support, resolved_support
  type(fusion_beam_birth_v1), target :: out_slot0, out_slot1, bad_out
  integer(c_int) :: status_slot0, status_slot1, bad_status
  integer :: i

  call require(c_sizeof(options) == 104_c_size_t, &
       'beam options have the 104-byte C ABI layout')
  call require(c_sizeof(out_slot0) == 376_c_size_t, &
       'beam result has the 376-byte C ABI layout')

  call make_options(options)
  do i = 1, cells + 1
     edges(i) = 20.0_c_double * mev * real(i - 1, c_double) / &
          real(cells, c_double)
  end do

  birth_slot0 = -1.0_c_double
  call fusion_beam_birth_grid(FUSION_DT_ALPHAN, 0_c_int, 100.0_c_double * kev, &
       0.0_c_double, options, edges, birth_slot0, out_slot0, status_slot0)
  call require(status_slot0 == PB11_STATUS_OK .and. &
       out_slot0%spectrum%reactivity_m3_s > 0.0_c_double .and. &
       out_slot0%spectrum%reference_reactivity_m3_s > 0.0_c_double, &
       'cold DT beam source succeeds for canonical projectile slot 0')
  call require(out_slot0%spectrum%reference_reactant_energy_moment_J_m3_s(1) &
       > 0.0_c_double .and. &
       out_slot0%spectrum%reference_reactant_energy_moment_J_m3_s(2) == &
       0.0_c_double, 'slot 0 reference debit follows canonical reactant order')

  birth_slot1 = -1.0_c_double
  call fusion_beam_birth_grid(FUSION_DT_ALPHAN, 1_c_int, 100.0_c_double * kev, &
       0.0_c_double, options, edges, birth_slot1, out_slot1, status_slot1)
  call require(status_slot1 == PB11_STATUS_OK .and. &
       out_slot1%spectrum%reactivity_m3_s > 0.0_c_double .and. &
       out_slot1%spectrum%reference_reactivity_m3_s > 0.0_c_double, &
       'cold DT beam source succeeds for canonical projectile slot 1')
  call require(out_slot1%spectrum%reference_reactant_energy_moment_J_m3_s(1) &
       == 0.0_c_double .and. &
       out_slot1%spectrum%reference_reactant_energy_moment_J_m3_s(2) > &
       0.0_c_double, 'slot 1 reference debit follows canonical reactant order')

  allocate(zero_birth(0, FUSION_BEAM_BIRTH_SPECIES), zero_edges(0))
  call fill_bad_result(bad_out)
  call fusion_beam_birth_grid(FUSION_DT_ALPHAN, 0_c_int, 100.0_c_double * kev, &
       0.0_c_double, options, zero_edges, zero_birth, bad_out, bad_status)
  call require(bad_status == PB11_STATUS_INVALID_ARGUMENT .and. &
       result_is_zero(bad_out), 'zero-cell extent is rejected and result cleared')

  ! The resolver builds only relative-energy support. It does not validate
  ! the remaining controls or certify a physical grid/table evaluation.
  base_support = options
  base_support%relative_max_J = 10.0_c_double * mev
  call fusion_beam_birth_resolve_support(30.0_c_double * mev, base_support, &
       resolved_support, bad_status)
  call require(bad_status == PB11_STATUS_OK .and. &
       resolved_support%relative_max_J == 30.0_c_double * mev .and. &
       base_support%relative_max_J == 10.0_c_double * mev .and. &
       resolved_support%angular_order == base_support%angular_order, &
       'support wrapper widens only the output cap and preserves base')
  call fusion_beam_birth_resolve_support(5.0_c_double * mev, base_support, &
       resolved_support, bad_status)
  call require(bad_status == PB11_STATUS_OK .and. &
       resolved_support%relative_max_J == 10.0_c_double * mev, &
       'support wrapper starts again from base, not prior resolved cap')

  resolved_support = base_support
  call fusion_beam_birth_resolve_support(-1.0_c_double * mev, base_support, &
       resolved_support, bad_status)
  call require(bad_status == PB11_STATUS_OUT_OF_RANGE .and. &
       options_is_zero(resolved_support), &
       'negative projectile energy rejects and clears output')
  resolved_support = base_support
  call fusion_beam_birth_resolve_support(ieee_value(0.0_c_double, &
       ieee_quiet_nan), base_support, resolved_support, bad_status)
  call require(bad_status == PB11_STATUS_INVALID_ARGUMENT .and. &
       options_is_zero(resolved_support), &
       'NaN projectile energy rejects and clears output')

  write(*, '(A)') 'beam birth Fortran binding smoke test passed'

contains

  logical function options_is_zero(o)
    type(fusion_beam_birth_options_v1), intent(in) :: o

    options_is_zero = o%relative_max_J == 0.0_c_double .and. &
         o%angular_max_exponent == 0.0_c_double .and. &
         o%ground_state_q_J == 0.0_c_double .and. &
         o%cutoff_J == 0.0_c_double .and. &
         o%l1_fraction == 0.0_c_double .and. &
         o%relative_phase == 0.0_c_double .and. &
         o%narrow_peak_fraction == 0.0_c_double .and. &
         o%continuum_peak_scale == 0.0_c_double .and. &
         o%continuation == 0_c_int .and. o%pb_low == 0_c_int .and. &
         o%remainder_policy == 0_c_int .and. o%broad_mode == 0_c_int .and. &
         o%fsci_policy == 0_c_int .and. o%relative_order == 0_c_int .and. &
         o%angular_order == 0_c_int .and. o%nq == 0_c_int .and. &
         o%ncos == 0_c_int
  end function options_is_zero

  subroutine require(ok, message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message

    if (.not. ok) then
       write(*, '(A)') 'FAIL: ' // message
       error stop 1
    end if
  end subroutine require

  subroutine make_options(o)
    type(fusion_beam_birth_options_v1), intent(out) :: o

    o%relative_max_J = 200.0_c_double * kev
    o%angular_max_exponent = 40.0_c_double
    o%ground_state_q_J = 91.84_c_double * kev
    o%cutoff_J = 0.001_c_double * mev
    o%l1_fraction = 0.76_c_double
    o%relative_phase = 0.0_c_double
    o%narrow_peak_fraction = 0.051_c_double
    o%continuum_peak_scale = 1.0_c_double
    o%continuation = FUSION_ENDPOINT_S
    o%pb_low = FUSION_PB_LOW_TB
    o%remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY
    o%broad_mode = 13_c_int
    o%fsci_policy = 0_c_int
    o%relative_order = 4_c_int
    o%angular_order = 4_c_int
    o%nq = 4_c_int
    o%ncos = 4_c_int
  end subroutine make_options

  subroutine fill_bad_result(o)
    type(fusion_beam_birth_v1), intent(out) :: o

    call fill_bad_thermal_result(o%spectrum)
    o%relative_retained_probability = -7.0_c_double
    o%retained_pair_probability = -7.0_c_double
    o%angular_omitted_pair_probability = -7.0_c_double
  end subroutine fill_bad_result

  subroutine fill_bad_thermal_result(o)
    type(fusion_thermal_birth_v1), intent(out) :: o

    o%reactivity_m3_s = -7.0_c_double
    o%reference_reactivity_m3_s = -7.0_c_double
    o%reactant_energy_moment_J_m3_s = -7.0_c_double
    o%reference_reactant_energy_moment_J_m3_s = -7.0_c_double
    o%product_energy_moment_J_m3_s = -7.0_c_double
    o%number_residual_m3_s = -7.0_c_double
    o%energy_residual_J_m3_s = -7.0_c_double
    o%relative_rate_discrepancy = -7.0_c_double
    o%relative_reactant_energy_discrepancy = -7.0_c_double
    o%cm_retained_probability = -7.0_c_double
    o%cm_tail_probability = -7.0_c_double
    o%cm_tail_energy_moment_J = -7.0_c_double
    o%max_cm_energy_shift_fraction = -7.0_c_double
    o%max_shell_remap_fraction = -7.0_c_double
    o%below_number_m3_s = -7.0_c_double
    o%below_energy_J_m3_s = -7.0_c_double
    o%above_number_m3_s = -7.0_c_double
    o%above_energy_J_m3_s = -7.0_c_double
  end subroutine fill_bad_thermal_result

  logical function result_is_zero(o)
    type(fusion_beam_birth_v1), intent(in) :: o

    result_is_zero = o%relative_retained_probability == 0.0_c_double .and. &
         o%retained_pair_probability == 0.0_c_double .and. &
         o%angular_omitted_pair_probability == 0.0_c_double .and. &
         o%spectrum%reactivity_m3_s == 0.0_c_double .and. &
         o%spectrum%reference_reactivity_m3_s == 0.0_c_double .and. &
         all(o%spectrum%reactant_energy_moment_J_m3_s == 0.0_c_double) .and. &
         all(o%spectrum%reference_reactant_energy_moment_J_m3_s == &
             0.0_c_double) .and. &
         o%spectrum%product_energy_moment_J_m3_s == 0.0_c_double .and. &
         o%spectrum%number_residual_m3_s == 0.0_c_double .and. &
         o%spectrum%energy_residual_J_m3_s == 0.0_c_double .and. &
         o%spectrum%relative_rate_discrepancy == 0.0_c_double .and. &
         o%spectrum%relative_reactant_energy_discrepancy == 0.0_c_double .and. &
         o%spectrum%cm_retained_probability == 0.0_c_double .and. &
         o%spectrum%cm_tail_probability == 0.0_c_double .and. &
         o%spectrum%cm_tail_energy_moment_J == 0.0_c_double .and. &
         o%spectrum%max_cm_energy_shift_fraction == 0.0_c_double .and. &
         o%spectrum%max_shell_remap_fraction == 0.0_c_double .and. &
         all(o%spectrum%below_number_m3_s == 0.0_c_double) .and. &
         all(o%spectrum%below_energy_J_m3_s == 0.0_c_double) .and. &
         all(o%spectrum%above_number_m3_s == 0.0_c_double) .and. &
         all(o%spectrum%above_energy_J_m3_s == 0.0_c_double)
  end function result_is_zero

end program beam_birth_binding_test
