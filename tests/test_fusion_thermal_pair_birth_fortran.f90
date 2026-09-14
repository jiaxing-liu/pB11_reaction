program pair_birth_binding_test
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use fusion_thermal_birth_fortran
  implicit none

  real(c_double), parameter :: kev = 1.602176634e-16_c_double
  integer(c_int), parameter :: cells = 12_c_int
  real(c_double), target :: edges(cells + 1), old_birth(cells, 7)
  real(c_double), target :: pair_birth(cells, 7), unequal_birth(cells, 7)
  real(c_double), target :: bad_birth(cells, 6)
  type(fusion_thermal_birth_options_v1), target :: options
  type(fusion_thermal_birth_v1), target :: old_out, pair_out, unequal_out
  type(fusion_thermal_birth_v1), target :: bad_out
  integer(c_int) :: old_status, pair_status, unequal_status, bad_status
  integer :: i

  call make_options(options)
  do i = 1, cells + 1
     edges(i) = 20.0_c_double * 1.602176634e-13_c_double * &
          real(i - 1, c_double) / real(cells, c_double)
  end do

  old_birth = -1.0_c_double
  pair_birth = -1.0_c_double
  call fusion_thermal_birth_grid(3_c_int, kev, options, edges, old_birth, &
       old_out, old_status)
  call fusion_thermal_pair_birth_grid(3_c_int, kev, kev, 4_c_int, options, &
       edges, pair_birth, pair_out, pair_status)
  call require(old_status == PB11_STATUS_OK .and. pair_status == PB11_STATUS_OK, &
       'equal-temperature wrappers return OK')
  call require(maxval(abs(old_birth - pair_birth)) == 0.0_c_double, &
       'equal-temperature birth grids have parity')
  call require(old_out%reactivity_m3_s == pair_out%reactivity_m3_s, &
       'equal-temperature rates have parity')

  unequal_birth = -1.0_c_double
  call fusion_thermal_pair_birth_grid(3_c_int, kev, 2.0_c_double * kev, &
       8_c_int, options, edges, unequal_birth, unequal_out, &
       unequal_status)
  call require(unequal_status == PB11_STATUS_OK .and. &
       unequal_out%reactivity_m3_s > 0.0_c_double, &
       'unequal-temperature DT wrapper returns a positive result')

  bad_birth = 1.0_c_double
  call fill_bad_result(bad_out)
  call fusion_thermal_pair_birth_grid(3_c_int, kev, kev, 4_c_int, options, &
       edges, bad_birth, bad_out, bad_status)
  call require(bad_status == PB11_STATUS_INVALID_ARGUMENT .and. &
       all(bad_birth == 0.0_c_double) .and. result_is_zero(bad_out), &
       'malformed birth extent is rejected and cleared')

  write(*, '(A)') 'pair birth Fortran binding smoke test passed'

contains

  subroutine require(ok, message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    if (.not. ok) then
       write(*, '(A)') 'FAIL: ' // message
       error stop 1
    end if
  end subroutine require

  subroutine make_options(o)
    type(fusion_thermal_birth_options_v1), intent(out) :: o
    o%relative_max_J = 100.0_c_double * kev
    o%cm_max_kT = 40.0_c_double
    o%ground_state_q_J = 91.84_c_double * kev
    o%cutoff_J = 0.001_c_double * 1.602176634e-13_c_double
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
    o%cm_order = 4_c_int
    o%nq = 4_c_int
    o%ncos = 4_c_int
  end subroutine make_options

  subroutine fill_bad_result(o)
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
  end subroutine fill_bad_result

  logical function result_is_zero(o)
    type(fusion_thermal_birth_v1), intent(in) :: o
    result_is_zero = o%reactivity_m3_s == 0.0_c_double .and. &
         o%reference_reactivity_m3_s == 0.0_c_double .and. &
         all(o%reactant_energy_moment_J_m3_s == 0.0_c_double) .and. &
         all(o%reference_reactant_energy_moment_J_m3_s == 0.0_c_double) .and. &
         o%product_energy_moment_J_m3_s == 0.0_c_double .and. &
         o%number_residual_m3_s == 0.0_c_double .and. &
         o%energy_residual_J_m3_s == 0.0_c_double .and. &
         o%relative_rate_discrepancy == 0.0_c_double .and. &
         o%relative_reactant_energy_discrepancy == 0.0_c_double .and. &
         o%cm_retained_probability == 0.0_c_double .and. &
         o%cm_tail_probability == 0.0_c_double .and. &
         o%cm_tail_energy_moment_J == 0.0_c_double .and. &
         o%max_cm_energy_shift_fraction == 0.0_c_double .and. &
         o%max_shell_remap_fraction == 0.0_c_double .and. &
         all(o%below_number_m3_s == 0.0_c_double) .and. &
         all(o%below_energy_J_m3_s == 0.0_c_double) .and. &
         all(o%above_number_m3_s == 0.0_c_double) .and. &
         all(o%above_energy_J_m3_s == 0.0_c_double)
  end function result_is_zero

end program pair_birth_binding_test
