program test_fusion_thermal_birth_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_thermal_birth_fortran
  implicit none

  real(c_double), parameter :: kev_j = 1.602176634e-16_c_double
  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  real(c_double), parameter :: kT_j = 1.0_c_double * kev_j
  integer(c_int), parameter :: cells = 120_c_int
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_dt_grid()
  call verify_pB_grid()
  call verify_shape_and_policy_errors()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran thermal-birth test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran thermal-birth tests passed'

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

  subroutine clear_result_for_test(out)
    type(fusion_thermal_birth_v1), intent(out) :: out

    out%reactivity_m3_s = -7.0_c_double
    out%reference_reactivity_m3_s = -7.0_c_double
    out%reactant_energy_moment_J_m3_s = -7.0_c_double
    out%reference_reactant_energy_moment_J_m3_s = -7.0_c_double
    out%product_energy_moment_J_m3_s = -7.0_c_double
    out%number_residual_m3_s = -7.0_c_double
    out%energy_residual_J_m3_s = -7.0_c_double
    out%relative_rate_discrepancy = -7.0_c_double
    out%relative_reactant_energy_discrepancy = -7.0_c_double
    out%cm_retained_probability = -7.0_c_double
    out%cm_tail_probability = -7.0_c_double
    out%cm_tail_energy_moment_J = -7.0_c_double
    out%max_cm_energy_shift_fraction = -7.0_c_double
    out%max_shell_remap_fraction = -7.0_c_double
    out%below_number_m3_s = -7.0_c_double
    out%below_energy_J_m3_s = -7.0_c_double
    out%above_number_m3_s = -7.0_c_double
    out%above_energy_J_m3_s = -7.0_c_double
  end subroutine clear_result_for_test

  logical function result_is_zero(out)
    type(fusion_thermal_birth_v1), intent(in) :: out

    result_is_zero = out%reactivity_m3_s == 0.0_c_double .and. &
         out%reference_reactivity_m3_s == 0.0_c_double .and. &
         all(out%reactant_energy_moment_J_m3_s == 0.0_c_double) .and. &
         all(out%reference_reactant_energy_moment_J_m3_s == 0.0_c_double) .and. &
         out%product_energy_moment_J_m3_s == 0.0_c_double .and. &
         out%number_residual_m3_s == 0.0_c_double .and. &
         out%energy_residual_J_m3_s == 0.0_c_double .and. &
         out%relative_rate_discrepancy == 0.0_c_double .and. &
         out%relative_reactant_energy_discrepancy == 0.0_c_double .and. &
         out%cm_retained_probability == 0.0_c_double .and. &
         out%cm_tail_probability == 0.0_c_double .and. &
         out%cm_tail_energy_moment_J == 0.0_c_double .and. &
         out%max_cm_energy_shift_fraction == 0.0_c_double .and. &
         out%max_shell_remap_fraction == 0.0_c_double .and. &
         all(out%below_number_m3_s == 0.0_c_double) .and. &
         all(out%below_energy_J_m3_s == 0.0_c_double) .and. &
         all(out%above_number_m3_s == 0.0_c_double) .and. &
         all(out%above_energy_J_m3_s == 0.0_c_double)
  end function result_is_zero

  logical function result_is_finite(out)
    type(fusion_thermal_birth_v1), intent(in) :: out

    result_is_finite = ieee_is_finite(out%reactivity_m3_s) .and. &
         ieee_is_finite(out%reference_reactivity_m3_s) .and. &
         all(ieee_is_finite(out%reactant_energy_moment_J_m3_s)) .and. &
         all(ieee_is_finite(out%reference_reactant_energy_moment_J_m3_s)) .and. &
         ieee_is_finite(out%product_energy_moment_J_m3_s) .and. &
         ieee_is_finite(out%number_residual_m3_s) .and. &
         ieee_is_finite(out%energy_residual_J_m3_s) .and. &
         ieee_is_finite(out%relative_rate_discrepancy) .and. &
         ieee_is_finite(out%relative_reactant_energy_discrepancy) .and. &
         ieee_is_finite(out%cm_retained_probability) .and. &
         ieee_is_finite(out%cm_tail_probability) .and. &
         ieee_is_finite(out%cm_tail_energy_moment_J) .and. &
         ieee_is_finite(out%max_cm_energy_shift_fraction) .and. &
         ieee_is_finite(out%max_shell_remap_fraction) .and. &
         all(ieee_is_finite(out%below_number_m3_s)) .and. &
         all(ieee_is_finite(out%below_energy_J_m3_s)) .and. &
         all(ieee_is_finite(out%above_number_m3_s)) .and. &
         all(ieee_is_finite(out%above_energy_J_m3_s))
  end function result_is_finite

  subroutine make_options(options)
    type(fusion_thermal_birth_options_v1), intent(out) :: options

    options%relative_max_J = 100.0_c_double * kev_j
    options%cm_max_kT = 40.0_c_double
    options%ground_state_q_J = 91.84_c_double * kev_j
    ! Spell 1 keV as the existing API's canonical 0.001 MeV endpoint.  The
    ! separately rounded SI literal can land one ULP below that lower bound.
    options%cutoff_J = 0.001_c_double * mev_j
    options%l1_fraction = 0.76_c_double
    options%relative_phase = 0.0_c_double
    options%narrow_peak_fraction = 0.051_c_double
    options%continuum_peak_scale = 1.0_c_double
    options%continuation = FUSION_ENDPOINT_S
    options%pb_low = FUSION_PB_LOW_TB
    options%remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY
    options%broad_mode = 13_c_int
    options%fsci_policy = 0_c_int
    options%relative_order = 8_c_int
    options%cm_order = 8_c_int
    options%nq = 4_c_int
    options%ncos = 8_c_int
  end subroutine make_options

  subroutine make_edges(edges, upper_J)
    real(c_double), intent(out) :: edges(:)
    real(c_double), intent(in) :: upper_J
    integer :: i

    do i = 1, size(edges)
       edges(i) = upper_J * real(i - 1, c_double) / &
            real(size(edges) - 1, c_double)
    end do
  end subroutine make_edges

  subroutine verify_layout()
    type(fusion_thermal_birth_options_v1) :: options
    type(fusion_thermal_birth_v1) :: out
    real(c_double) :: dummy

    call check(c_sizeof(options) == 104_c_size_t, &
         'thermal-birth options have the 104-byte C ABI layout')
    call check(c_sizeof(out) == 352_c_size_t .and. &
         c_sizeof(out) == 44_c_size_t * c_sizeof(dummy), &
         'thermal-birth result has the 352-byte/44-double C ABI layout')
  end subroutine verify_layout

  subroutine verify_dt_grid()
    real(c_double), target :: edges(cells + 1), birth(cells, 7)
    type(fusion_thermal_birth_options_v1), target :: options
    type(fusion_thermal_birth_v1), target :: out
    real(c_double) :: mapped_product_energy, total_product_energy
    real(c_double) :: alpha_number, neutron_number, event_rate
    real(c_double) :: center
    integer(c_int) :: status
    integer :: i, species

    call make_options(options)
    call make_edges(edges, 20.0_c_double * mev_j)
    birth = -1.0_c_double
    call clear_result_for_test(out)

    call fusion_thermal_birth_grid(FUSION_DT_ALPHAN, kT_j, options, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_OK, 'DT thermal-birth grid returns OK')
    call check(result_is_finite(out), 'DT thermal-birth result is finite')
    call check(out%reactivity_m3_s > 0.0_c_double .and. &
         out%reference_reactivity_m3_s > 0.0_c_double, &
         'DT grid supplies positive actual and reference rates')
    call check(all(ieee_is_finite(birth)) .and. all(birth >= 0.0_c_double), &
         'DT grid births are finite and nonnegative')
    call check(out%cm_retained_probability >= -2.0e-12_c_double .and. &
         out%cm_retained_probability <= 1.0_c_double + 2.0e-12_c_double .and. &
         out%cm_tail_probability >= -2.0e-12_c_double .and. &
         out%cm_tail_probability <= 1.0_c_double + 2.0e-12_c_double .and. &
         close_scaled(out%cm_retained_probability + out%cm_tail_probability, &
         1.0_c_double, 2.0e-12_c_double, 2.0e-14_c_double), &
         'DT CM retained and tail diagnostics are probabilities')
    call check(all(out%reference_reactant_energy_moment_J_m3_s > 0.0_c_double), &
         'DT reference reactant-energy diagnostics are positive')

    event_rate = out%reactivity_m3_s
    alpha_number = sum(birth(:, 5)) + out%below_number_m3_s(5) + &
         out%above_number_m3_s(5)
    neutron_number = sum(birth(:, 7)) + out%below_number_m3_s(7) + &
         out%above_number_m3_s(7)
    call check(close_scaled(alpha_number, event_rate, 2.0e-9_c_double, &
         1.0e-300_c_double), &
         'DT alpha number closes across mapped and spill ledgers')
    call check(close_scaled(neutron_number, event_rate, 2.0e-9_c_double, &
         1.0e-300_c_double), &
         'DT neutron number closes across mapped and spill ledgers')
    call check(close_scaled(alpha_number + neutron_number, 2.0_c_double * &
         event_rate, 2.0e-9_c_double, 1.0e-300_c_double), &
         'DT product count closes at two products per event')

    mapped_product_energy = 0.0_c_double
    do i = 1, int(cells)
       center = 0.5_c_double * (edges(i) + edges(i + 1))
       do species = 1, 7
          mapped_product_energy = mapped_product_energy + &
               center * birth(i, species)
       end do
    end do
    total_product_energy = mapped_product_energy + &
         sum(out%below_energy_J_m3_s) + sum(out%above_energy_J_m3_s)
    call check(close_scaled(total_product_energy, &
         out%product_energy_moment_J_m3_s, 2.0e-8_c_double, 1.0e-300_c_double), &
         'DT product energy closes across mapped and spill ledgers')
    call check(ieee_is_finite(out%relative_rate_discrepancy) .and. &
         ieee_is_finite(out%relative_reactant_energy_discrepancy) .and. &
         ieee_is_finite(out%cm_tail_energy_moment_J) .and. &
         ieee_is_finite(out%max_cm_energy_shift_fraction) .and. &
         ieee_is_finite(out%max_shell_remap_fraction), &
         'DT reference and remap diagnostics are finite')
  end subroutine verify_dt_grid

  subroutine verify_pB_grid()
    integer(c_int), parameter :: small_cells = 32_c_int
    real(c_double), target :: edges(small_cells + 1), birth(small_cells, 7)
    type(fusion_thermal_birth_options_v1), target :: options
    type(fusion_thermal_birth_v1), target :: out
    integer(c_int) :: status

    call make_options(options)
    call make_edges(edges, 20.0_c_double * mev_j)
    birth = -1.0_c_double
    call clear_result_for_test(out)
    call fusion_thermal_birth_grid(FUSION_PB11_3ALPHA, kT_j, options, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_OK, 'pB thermal-birth grid returns OK')
    call check(result_is_finite(out), 'pB thermal-birth result is finite')
    call check(all(ieee_is_finite(birth)) .and. all(birth >= 0.0_c_double), &
         'pB thermal-birth births are finite and nonnegative')
  end subroutine verify_pB_grid

  subroutine verify_shape_and_policy_errors()
    real(c_double), target :: edges(cells + 1), short_edges(cells)
    real(c_double), target :: birth(cells, 7), wrong_species(cells, 6)
    type(fusion_thermal_birth_options_v1), target :: options
    type(fusion_thermal_birth_v1), target :: out
    integer(c_int) :: status

    call make_options(options)
    call make_edges(edges, 20.0_c_double * mev_j)
    call make_edges(short_edges, 20.0_c_double * mev_j)

    wrong_species = 1.0_c_double
    call clear_result_for_test(out)
    call fusion_thermal_birth_grid(FUSION_DT_ALPHAN, kT_j, options, edges, &
         wrong_species, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(wrong_species == 0.0_c_double) .and. result_is_zero(out), &
         'wrong seven-species extent is rejected before C_LOC and cleared')

    birth = 1.0_c_double
    call clear_result_for_test(out)
    call fusion_thermal_birth_grid(FUSION_DT_ALPHAN, kT_j, options, &
         short_edges, birth, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(birth == 0.0_c_double) .and. result_is_zero(out), &
         'wrong edge extent is rejected before C_LOC and cleared')

    options%remainder_policy = 99_c_int
    birth = 1.0_c_double
    call clear_result_for_test(out)
    call fusion_thermal_birth_grid(FUSION_PB11_3ALPHA, kT_j, options, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(birth == 0.0_c_double) .and. result_is_zero(out), &
         'invalid remainder policy clears every C output')
  end subroutine verify_shape_and_policy_errors

end program test_fusion_thermal_birth_fortran
