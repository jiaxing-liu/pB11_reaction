program test_fusion_fsci_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_spectrum_fortran, only : fusion_nuclear_coulomb_v1, &
       fusion_alpha_amplitudes_v1, fusion_alpha_spectrum_v1, &
       FUSION_ALPHA_ALPHA_L2, fusion_alpha_amplitudes_cutoff, &
       fusion_alpha_spectrum_grid
  use fusion_fsci_fortran
  implicit none

  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  integer, parameter :: cells = 20
  integer(c_int), parameter :: nq = 16_c_int, nc = 16_c_int
  real(c_double), parameter :: available = 8.84_c_double * mev_j
  real(c_double), parameter :: intermediate = 3.129_c_double * mev_j
  real(c_double), parameter :: cutoff = 0.001_c_double * mev_j
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_radius16()
  call verify_policy0_parity()
  call verify_policy1_model()
  call verify_errors_and_clearing()
  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran FSCI test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran FSCI tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function close_value(actual, expected, tolerance)
    real(c_double), intent(in) :: actual, expected, tolerance
    real(c_double) :: scale
    if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
       close_value = .false.
    else if (actual == expected) then
       close_value = .true.
    else
       scale = max(abs(actual), abs(expected))
       close_value = scale > 0.0_c_double .and. &
            abs(actual - expected) <= tolerance * scale
    end if
  end function close_value

  logical function nuclear_is_zero(out)
    type(fusion_nuclear_coulomb_v1), intent(in) :: out
    nuclear_is_zero = all([out%log_penetrability, out%shift, out%phase_real, &
         out%phase_imag, out%rho] == 0.0_c_double)
  end function nuclear_is_zero

  logical function amplitudes_are_zero(out)
    type(fusion_alpha_amplitudes_v1), intent(in) :: out
    amplitudes_are_zero = all(out%unsym_real == 0.0_c_double) .and. &
         all(out%unsym_imag == 0.0_c_double) .and. &
         all(out%sym_real == 0.0_c_double) .and. &
         all(out%sym_imag == 0.0_c_double) .and. &
         out%phase_space_J == 0.0_c_double
  end function amplitudes_are_zero

  logical function spectrum_is_zero(out)
    type(fusion_alpha_spectrum_v1), intent(in) :: out
    spectrum_is_zero = all([out%mapped_number, out%mapped_energy_J, &
         out%below_number, out%below_energy_J, out%above_number, &
         out%above_energy_J, out%number_residual, out%energy_residual_J, &
         out%normalization_J2, out%l1_normalization_J2, &
         out%l3_normalization_J2] == 0.0_c_double) .and. &
         out%quadrature_events == 0_c_int .and. out%pruned_events == 0_c_int
  end function spectrum_is_zero

  subroutine make_edges(edges)
    real(c_double), intent(out) :: edges(:)
    integer :: i
    do i = 1, size(edges)
       edges(i) = 10.0_c_double * mev_j * real(i - 1, c_double) / &
            real(size(edges) - 1, c_double)
    end do
  end subroutine make_edges

  subroutine verify_layout()
    type(fusion_nuclear_coulomb_v1) :: coulomb
    type(fusion_alpha_amplitudes_v1) :: amplitudes
    type(fusion_alpha_spectrum_v1) :: spectrum
    real(c_double) :: d
    integer(c_int) :: j
    call check(c_sizeof(coulomb) == 5_c_size_t * c_sizeof(d), &
         'FSCI Coulomb layout is five doubles')
    call check(c_sizeof(amplitudes) == 21_c_size_t * c_sizeof(d), &
         'FSCI amplitude layout is twenty-one doubles')
    call check(c_sizeof(spectrum) == 11_c_size_t * c_sizeof(d) + &
         2_c_size_t * c_sizeof(j), 'FSCI spectrum layout is unchanged')
  end subroutine verify_layout

  subroutine verify_radius16()
    type(fusion_nuclear_coulomb_v1) :: out
    real(c_double) :: phase_norm
    integer(c_int) :: status
    call fusion_nuclear_coulomb_radius16(FUSION_ALPHA_ALPHA_L2, &
         intermediate, out, status)
    phase_norm = hypot(out%phase_real, out%phase_imag)
    call check(status == PB11_STATUS_OK .and. &
         all(ieee_is_finite([out%log_penetrability, out%shift, &
         out%phase_real, out%phase_imag, out%rho])) .and. out%rho > 0.0_c_double &
         .and. abs(phase_norm - 1.0_c_double) < 5.0e-6_c_double, &
         'radius-16 wrapper returns a finite unit phase')
  end subroutine verify_radius16

  subroutine verify_policy0_parity()
    real(c_double) :: edges(cells + 1), old_birth(cells), new_birth(cells)
    type(fusion_alpha_amplitudes_v1) :: old_a, new_a
    type(fusion_alpha_spectrum_v1) :: old_s, new_s
    integer(c_int) :: old_pruned, new_pruned, old_status, new_status
    call make_edges(edges)
    call fusion_alpha_amplitudes_cutoff(2_c_int, available, intermediate, &
         0.0_c_double, cutoff, old_a, old_pruned, old_status)
    call fusion_alpha_amplitudes_fsci_cutoff(2_c_int, FUSION_ALPHA_FSCI_NONE, &
         available, intermediate, 0.0_c_double, cutoff, new_a, new_pruned, &
         new_status)
    call check(old_status == PB11_STATUS_OK .and. &
         new_status == PB11_STATUS_OK .and. old_pruned == new_pruned .and. &
         maxval(abs(old_a%sym_real - new_a%sym_real)) == 0.0_c_double .and. &
         maxval(abs(old_a%sym_imag - new_a%sym_imag)) == 0.0_c_double, &
         'policy-0 amplitudes preserve the existing API exactly')

    old_birth = -1.0_c_double
    new_birth = -1.0_c_double
    call fusion_alpha_spectrum_grid(2_c_int, available, cutoff, 0.0_c_double, &
         0.0_c_double, nq, nc, cells, edges, old_birth, old_s, old_status)
    call fusion_alpha_spectrum_model_grid(2_c_int, FUSION_ALPHA_FSCI_NONE, &
         available, cutoff, 0.0_c_double, 0.0_c_double, nq, nc, cells, edges, &
         new_birth, new_s, new_status)
    call check(old_status == PB11_STATUS_OK .and. new_status == PB11_STATUS_OK &
         .and. maxval(abs(old_birth - new_birth)) == 0.0_c_double .and. &
         close_value(old_s%mapped_number, new_s%mapped_number, 1.0e-12_c_double) &
         .and. close_value(old_s%mapped_energy_J, new_s%mapped_energy_J, &
         1.0e-12_c_double) .and. old_s%quadrature_events == &
         new_s%quadrature_events .and. old_s%pruned_events == new_s%pruned_events, &
         'policy-0 spectrum preserves the existing grid API')
  end subroutine verify_policy0_parity

  subroutine verify_policy1_model()
    real(c_double) :: edges(cells + 1), birth(cells), mapped_number, mapped_energy
    type(fusion_alpha_amplitudes_v1) :: amplitudes
    type(fusion_alpha_spectrum_v1) :: spectrum
    integer(c_int) :: pruned, status
    call make_edges(edges)
    call fusion_alpha_amplitudes_fsci_cutoff(2_c_int, &
         FUSION_ALPHA_FSCI_REFS2018_R16, available, intermediate, 0.0_c_double, &
         cutoff, amplitudes, pruned, status)
    call check(status == PB11_STATUS_OK .and. &
         all(ieee_is_finite(amplitudes%sym_real)) .and. &
         all(ieee_is_finite(amplitudes%sym_imag)) .and. &
         pruned >= 0_c_int .and. pruned <= 3_c_int, &
         'policy-1 amplitudes are finite with a valid prune count')
    birth = -1.0_c_double
    call fusion_alpha_spectrum_model_grid(2_c_int, FUSION_ALPHA_FSCI_REFS2018_R16, &
         available, cutoff, 0.0_c_double, 0.0_c_double, nq, nc, cells, edges, &
         birth, spectrum, status)
    mapped_number = sum(birth)
    mapped_energy = sum(birth * 0.5_c_double * &
         (edges(1:cells) + edges(2:cells + 1)))
    call check(status == PB11_STATUS_OK .and. all(birth >= 0.0_c_double) .and. &
         all(ieee_is_finite(birth)) .and. &
         all(ieee_is_finite([spectrum%mapped_number, spectrum%mapped_energy_J, &
         spectrum%below_number, spectrum%below_energy_J, spectrum%above_number, &
         spectrum%above_energy_J, spectrum%number_residual, &
         spectrum%energy_residual_J, spectrum%normalization_J2, &
         spectrum%l1_normalization_J2, spectrum%l3_normalization_J2])) .and. &
         abs(mapped_number - spectrum%mapped_number) < 4.0e-12_c_double .and. &
         abs(mapped_energy - spectrum%mapped_energy_J) < 3.0e-11_c_double * available &
         .and. abs(mapped_number + spectrum%below_number + spectrum%above_number - &
         3.0_c_double) < 4.0e-12_c_double .and. &
         abs(mapped_energy + spectrum%below_energy_J + spectrum%above_energy_J - &
         available) < 3.0e-11_c_double * available, &
         'policy-1 spectrum is finite and closes number/energy including spill')
  end subroutine verify_policy1_model

  subroutine verify_errors_and_clearing()
    real(c_double) :: edges(cells + 1), short_edges(cells), birth(cells)
    type(fusion_nuclear_coulomb_v1) :: coulomb
    type(fusion_alpha_amplitudes_v1) :: amplitudes
    type(fusion_alpha_spectrum_v1) :: spectrum
    integer(c_int) :: pruned, status
    call make_edges(edges)
    call make_edges(short_edges)
    call fusion_nuclear_coulomb_radius16(-1_c_int, intermediate, coulomb, status)
    call check(status /= PB11_STATUS_OK .and. nuclear_is_zero(coulomb), &
         'invalid radius-16 call clears output')
    amplitudes%unsym_real = -1.0_c_double
    amplitudes%unsym_imag = -1.0_c_double
    amplitudes%sym_real = -1.0_c_double
    amplitudes%sym_imag = -1.0_c_double
    amplitudes%phase_space_J = -1.0_c_double
    pruned = -1_c_int
    call fusion_alpha_amplitudes_fsci_cutoff(2_c_int, 99_c_int, available, &
         intermediate, 0.0_c_double, cutoff, amplitudes, pruned, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. pruned == 0_c_int &
         .and. amplitudes_are_zero(amplitudes), &
         'invalid amplitude policy clears output and prune count')
    birth = 1.0_c_double
    spectrum%mapped_number = -1.0_c_double
    call fusion_alpha_spectrum_model_grid(2_c_int, 99_c_int, available, cutoff, &
         0.0_c_double, 0.0_c_double, nq, nc, cells, edges, birth, spectrum, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(spectrum), 'invalid grid policy clears outputs')
    birth = 1.0_c_double
    call fusion_alpha_spectrum_model_grid(2_c_int, FUSION_ALPHA_FSCI_NONE, available, &
         cutoff, 0.0_c_double, 0.0_c_double, nq, nc, cells, short_edges, birth, &
         spectrum, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(spectrum), 'short edge array is rejected safely')
  end subroutine verify_errors_and_clearing

end program test_fusion_fsci_fortran
