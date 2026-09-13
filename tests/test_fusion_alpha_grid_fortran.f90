program test_fusion_alpha_grid_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_spectrum_fortran
  implicit none

  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  integer, parameter :: cells = 20
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_cutoff_wrapper()
  call verify_normalization_and_spill()
  call verify_mixture_endpoints()
  call verify_errors_and_clearing()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran alpha-grid test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran alpha-grid tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function spectrum_is_zero(out)
    type(fusion_alpha_spectrum_v1), intent(in) :: out

    spectrum_is_zero = out%mapped_number == 0.0_c_double .and. &
         out%mapped_energy_J == 0.0_c_double .and. &
         out%below_number == 0.0_c_double .and. &
         out%below_energy_J == 0.0_c_double .and. &
         out%above_number == 0.0_c_double .and. &
         out%above_energy_J == 0.0_c_double .and. &
         out%number_residual == 0.0_c_double .and. &
         out%energy_residual_J == 0.0_c_double .and. &
         out%normalization_J2 == 0.0_c_double .and. &
         out%l1_normalization_J2 == 0.0_c_double .and. &
         out%l3_normalization_J2 == 0.0_c_double .and. &
         out%quadrature_events == 0_c_int .and. &
         out%pruned_events == 0_c_int
  end function spectrum_is_zero

  logical function amplitudes_are_zero(out)
    type(fusion_alpha_amplitudes_v1), intent(in) :: out

    amplitudes_are_zero = all(out%unsym_real == 0.0_c_double) .and. &
         all(out%unsym_imag == 0.0_c_double) .and. &
         all(out%sym_real == 0.0_c_double) .and. &
         all(out%sym_imag == 0.0_c_double) .and. &
         out%phase_space_J == 0.0_c_double
  end function amplitudes_are_zero

  real(c_double) function amplitude_norm(out)
    type(fusion_alpha_amplitudes_v1), intent(in) :: out

    amplitude_norm = sum(out%sym_real**2 + out%sym_imag**2)
  end function amplitude_norm

  subroutine make_edges(edges, lo_mev, hi_mev)
    real(c_double), intent(out) :: edges(:)
    real(c_double), intent(in) :: lo_mev, hi_mev
    integer :: i

    do i = 1, size(edges)
       edges(i) = (lo_mev + (hi_mev - lo_mev) * &
            real(i - 1, c_double) / real(size(edges) - 1, c_double)) * mev_j
    end do
  end subroutine make_edges

  subroutine check_closure(edges, birth, out, available, label)
    real(c_double), intent(in) :: edges(:), birth(:), available
    type(fusion_alpha_spectrum_v1), intent(in) :: out
    character(len=*), intent(in) :: label
    real(c_double) :: mapped_number, mapped_energy
    real(c_double) :: total_number, total_energy
    integer :: i

    mapped_number = 0.0_c_double
    mapped_energy = 0.0_c_double
    do i = 1, size(birth)
       call check(birth(i) >= 0.0_c_double .and. &
            ieee_is_finite(birth(i)), trim(label) // ' has finite nonnegative bins')
       mapped_number = mapped_number + birth(i)
       mapped_energy = mapped_energy + birth(i) * 0.5_c_double * &
            (edges(i) + edges(i + 1))
    end do
    total_number = mapped_number + out%below_number + out%above_number
    total_energy = mapped_energy + out%below_energy_J + out%above_energy_J
    call check(abs(mapped_number - out%mapped_number) < 4.0e-12_c_double, &
         trim(label) // ' mapped number diagnostic')
    call check(abs(mapped_energy - out%mapped_energy_J) < 3.0e-11_c_double * &
         available, trim(label) // ' mapped energy diagnostic')
    call check(abs(total_number - 3.0_c_double) < 4.0e-12_c_double, &
         trim(label) // ' three-alpha closure including spill')
    call check(abs(total_energy - available) < 3.0e-11_c_double * available, &
         trim(label) // ' CM-energy closure including spill')
    call check(abs(out%number_residual) < 4.0e-12_c_double .and. &
         abs(out%energy_residual_J) < 3.0e-11_c_double * available, &
         trim(label) // ' reported residuals are closed')
  end subroutine check_closure

  subroutine verify_layout()
    type(fusion_alpha_spectrum_v1) :: out
    real(c_double) :: d
    integer(c_int) :: j

    call check(c_sizeof(out) == 11 * c_sizeof(d) + 2 * c_sizeof(j), &
         'alpha spectrum struct has eleven c_double and two c_int fields')
  end subroutine verify_layout

  subroutine verify_normalization_and_spill()
    real(c_double) :: edges(cells + 1), birth(cells), available, cutoff
    real(c_double) :: mapped_number
    type(fusion_alpha_spectrum_v1) :: out
    integer(c_int) :: status

    available = 9.3_c_double * mev_j
    cutoff = 0.001_c_double * mev_j
    call make_edges(edges, 2.0_c_double, 4.0_c_double)
    birth = -1.0_c_double
    out%mapped_number = -1.0_c_double
    call fusion_alpha_spectrum_grid(1_c_int, available, cutoff, &
         0.76_c_double, 0.67_c_double, 20_c_int, 20_c_int, cells, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_OK, 'pure l=1 grid call returns OK')
    call check(out%normalization_J2 > 0.0_c_double .and. &
         ieee_is_finite(out%normalization_J2), &
         'pure l=1 normalization is finite and positive')
    call check(out%quadrature_events == 400_c_int, &
         'quadrature event count is nq*ncos')
    call check(out%pruned_events >= 0_c_int, 'pruned event count is nonnegative')
    call check_closure(edges, birth, out, available, 'pure l=1')
    mapped_number = sum(birth)
    call check(mapped_number > 0.0_c_double .and. &
         out%below_number + out%above_number > 0.0_c_double, &
         'grid has mapped and explicit spill particles')
  end subroutine verify_normalization_and_spill

  subroutine verify_cutoff_wrapper()
    type(fusion_alpha_amplitudes_v1) :: strict, cutoff_out
    real(c_double) :: available, intermediate, cutoff
    integer(c_int) :: status, pruned
    integer :: i

    available = 9.3_c_double * mev_j
    intermediate = 3.129_c_double * mev_j
    cutoff = 0.001_c_double * mev_j

    call fusion_alpha_amplitudes(2_c_int, available, intermediate, &
         0.3_c_double, strict, status)
    call check(status == PB11_STATUS_OK, &
         'strict l=2 amplitude reference returns OK')
    call fusion_alpha_amplitudes_cutoff(2_c_int, available, intermediate, &
         0.3_c_double, cutoff, cutoff_out, pruned, status)
    call check(status == PB11_STATUS_OK .and. pruned == 0_c_int, &
         'cutoff l=2 interior event returns OK and is unpruned')
    do i = 1, 5
       call check(cutoff_out%sym_real(i) == strict%sym_real(i) .and. &
            cutoff_out%sym_imag(i) == strict%sym_imag(i), &
            'cutoff interior amplitudes match strict API')
    end do

    call fusion_alpha_amplitudes_cutoff(2_c_int, available, 0.0_c_double, &
         0.3_c_double, cutoff, cutoff_out, pruned, status)
    call check(status == PB11_STATUS_OK .and. pruned == 1_c_int .and. &
         amplitude_norm(cutoff_out) > 0.0_c_double .and. &
         cutoff_out%phase_space_J == 0.0_c_double, &
         'cutoff q=0 endpoint retains other permutations')
    call fusion_alpha_amplitudes_cutoff(2_c_int, available, available, &
         0.3_c_double, cutoff, cutoff_out, pruned, status)
    call check(status == PB11_STATUS_OK .and. pruned == 1_c_int .and. &
         amplitude_norm(cutoff_out) > 0.0_c_double .and. &
         cutoff_out%phase_space_J == 0.0_c_double, &
         'cutoff q=A endpoint retains other permutations')

    cutoff_out%phase_space_J = -1.0_c_double
    call fusion_alpha_amplitudes_cutoff(2_c_int, available, intermediate, &
         0.3_c_double, 0.0009_c_double * mev_j, cutoff_out, pruned, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. pruned == 0_c_int .and. &
         amplitudes_are_zero(cutoff_out), &
         'out-of-domain cutoff clears amplitude output and count')
  end subroutine verify_cutoff_wrapper

  subroutine verify_mixture_endpoints()
    real(c_double) :: edges(cells + 1), birth_l1(cells), birth_l3(cells)
    real(c_double) :: birth_m1(cells), birth_m3(cells), birth_mix(cells)
    real(c_double) :: available, cutoff
    type(fusion_alpha_spectrum_v1) :: l1, l3, m1, m3, mix
    integer(c_int) :: status
    integer :: i

    available = 9.3_c_double * mev_j
    cutoff = 0.001_c_double * mev_j
    call make_edges(edges, 0.0_c_double, 8.0_c_double)

    call fusion_alpha_spectrum_grid(1_c_int, available, cutoff, &
         0.2_c_double, 0.0_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth_l1, l1, status)
    call check(status == PB11_STATUS_OK, 'pure l=1 endpoint reference returns OK')
    call fusion_alpha_spectrum_grid(3_c_int, available, cutoff, &
         0.2_c_double, 0.0_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth_l3, l3, status)
    call check(status == PB11_STATUS_OK, 'pure l=3 endpoint reference returns OK')

    call fusion_alpha_spectrum_grid(13_c_int, available, cutoff, &
         1.0_c_double, 2.1_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth_m1, m1, status)
    call check(status == PB11_STATUS_OK, 'mode-13 l=1 endpoint returns OK')
    call check(m1%l1_normalization_J2 > 0.0_c_double .and. &
         m1%l3_normalization_J2 > 0.0_c_double, &
         'mode-13 records both unit-normalized bases')
    do i = 1, cells
       call check(abs(birth_m1(i) - birth_l1(i)) < 2.0e-12_c_double, &
            'mode-13 fraction one matches pure l=1 bins')
    end do
    call check(abs(m1%mapped_number - l1%mapped_number) < 2.0e-12_c_double .and. &
         abs(m1%mapped_energy_J - l1%mapped_energy_J) < 3.0e-11_c_double * available, &
         'mode-13 fraction one matches pure l=1 ledger')

    call fusion_alpha_spectrum_grid(13_c_int, available, cutoff, &
         0.0_c_double, -1.4_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth_m3, m3, status)
    call check(status == PB11_STATUS_OK, 'mode-13 l=3 endpoint returns OK')
    do i = 1, cells
       call check(abs(birth_m3(i) - birth_l3(i)) < 2.0e-12_c_double, &
            'mode-13 fraction zero matches pure l=3 bins')
    end do
    call check(abs(m3%mapped_number - l3%mapped_number) < 2.0e-12_c_double .and. &
         abs(m3%mapped_energy_J - l3%mapped_energy_J) < 3.0e-11_c_double * available, &
         'mode-13 fraction zero matches pure l=3 ledger')

    call fusion_alpha_spectrum_grid(13_c_int, available, cutoff, &
         0.5_c_double, 0.3_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth_mix, mix, status)
    call check(status == PB11_STATUS_OK .and. mix%normalization_J2 > 0.0_c_double, &
         'interior coherent mixture returns normalized spectrum')
    call check_closure(edges, birth_mix, mix, available, 'mode-13 interior')
  end subroutine verify_mixture_endpoints

  subroutine verify_errors_and_clearing()
    real(c_double) :: edges(cells + 1), short_edges(cells), birth(cells)
    real(c_double) :: available, cutoff
    type(fusion_alpha_spectrum_v1) :: out
    integer(c_int) :: status

    available = 9.3_c_double * mev_j
    cutoff = 0.001_c_double * mev_j
    call make_edges(edges, 0.0_c_double, 8.0_c_double)
    call make_edges(short_edges, 0.0_c_double, 8.0_c_double)

    birth = 1.0_c_double
    out%mapped_number = -1.0_c_double
    call fusion_alpha_spectrum_grid(99_c_int, available, cutoff, &
         0.5_c_double, 0.0_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(out), 'invalid mode clears all outputs')

    birth = 1.0_c_double
    call fusion_alpha_spectrum_grid(1_c_int, 0.0_c_double, cutoff, &
         0.5_c_double, 0.0_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(out), 'nonpositive available energy clears outputs')

    birth = 1.0_c_double
    call fusion_alpha_spectrum_grid(1_c_int, available, 0.0009_c_double * mev_j, &
         0.5_c_double, 0.0_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(out), 'out-of-domain cutoff clears outputs')

    birth = 1.0_c_double
    call fusion_alpha_spectrum_grid(13_c_int, available, cutoff, &
         -0.01_c_double, 0.0_c_double, 16_c_int, 16_c_int, cells, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(out), 'invalid mixture fraction clears outputs')

    birth = 1.0_c_double
    call fusion_alpha_spectrum_grid(1_c_int, available, cutoff, &
         0.5_c_double, 0.0_c_double, 3_c_int, 16_c_int, cells, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(out), 'too-small quadrature clears outputs')

    birth = 1.0_c_double
    call fusion_alpha_spectrum_grid(1_c_int, available, cutoff, &
         0.5_c_double, 0.0_c_double, 16_c_int, 1025_c_int, cells, edges, &
         birth, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(out), 'too-large quadrature clears outputs')

    birth = 1.0_c_double
    call fusion_alpha_spectrum_grid(1_c_int, available, cutoff, &
         0.5_c_double, 0.0_c_double, 16_c_int, 16_c_int, cells, short_edges, &
         birth, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. all(birth == 0.0_c_double) &
         .and. spectrum_is_zero(out), 'edge extent mismatch clears outputs')
  end subroutine verify_errors_and_clearing

end program test_fusion_alpha_grid_fortran
