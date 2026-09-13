program test_fusion_spectrum_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_spectrum_fortran
  implicit none

  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_nuclear_coulomb_boundary()
  call verify_amplitudes()
  call verify_failure_clears()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran fusion-spectrum test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran fusion-spectrum tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function nuclear_is_zero(out)
    type(fusion_nuclear_coulomb_v1), intent(in) :: out

    nuclear_is_zero = out%log_penetrability == 0.0_c_double .and. &
         out%shift == 0.0_c_double .and. &
         out%phase_real == 0.0_c_double .and. &
         out%phase_imag == 0.0_c_double .and. &
         out%rho == 0.0_c_double
  end function nuclear_is_zero

  logical function amplitudes_are_zero(out)
    type(fusion_alpha_amplitudes_v1), intent(in) :: out

    amplitudes_are_zero = all(out%unsym_real == 0.0_c_double) .and. &
         all(out%unsym_imag == 0.0_c_double) .and. &
         all(out%sym_real == 0.0_c_double) .and. &
         all(out%sym_imag == 0.0_c_double) .and. &
         out%phase_space_J == 0.0_c_double
  end function amplitudes_are_zero

  subroutine verify_layout()
    type(fusion_nuclear_coulomb_v1) :: coulomb
    type(fusion_alpha_amplitudes_v1) :: amplitudes

    call check(c_sizeof(coulomb) == 5 * c_sizeof(0.0_c_double), &
         'nuclear Coulomb struct has five c_double fields')
    call check(c_sizeof(amplitudes) == 21 * c_sizeof(0.0_c_double), &
         'alpha amplitudes struct has twenty-one c_double fields')
  end subroutine verify_layout

  subroutine verify_nuclear_coulomb_boundary()
    type(fusion_nuclear_coulomb_v1) :: out
    real(c_double) :: phase_norm
    integer(c_int) :: status

    call fusion_nuclear_coulomb(FUSION_ALPHA_ALPHA_L2, 3.129_c_double * &
         mev_j, out, status)
    call check(status == PB11_STATUS_OK, &
         '3.129 MeV nuclear Coulomb call returns OK')
    call check(ieee_is_finite(out%log_penetrability) .and. &
         ieee_is_finite(out%shift) .and. &
         ieee_is_finite(out%phase_real) .and. &
         ieee_is_finite(out%phase_imag) .and. &
         ieee_is_finite(out%rho), &
         '3.129 MeV nuclear Coulomb result is finite')
    phase_norm = hypot(out%phase_real, out%phase_imag)
    call check(abs(phase_norm - 1.0_c_double) < 3.0e-6_c_double, &
         'nuclear Coulomb phase has unit norm')
    call check(out%rho > 0.0_c_double, 'nuclear Coulomb rho is positive')
  end subroutine verify_nuclear_coulomb_boundary

  subroutine verify_amplitudes()
    type(fusion_alpha_amplitudes_v1) :: out
    real(c_double) :: unsym_norm, sym_norm
    integer(c_int) :: status

    call fusion_alpha_amplitudes(2_c_int, 9.3_c_double * mev_j, &
         3.129_c_double * mev_j, 0.0_c_double, out, status)
    call check(status == PB11_STATUS_OK, &
         'valid alpha amplitude call returns OK')
    call check(all(ieee_is_finite(out%unsym_real)) .and. &
         all(ieee_is_finite(out%unsym_imag)) .and. &
         all(ieee_is_finite(out%sym_real)) .and. &
         all(ieee_is_finite(out%sym_imag)) .and. &
         ieee_is_finite(out%phase_space_J), &
         'alpha amplitudes are finite')
    unsym_norm = sum(out%unsym_real**2 + out%unsym_imag**2)
    sym_norm = sum(out%sym_real**2 + out%sym_imag**2)
    call check(unsym_norm > 0.0_c_double .and. sym_norm > 0.0_c_double, &
         'alpha amplitude norms are positive')
    call check(out%phase_space_J > 0.0_c_double, &
         'alpha amplitude phase-space factor is positive')
  end subroutine verify_amplitudes

  subroutine verify_failure_clears()
    type(fusion_nuclear_coulomb_v1) :: coulomb
    type(fusion_alpha_amplitudes_v1) :: amplitudes
    integer(c_int) :: status

    coulomb%log_penetrability = -1.0_c_double
    coulomb%shift = -1.0_c_double
    coulomb%phase_real = -1.0_c_double
    coulomb%phase_imag = -1.0_c_double
    coulomb%rho = -1.0_c_double
    call fusion_nuclear_coulomb(FUSION_ALPHA_ALPHA_L2, 0.0009_c_double * &
         mev_j, coulomb, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         nuclear_is_zero(coulomb), &
         'out-of-domain nuclear Coulomb call clears output')

    amplitudes%unsym_real = -1.0_c_double
    amplitudes%unsym_imag = -1.0_c_double
    amplitudes%sym_real = -1.0_c_double
    amplitudes%sym_imag = -1.0_c_double
    amplitudes%phase_space_J = -1.0_c_double
    call fusion_alpha_amplitudes(2_c_int, 9.3_c_double * mev_j, &
         0.0_c_double, 0.0_c_double, amplitudes, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         amplitudes_are_zero(amplitudes), &
         'invalid alpha amplitude call clears output')
  end subroutine verify_failure_clears

end program test_fusion_spectrum_fortran
