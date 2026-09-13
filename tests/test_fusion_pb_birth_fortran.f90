program test_fusion_pb_birth_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_pb_birth_fortran
  implicit none

  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  real(c_double), parameter :: available = 8.68_c_double * mev_j
  real(c_double), parameter :: q_gs = 0.09184_c_double * mev_j
  integer(c_int), parameter :: cells = 64_c_int
  integer(c_int), parameter :: nq = 16_c_int, ncos = 16_c_int
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_alpha0()
  call verify_source_mixture()
  call verify_errors_and_clearing()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran pB birth test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran pB birth tests passed'

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
       return
    end if
    scale = max(abs(actual), abs(expected))
    if (scale == 0.0_c_double) then
       close_value = .true.
    else
       close_value = abs(actual - expected) <= tolerance * scale
    end if
  end function close_value

  logical function birth_is_zero(out)
    type(fusion_pb_birth_v1), intent(in) :: out

    birth_is_zero = all([out%mapped_number, out%mapped_energy_J, &
         out%below_number, out%below_energy_J, out%above_number, &
         out%above_energy_J, out%number_residual, out%energy_residual_J, &
         out%alpha0_fraction, out%low_alpha1_fraction, &
         out%broad_alpha1_fraction, out%primary_alpha0_energy_J, &
         out%secondary_min_J, out%secondary_max_J] == 0.0_c_double)
  end function birth_is_zero

  logical function birth_is_finite(out)
    type(fusion_pb_birth_v1), intent(in) :: out

    birth_is_finite = all(ieee_is_finite([out%mapped_number, &
         out%mapped_energy_J, out%below_number, out%below_energy_J, &
         out%above_number, out%above_energy_J, out%number_residual, &
         out%energy_residual_J, out%alpha0_fraction, &
         out%low_alpha1_fraction, out%broad_alpha1_fraction, &
         out%primary_alpha0_energy_J, out%secondary_min_J, &
         out%secondary_max_J]))
  end function birth_is_finite

  subroutine make_edges(edges)
    real(c_double), intent(out), target :: edges(:)
    integer :: i

    do i = 1, size(edges)
       edges(i) = 12.0_c_double * mev_j * real(i - 1, c_double) / &
            real(size(edges) - 1, c_double)
    end do
  end subroutine make_edges

  subroutine check_closure(birth, out, label)
    real(c_double), intent(in) :: birth(:)
    type(fusion_pb_birth_v1), intent(in) :: out
    character(len=*), intent(in) :: label
    real(c_double) :: number_total, energy_total

    number_total = out%mapped_number + out%below_number + out%above_number
    energy_total = out%mapped_energy_J + out%below_energy_J + &
         out%above_energy_J
    call check(close_value(number_total, 3.0_c_double, 2.0e-11_c_double), &
         label // ' closes three-particle number')
    call check(close_value(energy_total, available, 2.0e-11_c_double), &
         label // ' closes available CM energy')
    call check(all(birth >= 0.0_c_double) .and. all(ieee_is_finite(birth)), &
         label // ' has finite nonnegative cell births')
  end subroutine check_closure

  subroutine verify_layout()
    type(fusion_pb_birth_v1) :: out
    real(c_double) :: d

    call check(c_sizeof(out) == 14_c_size_t * c_sizeof(d), &
         'pB birth result layout is fourteen c_doubles')
  end subroutine verify_layout

  subroutine verify_alpha0()
    real(c_double), target :: edges(cells + 1), birth(cells)
    type(fusion_pb_birth_v1), target :: out
    integer(c_int) :: status

    call make_edges(edges)
    birth = -1.0_c_double
    call fusion_pb_alpha0_grid(available, q_gs, cells, edges, birth, out, &
         status)
    call check(status == PB11_STATUS_OK, 'alpha0 grid returns OK')
    call check(birth_is_finite(out), 'alpha0 result is finite')
    call check(out%alpha0_fraction == 1.0_c_double .and. &
         out%low_alpha1_fraction == 0.0_c_double .and. &
         out%broad_alpha1_fraction == 0.0_c_double, &
         'alpha0 grid reports its unit source fraction')
    call check(out%primary_alpha0_energy_J >= 0.0_c_double .and. &
         out%primary_alpha0_energy_J <= available .and. &
         out%secondary_min_J >= 0.0_c_double .and. &
         out%secondary_max_J >= out%secondary_min_J, &
         'alpha0 endpoint diagnostics are ordered and in range')
    call check_closure(birth, out, 'alpha0 grid')
  end subroutine verify_alpha0

  subroutine verify_source_mixture()
    real(c_double), target :: edges(cells + 1), birth(cells)
    type(fusion_pb_birth_v1), target :: out
    real(c_double), parameter :: f0 = 0.05_c_double
    real(c_double), parameter :: flow = 0.65_c_double
    integer(c_int) :: status

    call make_edges(edges)
    birth = -1.0_c_double
    call fusion_pb_cm_source_grid(available, q_gs, f0, flow, 13_c_int, &
         0_c_int, 0.001_c_double * mev_j, 0.76_c_double, 0.2_c_double, &
         nq, ncos, cells, edges, birth, out, status)
    call check(status == PB11_STATUS_OK, 'CM source mixture returns OK')
    call check(birth_is_finite(out), 'CM source result is finite')
    call check(close_value(out%alpha0_fraction, f0, 2.0e-14_c_double) .and. &
         close_value(out%low_alpha1_fraction, flow, 2.0e-14_c_double) .and. &
         close_value(out%broad_alpha1_fraction, 1.0_c_double - f0 - flow, &
         2.0e-14_c_double), 'CM source reports the supplied event fractions')
    call check(out%primary_alpha0_energy_J >= 0.0_c_double .and. &
         out%secondary_min_J >= 0.0_c_double .and. &
         out%secondary_max_J >= out%secondary_min_J, &
         'CM source always supplies alpha0 endpoint diagnostics')
    call check_closure(birth, out, 'CM source mixture')
  end subroutine verify_source_mixture

  subroutine verify_errors_and_clearing()
    real(c_double), target :: edges(cells + 1), short_edges(cells)
    real(c_double), target :: birth(cells), short_birth(cells - 1)
    type(fusion_pb_birth_v1), target :: out
    integer(c_int) :: status

    call make_edges(edges)
    call make_edges(short_edges)

    birth = 1.0_c_double
    out%mapped_number = -1.0_c_double
    call fusion_pb_alpha0_grid(available, q_gs, cells, short_edges, birth, &
         out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(birth == 0.0_c_double) .and. birth_is_zero(out), &
         'alpha0 extent mismatch is rejected before C_LOC and cleared')

    short_birth = 1.0_c_double
    call fusion_pb_alpha0_grid(available, q_gs, cells, edges, short_birth, &
         out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(short_birth == 0.0_c_double) .and. birth_is_zero(out), &
         'alpha0 birth extent mismatch is cleared')

    birth = 1.0_c_double
    call fusion_pb_alpha0_grid(0.0_c_double, q_gs, cells, edges, birth, out, &
         status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         all(birth == 0.0_c_double) .and. birth_is_zero(out), &
         'invalid alpha0 energy clears C outputs')

    birth = 1.0_c_double
    call fusion_pb_cm_source_grid(available, q_gs, -0.1_c_double, 0.5_c_double, &
         13_c_int, 0_c_int, 0.001_c_double * mev_j, 0.76_c_double, &
         0.0_c_double, nq, ncos, cells, edges, birth, out, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         all(birth == 0.0_c_double) .and. birth_is_zero(out), &
         'invalid source fraction clears C outputs')

    birth = 1.0_c_double
    call fusion_pb_cm_source_grid(available, q_gs, 0.0_c_double, 0.0_c_double, &
         99_c_int, 0_c_int, 0.001_c_double * mev_j, 0.76_c_double, &
         0.0_c_double, nq, ncos, cells, edges, birth, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(birth == 0.0_c_double) .and. birth_is_zero(out), &
         'invalid broad mode clears C outputs')
  end subroutine verify_errors_and_clearing

end program test_fusion_pb_birth_fortran
