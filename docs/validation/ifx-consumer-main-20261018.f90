program consumer
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use pb11_fortran
  implicit none
  real(c_double) :: rate
  integer(c_int) :: status
  call pb11_reactivity_fast(100.0_c_double, rate, status)
  if (status /= PB11_STATUS_OK) stop 2
  if (.not. ieee_is_finite(rate)) stop 3
  if (rate <= 0.0_c_double) stop 4
  write(*,'(a,es24.16)') 'pB reactivity at 100 keV [m3/s] = ', rate
end program consumer
