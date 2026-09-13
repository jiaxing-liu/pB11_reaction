program test_pb11_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use pb11_fortran
  implicit none

  type(pb11_instant_source_v1) :: source
  integer :: failures
  integer(c_int) :: status
  real(c_double) :: value, value_integral, value_fast

  failures = 0

  call pb11_reactivity(100.0_c_double, PB11_METHOD_INTEGRAL, value, status)
  call check(status == PB11_STATUS_OK, &
             "Fortran generic integral call returns OK")
  call check(value > 0.0_c_double .and. ieee_is_finite(value), &
             "Fortran generic integral result is finite and positive")
  value_integral = value

  call pb11_reactivity(100.0_c_double, PB11_METHOD_FAST, value, status)
  call check(status == PB11_STATUS_OK, &
             "Fortran generic fast call returns OK")
  call check(value > 0.0_c_double .and. ieee_is_finite(value), &
             "Fortran generic fast result is finite and positive")
  value_fast = value

  call pb11_reactivity_integral(100.0_c_double, value, status)
  call check(status == PB11_STATUS_OK .and. value == value_integral, &
             "Fortran explicit integral wrapper agrees with generic call")

  call pb11_reactivity_fast(100.0_c_double, value, status)
  call check(status == PB11_STATUS_OK .and. value == value_fast, &
             "Fortran explicit fast wrapper agrees with generic call")

  call pb11_sfactor(0.148_c_double, value, status)
  call check(status == PB11_STATUS_OK .and. value > 3000.0_c_double, &
             "Fortran S-factor wrapper returns the resonance")

  call pb11_reactivity(9.999_c_double, PB11_METHOD_FAST, value, status)
  call check(status == PB11_STATUS_OUT_OF_RANGE .and. value == 0.0_c_double, &
             "Fortran fast wrapper reports an out-of-range temperature")
  call check(ieee_is_finite(value), &
             "Fortran error result remains finite")

  call pb11_reactivity(100.0_c_double, 99_c_int, value, status)
  call check(status == PB11_STATUS_UNKNOWN_METHOD .and. value == 0.0_c_double, &
             "Fortran generic wrapper reports an unknown method")

  call pb11_instant_thermal_source(100._c_double*1.602176634e-16_c_double, &
       1.e20_c_double,2.e19_c_double,0.3_c_double,PB11_METHOD_INTEGRAL,source,status)
  call check(status==PB11_STATUS_OK,'SI instant-source status')
  call check(abs(source%reaction_rate_m3_s/(1.e20_c_double*2.e19_c_double*value_integral)-1)<1.e-12_c_double, &
       'source uses SI density and thermal energy units')
  call check(source%proton_source_m3_s==source%boron_source_m3_s,'equal fuel consumption')
  call check(abs(source%helium4_source_m3_s/(-3*source%proton_source_m3_s)-1)<1.e-14_c_double, &
       'C record layout and helium stoichiometry')
  call check(abs((source%electron_power_W_m3+source%ion_power_W_m3)/source%fusion_power_W_m3-1)<1.e-14_c_double, &
       'C record layout and deposited power closure')
  call pb11_instant_thermal_source(-1._c_double,1.e20_c_double,2.e19_c_double, &
       0.3_c_double,PB11_METHOD_INTEGRAL,source,status)
  call check(status/=PB11_STATUS_OK.and.source%reaction_rate_m3_s==0.and.source%fusion_power_W_m3==0, &
       'failed source resets outputs')

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All p-11B Fortran bindings tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

end program test_pb11_fortran
