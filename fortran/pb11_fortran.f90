module pb11_fortran
  !! Fortran ISO_C_BINDING wrappers for the p-11B thermal reaction library.
  !!
  !! All public routines return a status code.  On an error, VALUE is set to
  !! zero and STATUS is non-zero; callers should check STATUS before using it.
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  ! Fortran identifiers are case-insensitive, so method constants cannot have
  ! the same names as the explicit wrapper subroutines below.
  integer(c_int), parameter, public :: PB11_METHOD_INTEGRAL = 0_c_int
  integer(c_int), parameter, public :: PB11_METHOD_FAST = 1_c_int

  public :: pb11_sfactor
  public :: pb11_cross_section
  public :: pb11_reactivity
  public :: pb11_reactivity_integral
  public :: pb11_reactivity_fast

  interface
     function c_pb11_sfactor(energy_mev, value) bind(C, name="pb11_c_sfactor") &
          result(status)
       import :: c_double, c_int
       real(c_double), value :: energy_mev
       real(c_double) :: value
       integer(c_int) :: status
     end function c_pb11_sfactor

     function c_pb11_cross_section(energy_mev, value) &
          bind(C, name="pb11_c_cross_section") result(status)
       import :: c_double, c_int
       real(c_double), value :: energy_mev
       real(c_double) :: value
       integer(c_int) :: status
     end function c_pb11_cross_section

     function c_pb11_reactivity_integral(temperature_kev, value) &
          bind(C, name="pb11_c_reactivity_integral") result(status)
       import :: c_double, c_int
       real(c_double), value :: temperature_kev
       real(c_double) :: value
       integer(c_int) :: status
     end function c_pb11_reactivity_integral

     function c_pb11_reactivity_fast(temperature_kev, value) &
          bind(C, name="pb11_c_reactivity_fast") result(status)
       import :: c_double, c_int
       real(c_double), value :: temperature_kev
       real(c_double) :: value
       integer(c_int) :: status
     end function c_pb11_reactivity_fast

     function c_pb11_reactivity(temperature_kev, method, value) &
          bind(C, name="pb11_c_reactivity") result(status)
       import :: c_double, c_int
       real(c_double), value :: temperature_kev
       integer(c_int), value :: method
       real(c_double) :: value
       integer(c_int) :: status
     end function c_pb11_reactivity
  end interface

contains

  subroutine pb11_sfactor(energy_mev, value, status)
    real(c_double), intent(in) :: energy_mev
    real(c_double), intent(out) :: value
    integer(c_int), intent(out) :: status

    value = 0.0_c_double
    status = c_pb11_sfactor(energy_mev, value)
  end subroutine pb11_sfactor

  subroutine pb11_cross_section(energy_mev, value, status)
    real(c_double), intent(in) :: energy_mev
    real(c_double), intent(out) :: value
    integer(c_int), intent(out) :: status

    value = 0.0_c_double
    status = c_pb11_cross_section(energy_mev, value)
  end subroutine pb11_cross_section

  subroutine pb11_reactivity(temperature_kev, method, value, status)
    real(c_double), intent(in) :: temperature_kev
    integer(c_int), intent(in) :: method
    real(c_double), intent(out) :: value
    integer(c_int), intent(out) :: status

    value = 0.0_c_double
    status = c_pb11_reactivity(temperature_kev, method, value)
  end subroutine pb11_reactivity

  subroutine pb11_reactivity_integral(temperature_kev, value, status)
    real(c_double), intent(in) :: temperature_kev
    real(c_double), intent(out) :: value
    integer(c_int), intent(out) :: status

    value = 0.0_c_double
    status = c_pb11_reactivity_integral(temperature_kev, value)
  end subroutine pb11_reactivity_integral

  subroutine pb11_reactivity_fast(temperature_kev, value, status)
    real(c_double), intent(in) :: temperature_kev
    real(c_double), intent(out) :: value
    integer(c_int), intent(out) :: status

    value = 0.0_c_double
    status = c_pb11_reactivity_fast(temperature_kev, value)
  end subroutine pb11_reactivity_fast

end module pb11_fortran
