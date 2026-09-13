module fusion_pb_population_fortran
  !! ISO_C_BINDING wrappers for the scalar p-11B effective population API.
  !!
  !! The C implementation owns the effective population model and the
  !! source-weighted rate integrations.  This module only mirrors the frozen
  !! C layouts, passes scalar arguments with their explicit C kinds, and
  !! clears all result fields whenever a C call is rejected.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  use fusion_beam_fortran, only : fusion_beam_window_v1
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: FUSION_ENDPOINT_S = 1_c_int
  integer(c_int), parameter, public :: FUSION_HIGH_FLAT = 2_c_int

  integer(c_int), parameter, public :: FUSION_PB_LOW_TB = 0_c_int
  integer(c_int), parameter, public :: FUSION_PB_LOW_NS = 1_c_int
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_MINUS12 = 2_c_int
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_PLUS12 = 3_c_int

  ! Exact C layout: nine consecutive c_double values.
  type, bind(C), public :: fusion_pb_population_v1
     real(c_double) :: total_cross_section_m2
     real(c_double) :: narrow_fit_cross_section_m2
     real(c_double) :: narrow_fit_fraction
     real(c_double) :: alpha0_peak_fraction
     real(c_double) :: narrow_remainder_fraction
     real(c_double) :: other_remainder_fraction
     real(c_double) :: continuum_peak_fraction
     real(c_double) :: continuum_extrapolated_fraction
     real(c_double) :: equivalent_proton_lab_energy_J
  end type fusion_pb_population_v1

  ! Exact C layout: five consecutive fusion_beam_window_v1 values.  Reuse
  ! the type from fusion_beam_fortran so that all moment fields and the
  ! trailing C int have one canonical Fortran definition.
  type, bind(C), public :: fusion_pb_population_rates_v1
     type(fusion_beam_window_v1) :: total
     type(fusion_beam_window_v1) :: alpha0_peak
     type(fusion_beam_window_v1) :: narrow_remainder
     type(fusion_beam_window_v1) :: other_remainder
     type(fusion_beam_window_v1) :: continuum_extrapolated
  end type fusion_pb_population_rates_v1

  public :: fusion_pb_population
  public :: fusion_pb_thermal_population_rates
  public :: fusion_pb_beam_population_rates

  interface
     function c_fusion_pb_population(continuation, pb_low, E_J, narrow_peak, &
          continuum_scale, out) bind(C, name="fusion_c_pb_population") &
          result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: continuation
       integer(c_int), value :: pb_low
       real(c_double), value :: E_J
       real(c_double), value :: narrow_peak
       real(c_double), value :: continuum_scale
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_pb_population

     function c_fusion_pb_thermal_population_rates(continuation, pb_low, &
          mass_a_kg, mass_b_kg, kTa_J, kTb_J, narrow_peak, continuum_scale, &
          out) bind(C, name="fusion_c_pb_thermal_population_rates") &
          result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: continuation
       integer(c_int), value :: pb_low
       real(c_double), value :: mass_a_kg
       real(c_double), value :: mass_b_kg
       real(c_double), value :: kTa_J
       real(c_double), value :: kTb_J
       real(c_double), value :: narrow_peak
       real(c_double), value :: continuum_scale
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_pb_thermal_population_rates

     function c_fusion_pb_beam_population_rates(continuation, pb_low, &
          mass_a_kg, mass_b_kg, projectile_energy_J, target_kT_J, &
          narrow_peak, continuum_scale, out) bind(C, &
          name="fusion_c_pb_beam_population_rates") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: continuation
       integer(c_int), value :: pb_low
       real(c_double), value :: mass_a_kg
       real(c_double), value :: mass_b_kg
       real(c_double), value :: projectile_energy_J
       real(c_double), value :: target_kT_J
       real(c_double), value :: narrow_peak
       real(c_double), value :: continuum_scale
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_pb_beam_population_rates
  end interface

contains

  subroutine clear_population(out)
    type(fusion_pb_population_v1), intent(out) :: out

    out%total_cross_section_m2 = 0.0_c_double
    out%narrow_fit_cross_section_m2 = 0.0_c_double
    out%narrow_fit_fraction = 0.0_c_double
    out%alpha0_peak_fraction = 0.0_c_double
    out%narrow_remainder_fraction = 0.0_c_double
    out%other_remainder_fraction = 0.0_c_double
    out%continuum_peak_fraction = 0.0_c_double
    out%continuum_extrapolated_fraction = 0.0_c_double
    out%equivalent_proton_lab_energy_J = 0.0_c_double
  end subroutine clear_population

  subroutine clear_window(out)
    type(fusion_beam_window_v1), intent(out) :: out

    out%resolved_reactivity_m3_s = 0.0_c_double
    out%projectile_energy_reactivity_J_m3_s = 0.0_c_double
    out%target_energy_reactivity_J_m3_s = 0.0_c_double
    out%relative_energy_reactivity_J_m3_s = 0.0_c_double
    out%cm_energy_reactivity_J_m3_s = 0.0_c_double
    out%resolved_pair_probability = 0.0_c_double
    out%unresolved_pair_probability = 0.0_c_double
    out%unresolved_relative_speed_m_s = 0.0_c_double
    out%quadrature_error_m3_s = 0.0_c_double
    out%energy_identity_error_J_m3_s = 0.0_c_double
    out%domain_incomplete = 0_c_int
  end subroutine clear_window

  subroutine clear_population_rates(out)
    type(fusion_pb_population_rates_v1), intent(out) :: out

    call clear_window(out%total)
    call clear_window(out%alpha0_peak)
    call clear_window(out%narrow_remainder)
    call clear_window(out%other_remainder)
    call clear_window(out%continuum_extrapolated)
  end subroutine clear_population_rates

  subroutine fusion_pb_population(continuation, pb_low, E_J, narrow_peak, &
       continuum_scale, out, status)
    integer(c_int), intent(in) :: continuation, pb_low
    real(c_double), intent(in) :: E_J, narrow_peak, continuum_scale
    type(fusion_pb_population_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_population(out)
    status = c_fusion_pb_population(continuation, pb_low, E_J, narrow_peak, &
         continuum_scale, c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_population(out)
  end subroutine fusion_pb_population

  subroutine fusion_pb_thermal_population_rates(continuation, pb_low, &
       mass_a_kg, mass_b_kg, kTa_J, kTb_J, narrow_peak, continuum_scale, &
       out, status)
    integer(c_int), intent(in) :: continuation, pb_low
    real(c_double), intent(in) :: mass_a_kg, mass_b_kg, kTa_J, kTb_J
    real(c_double), intent(in) :: narrow_peak, continuum_scale
    type(fusion_pb_population_rates_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_population_rates(out)
    status = c_fusion_pb_thermal_population_rates(continuation, pb_low, &
         mass_a_kg, mass_b_kg, kTa_J, kTb_J, narrow_peak, continuum_scale, &
         c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_population_rates(out)
  end subroutine fusion_pb_thermal_population_rates

  subroutine fusion_pb_beam_population_rates(continuation, pb_low, &
       mass_a_kg, mass_b_kg, projectile_energy_J, target_kT_J, narrow_peak, &
       continuum_scale, out, status)
    integer(c_int), intent(in) :: continuation, pb_low
    real(c_double), intent(in) :: mass_a_kg, mass_b_kg
    real(c_double), intent(in) :: projectile_energy_J, target_kT_J
    real(c_double), intent(in) :: narrow_peak, continuum_scale
    type(fusion_pb_population_rates_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_population_rates(out)
    status = c_fusion_pb_beam_population_rates(continuation, pb_low, &
         mass_a_kg, mass_b_kg, projectile_energy_J, target_kT_J, narrow_peak, &
         continuum_scale, c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_population_rates(out)
  end subroutine fusion_pb_beam_population_rates

end module fusion_pb_population_fortran
