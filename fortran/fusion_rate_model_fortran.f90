module fusion_rate_model_fortran
  !! ISO_C_BINDING wrappers for the bounded continuation rate-model ABI.
  !!
  !! The segment fields reuse fusion_beam_fortran's C-compatible window type;
  !! this module does not introduce a second definition of that layout.
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

  ! Exact C layout: four consecutive fusion_beam_window_v1 values.
  type, bind(C), public :: fusion_rate_model_v1
     type(fusion_beam_window_v1) :: total
     type(fusion_beam_window_v1) :: fit
     type(fusion_beam_window_v1) :: below
     type(fusion_beam_window_v1) :: above
  end type fusion_rate_model_v1

  public :: fusion_cross_section_model
  public :: fusion_beam_maxwellian_model
  public :: fusion_thermal_pair_maxwellian_model

  interface
     function c_fusion_cross_section_model(channel, continuation, pb_low, &
          relative_energy_J, cross_section_m2) bind(C, &
          name="fusion_c_cross_section_model") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       integer(c_int), value :: continuation
       integer(c_int), value :: pb_low
       real(c_double), value :: relative_energy_J
       type(c_ptr), value :: cross_section_m2
       integer(c_int) :: status
     end function c_fusion_cross_section_model

     function c_fusion_beam_maxwellian_model(channel, continuation, pb_low, &
          projectile_mass_kg, target_mass_kg, projectile_energy_J, &
          target_kT_J, out) bind(C, &
          name="fusion_c_beam_maxwellian_model") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       integer(c_int), value :: continuation
       integer(c_int), value :: pb_low
       real(c_double), value :: projectile_mass_kg
       real(c_double), value :: target_mass_kg
       real(c_double), value :: projectile_energy_J
       real(c_double), value :: target_kT_J
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_maxwellian_model

     function c_fusion_thermal_pair_maxwellian_model(channel, continuation, &
          pb_low, reactant_a_mass_kg, reactant_b_mass_kg, reactant_a_kT_J, &
          reactant_b_kT_J, out) bind(C, &
          name="fusion_c_thermal_pair_maxwellian_model") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       integer(c_int), value :: continuation
       integer(c_int), value :: pb_low
       real(c_double), value :: reactant_a_mass_kg
       real(c_double), value :: reactant_b_mass_kg
       real(c_double), value :: reactant_a_kT_J
       real(c_double), value :: reactant_b_kT_J
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_thermal_pair_maxwellian_model
  end interface

contains

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

  subroutine clear_rate_model(out)
    type(fusion_rate_model_v1), intent(out) :: out

    call clear_window(out%total)
    call clear_window(out%fit)
    call clear_window(out%below)
    call clear_window(out%above)
  end subroutine clear_rate_model

  subroutine fusion_cross_section_model(channel, continuation, pb_low, &
       relative_energy_J, cross_section_m2, status)
    integer(c_int), intent(in) :: channel, continuation, pb_low
    real(c_double), intent(in) :: relative_energy_J
    real(c_double), intent(out), target :: cross_section_m2
    integer(c_int), intent(out) :: status

    cross_section_m2 = 0.0_c_double
    status = c_fusion_cross_section_model(channel, continuation, pb_low, &
         relative_energy_J, c_loc(cross_section_m2))
    if (status /= PB11_STATUS_OK) cross_section_m2 = 0.0_c_double
  end subroutine fusion_cross_section_model

  subroutine fusion_beam_maxwellian_model(channel, continuation, pb_low, &
       projectile_mass_kg, target_mass_kg, projectile_energy_J, target_kT_J, &
       out, status)
    integer(c_int), intent(in) :: channel, continuation, pb_low
    real(c_double), intent(in) :: projectile_mass_kg, target_mass_kg
    real(c_double), intent(in) :: projectile_energy_J, target_kT_J
    type(fusion_rate_model_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_rate_model(out)
    status = c_fusion_beam_maxwellian_model(channel, continuation, pb_low, &
         projectile_mass_kg, target_mass_kg, projectile_energy_J, target_kT_J, &
         c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_rate_model(out)
  end subroutine fusion_beam_maxwellian_model

  subroutine fusion_thermal_pair_maxwellian_model(channel, continuation, &
       pb_low, reactant_a_mass_kg, reactant_b_mass_kg, reactant_a_kT_J, &
       reactant_b_kT_J, out, status)
    integer(c_int), intent(in) :: channel, continuation, pb_low
    real(c_double), intent(in) :: reactant_a_mass_kg, reactant_b_mass_kg
    real(c_double), intent(in) :: reactant_a_kT_J, reactant_b_kT_J
    type(fusion_rate_model_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_rate_model(out)
    status = c_fusion_thermal_pair_maxwellian_model(channel, continuation, &
         pb_low, reactant_a_mass_kg, reactant_b_mass_kg, reactant_a_kT_J, &
         reactant_b_kT_J, c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_rate_model(out)
  end subroutine fusion_thermal_pair_maxwellian_model

end module fusion_rate_model_fortran
