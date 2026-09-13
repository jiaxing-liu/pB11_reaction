module fusion_beam_fortran
  !! ISO_C_BINDING wrappers for the scalar beam-window ABI and its
  !! cross-section domain query.  The beam result keeps every C field,
  !! including domain_incomplete, visible to the caller.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = 0_c_int
  integer(c_int), parameter, public :: FUSION_DD_TP = 1_c_int
  integer(c_int), parameter, public :: FUSION_DD_HE3N = 2_c_int
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = 3_c_int
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = 4_c_int
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = 5_c_int

  ! Exact C layout: ten c_double fields followed by one C int.
  type, bind(C), public :: fusion_beam_window_v1
     real(c_double) :: resolved_reactivity_m3_s
     real(c_double) :: projectile_energy_reactivity_J_m3_s
     real(c_double) :: target_energy_reactivity_J_m3_s
     real(c_double) :: relative_energy_reactivity_J_m3_s
     real(c_double) :: cm_energy_reactivity_J_m3_s
     real(c_double) :: resolved_pair_probability
     real(c_double) :: unresolved_pair_probability
     real(c_double) :: unresolved_relative_speed_m_s
     real(c_double) :: quadrature_error_m3_s
     real(c_double) :: energy_identity_error_J_m3_s
     integer(c_int) :: domain_incomplete
  end type fusion_beam_window_v1

  public :: fusion_beam_maxwellian_window
  public :: fusion_thermal_pair_maxwellian_window
  public :: fusion_cross_section_domain

  interface
     function c_fusion_beam_maxwellian_window(channel, projectile_mass_kg, &
          target_mass_kg, projectile_energy_J, target_kT_J, out) &
          bind(C, name="fusion_c_beam_maxwellian_window") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       real(c_double), value :: projectile_mass_kg
       real(c_double), value :: target_mass_kg
       real(c_double), value :: projectile_energy_J
       real(c_double), value :: target_kT_J
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_maxwellian_window

     function c_fusion_thermal_pair_maxwellian_window(channel, &
          reactant_a_mass_kg, reactant_b_mass_kg, reactant_a_kT_J, &
          reactant_b_kT_J, out) bind(C, &
          name="fusion_c_thermal_pair_maxwellian_window") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       real(c_double), value :: reactant_a_mass_kg
       real(c_double), value :: reactant_b_mass_kg
       real(c_double), value :: reactant_a_kT_J
       real(c_double), value :: reactant_b_kT_J
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_thermal_pair_maxwellian_window

     function c_fusion_cross_section_domain(channel, minimum_J, maximum_J) &
          bind(C, name="fusion_c_cross_section_domain") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       type(c_ptr), value :: minimum_J
       type(c_ptr), value :: maximum_J
       integer(c_int) :: status
     end function c_fusion_cross_section_domain
  end interface

contains

  subroutine clear_beam_window(out)
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
  end subroutine clear_beam_window

  subroutine fusion_beam_maxwellian_window(channel, projectile_mass_kg, &
       target_mass_kg, projectile_energy_J, target_kT_J, out, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: projectile_mass_kg, target_mass_kg
    real(c_double), intent(in) :: projectile_energy_J, target_kT_J
    type(fusion_beam_window_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_beam_window(out)
    status = c_fusion_beam_maxwellian_window(channel, projectile_mass_kg, &
         target_mass_kg, projectile_energy_J, target_kT_J, c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_beam_window(out)
  end subroutine fusion_beam_maxwellian_window

  subroutine fusion_thermal_pair_maxwellian_window(channel, &
       reactant_a_mass_kg, reactant_b_mass_kg, reactant_a_kT_J, &
       reactant_b_kT_J, out, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: reactant_a_mass_kg, reactant_b_mass_kg
    real(c_double), intent(in) :: reactant_a_kT_J, reactant_b_kT_J
    type(fusion_beam_window_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_beam_window(out)
    status = c_fusion_thermal_pair_maxwellian_window(channel, &
         reactant_a_mass_kg, reactant_b_mass_kg, reactant_a_kT_J, &
         reactant_b_kT_J, c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_beam_window(out)
  end subroutine fusion_thermal_pair_maxwellian_window

  subroutine fusion_cross_section_domain(channel, minimum_relative_energy_J, &
       maximum_relative_energy_J, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(out), target :: minimum_relative_energy_J
    real(c_double), intent(out), target :: maximum_relative_energy_J
    integer(c_int), intent(out) :: status

    minimum_relative_energy_J = 0.0_c_double
    maximum_relative_energy_J = 0.0_c_double
    status = c_fusion_cross_section_domain(channel, &
         c_loc(minimum_relative_energy_J), c_loc(maximum_relative_energy_J))
    if (status /= PB11_STATUS_OK) then
       minimum_relative_energy_J = 0.0_c_double
       maximum_relative_energy_J = 0.0_c_double
    end if
  end subroutine fusion_cross_section_domain

end module fusion_beam_fortran
