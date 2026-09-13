module fusion_thermal_burn_fortran
  !! ISO_C_BINDING wrapper for the fixed six-species/five-channel thermal
  !! reaction-only backward-Euler trial.
  !!
  !! The C ABI owns the physical solve.  This layer only provides the
  !! interoperable ledger, checks every fixed extent before C_LOC, and clears
  !! all result storage when a call is rejected.
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

  integer(c_int), parameter, public :: FUSION_PROTON = 0_c_int
  integer(c_int), parameter, public :: FUSION_DEUTERON = 1_c_int
  integer(c_int), parameter, public :: FUSION_TRITON = 2_c_int
  integer(c_int), parameter, public :: FUSION_HELIUM3 = 3_c_int
  integer(c_int), parameter, public :: FUSION_HELIUM4 = 4_c_int
  integer(c_int), parameter, public :: FUSION_BORON11 = 5_c_int
  integer(c_int), parameter, public :: FUSION_SPECIES_COUNT = 6_c_int

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = 0_c_int
  integer(c_int), parameter, public :: FUSION_DD_TP = 1_c_int
  integer(c_int), parameter, public :: FUSION_DD_HE3N = 2_c_int
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = 3_c_int
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = 4_c_int
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = 5_c_int

  ! Exact C layout: 5 + 6 + 6 + 6 + 1 + 6 + 6 c_double values.
  type, bind(C), public :: fusion_thermal_burn_v1
     real(c_double) :: events_m3(5)
     real(c_double) :: reactant_removed_m3(6)
     real(c_double) :: reactant_removed_energy_J_m3(6)
     real(c_double) :: fast_product_birth_m3(6)
     real(c_double) :: neutron_birth_m3
     real(c_double) :: number_residual_m3(6)
     real(c_double) :: energy_residual_J_m3(6)
  end type fusion_thermal_burn_v1

  public :: fusion_thermal_burn_trial

  interface
     function c_fusion_thermal_burn_trial(dt_s, old_number, old_energy, &
          reactivity, mean_a, mean_b, trial_number, trial_energy, out) &
          bind(C, name="fusion_c_thermal_burn_trial") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: old_number
       type(c_ptr), value :: old_energy
       type(c_ptr), value :: reactivity
       type(c_ptr), value :: mean_a
       type(c_ptr), value :: mean_b
       type(c_ptr), value :: trial_number
       type(c_ptr), value :: trial_energy
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_thermal_burn_trial
  end interface

contains

  subroutine clear_thermal_burn(out)
    type(fusion_thermal_burn_v1), intent(out) :: out

    out%events_m3 = 0.0_c_double
    out%reactant_removed_m3 = 0.0_c_double
    out%reactant_removed_energy_J_m3 = 0.0_c_double
    out%fast_product_birth_m3 = 0.0_c_double
    out%neutron_birth_m3 = 0.0_c_double
    out%number_residual_m3 = 0.0_c_double
    out%energy_residual_J_m3 = 0.0_c_double
  end subroutine clear_thermal_burn

  subroutine fusion_thermal_burn_trial(dt_s, old_number_m3, &
       old_energy_J_m3, reactivity_m3_s, mean_a_energy_J, mean_b_energy_J, &
       trial_number_m3, trial_energy_J_m3, out, status)
    real(c_double), intent(in) :: dt_s
    real(c_double), intent(in), target, contiguous :: old_number_m3(:)
    real(c_double), intent(in), target, contiguous :: old_energy_J_m3(:)
    real(c_double), intent(in), target, contiguous :: reactivity_m3_s(:)
    real(c_double), intent(in), target, contiguous :: mean_a_energy_J(:)
    real(c_double), intent(in), target, contiguous :: mean_b_energy_J(:)
    real(c_double), intent(out), target, contiguous :: trial_number_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_energy_J_m3(:)
    type(fusion_thermal_burn_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_thermal_burn(out)
    trial_number_m3 = 0.0_c_double
    trial_energy_J_m3 = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! The C ABI has fixed-size vectors.  Validate all extents before taking
    ! C_LOC; this also makes malformed calls safe for zero-length actuals.
    if (size(old_number_m3, kind=c_int) /= FUSION_SPECIES_COUNT) return
    if (size(old_energy_J_m3, kind=c_int) /= FUSION_SPECIES_COUNT) return
    if (size(reactivity_m3_s, kind=c_int) /= FUSION_CHANNEL_COUNT) return
    if (size(mean_a_energy_J, kind=c_int) /= FUSION_CHANNEL_COUNT) return
    if (size(mean_b_energy_J, kind=c_int) /= FUSION_CHANNEL_COUNT) return
    if (size(trial_number_m3, kind=c_int) /= FUSION_SPECIES_COUNT) return
    if (size(trial_energy_J_m3, kind=c_int) /= FUSION_SPECIES_COUNT) return

    status = c_fusion_thermal_burn_trial(dt_s, c_loc(old_number_m3(1)), &
         c_loc(old_energy_J_m3(1)), c_loc(reactivity_m3_s(1)), &
         c_loc(mean_a_energy_J(1)), c_loc(mean_b_energy_J(1)), &
         c_loc(trial_number_m3(1)), c_loc(trial_energy_J_m3(1)), c_loc(out))
    if (status /= PB11_STATUS_OK) then
       trial_number_m3 = 0.0_c_double
       trial_energy_J_m3 = 0.0_c_double
       call clear_thermal_burn(out)
    end if
  end subroutine fusion_thermal_burn_trial

end module fusion_thermal_burn_fortran
