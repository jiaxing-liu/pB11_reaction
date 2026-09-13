module fusion_target_burn_fortran
  !! ISO_C_BINDING wrappers for the distinct fast/thermal target burn trial.
  !! Window diagnostics use fusion_beam_fortran's public C-compatible type;
  !! no second definition of fusion_beam_window_v1 is introduced here.
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

  ! Exact C layout: fourteen consecutive c_double fields.
  type, bind(C), public :: fusion_target_burn_v1
     real(c_double) :: initial_fast_number_m3
     real(c_double) :: final_fast_number_m3
     real(c_double) :: initial_fast_energy_J_m3
     real(c_double) :: final_fast_energy_J_m3
     real(c_double) :: initial_target_number_m3
     real(c_double) :: final_target_number_m3
     real(c_double) :: initial_target_energy_J_m3
     real(c_double) :: final_target_energy_J_m3
     real(c_double) :: reactions_m3
     real(c_double) :: removed_fast_energy_J_m3
     real(c_double) :: removed_target_energy_J_m3
     real(c_double) :: fast_number_residual_m3
     real(c_double) :: target_number_residual_m3
     real(c_double) :: energy_residual_J_m3
  end type fusion_target_burn_v1

  public :: fusion_target_burn_trial
  public :: fusion_beam_target_burn_window_trial

  interface
     function c_fusion_target_burn_trial(cells, dt_s, edges, old_fast, &
          target_number, target_energy, reactivity, target_energy_reactivity, &
          trial_fast, reaction_loss, out) bind(C, &
          name="fusion_c_target_burn_trial") result(status)
       import :: c_double, c_int, c_ptr, fusion_target_burn_v1
       integer(c_int), value :: cells
       real(c_double), value :: dt_s
       type(c_ptr), value :: edges
       type(c_ptr), value :: old_fast
       real(c_double), value :: target_number
       real(c_double), value :: target_energy
       type(c_ptr), value :: reactivity
       type(c_ptr), value :: target_energy_reactivity
       type(c_ptr), value :: trial_fast
       type(c_ptr), value :: reaction_loss
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_target_burn_trial

     function c_fusion_beam_target_burn_window_trial(channel, cells, dt_s, &
          projectile_mass, target_mass, edges, old_fast, target_number, &
          target_kT, trial_fast, reaction_loss, windows, out) bind(C, &
          name="fusion_c_beam_target_burn_window_trial") result(status)
       import :: c_double, c_int, c_ptr, fusion_beam_window_v1, &
            fusion_target_burn_v1
       integer(c_int), value :: channel
       integer(c_int), value :: cells
       real(c_double), value :: dt_s
       real(c_double), value :: projectile_mass
       real(c_double), value :: target_mass
       type(c_ptr), value :: edges
       type(c_ptr), value :: old_fast
       real(c_double), value :: target_number
       real(c_double), value :: target_kT
       type(c_ptr), value :: trial_fast
       type(c_ptr), value :: reaction_loss
       type(c_ptr), value :: windows
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_target_burn_window_trial
  end interface

contains

  subroutine clear_target_burn(out)
    type(fusion_target_burn_v1), intent(out) :: out

    out%initial_fast_number_m3 = 0.0_c_double
    out%final_fast_number_m3 = 0.0_c_double
    out%initial_fast_energy_J_m3 = 0.0_c_double
    out%final_fast_energy_J_m3 = 0.0_c_double
    out%initial_target_number_m3 = 0.0_c_double
    out%final_target_number_m3 = 0.0_c_double
    out%initial_target_energy_J_m3 = 0.0_c_double
    out%final_target_energy_J_m3 = 0.0_c_double
    out%reactions_m3 = 0.0_c_double
    out%removed_fast_energy_J_m3 = 0.0_c_double
    out%removed_target_energy_J_m3 = 0.0_c_double
    out%fast_number_residual_m3 = 0.0_c_double
    out%target_number_residual_m3 = 0.0_c_double
    out%energy_residual_J_m3 = 0.0_c_double
  end subroutine clear_target_burn

  subroutine clear_windows(windows)
    type(fusion_beam_window_v1), intent(out) :: windows(:)

    windows%resolved_reactivity_m3_s = 0.0_c_double
    windows%projectile_energy_reactivity_J_m3_s = 0.0_c_double
    windows%target_energy_reactivity_J_m3_s = 0.0_c_double
    windows%relative_energy_reactivity_J_m3_s = 0.0_c_double
    windows%cm_energy_reactivity_J_m3_s = 0.0_c_double
    windows%resolved_pair_probability = 0.0_c_double
    windows%unresolved_pair_probability = 0.0_c_double
    windows%unresolved_relative_speed_m_s = 0.0_c_double
    windows%quadrature_error_m3_s = 0.0_c_double
    windows%energy_identity_error_J_m3_s = 0.0_c_double
    windows%domain_incomplete = 0_c_int
  end subroutine clear_windows

  subroutine fusion_target_burn_trial(cells, dt_s, edges_J, old_fast_m3, &
       target_number_m3, target_energy_J_m3, reactivity_m3_s, &
       target_energy_reactivity_J_m3_s, trial_fast_m3, reaction_loss_m3, &
       out, status)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in) :: dt_s, target_number_m3, target_energy_J_m3
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: old_fast_m3(:)
    real(c_double), intent(in), target, contiguous :: reactivity_m3_s(:)
    real(c_double), intent(in), target, contiguous :: &
         target_energy_reactivity_J_m3_s(:)
    real(c_double), intent(out), target, contiguous :: trial_fast_m3(:)
    real(c_double), intent(out), target, contiguous :: reaction_loss_m3(:)
    type(fusion_target_burn_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges

    call clear_target_burn(out)
    trial_fast_m3 = 0.0_c_double
    reaction_loss_m3 = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(old_fast_m3, kind=c_int) /= cells) return
    if (size(reactivity_m3_s, kind=c_int) /= cells) return
    if (size(target_energy_reactivity_J_m3_s, kind=c_int) /= cells) return
    if (size(trial_fast_m3, kind=c_int) /= cells) return
    if (size(reaction_loss_m3, kind=c_int) /= cells) return

    status = c_fusion_target_burn_trial(cells, dt_s, c_loc(edges_J(1)), &
         c_loc(old_fast_m3(1)), target_number_m3, target_energy_J_m3, &
         c_loc(reactivity_m3_s(1)), c_loc(target_energy_reactivity_J_m3_s(1)), &
         c_loc(trial_fast_m3(1)), c_loc(reaction_loss_m3(1)), c_loc(out))
    if (status /= PB11_STATUS_OK) then
       trial_fast_m3 = 0.0_c_double
       reaction_loss_m3 = 0.0_c_double
       call clear_target_burn(out)
    end if
  end subroutine fusion_target_burn_trial

  subroutine fusion_beam_target_burn_window_trial(channel, cells, dt_s, &
       projectile_mass_kg, target_mass_kg, edges_J, old_fast_m3, &
       target_number_m3, target_kT_J, trial_fast_m3, reaction_loss_m3, &
       windows, out, status)
    integer(c_int), intent(in) :: channel, cells
    real(c_double), intent(in) :: dt_s, projectile_mass_kg, target_mass_kg
    real(c_double), intent(in) :: target_number_m3, target_kT_J
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: old_fast_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_fast_m3(:)
    real(c_double), intent(out), target, contiguous :: reaction_loss_m3(:)
    type(fusion_beam_window_v1), intent(out), target, contiguous :: windows(:)
    type(fusion_target_burn_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges

    call clear_target_burn(out)
    trial_fast_m3 = 0.0_c_double
    reaction_loss_m3 = 0.0_c_double
    call clear_windows(windows)
    status = PB11_STATUS_INVALID_ARGUMENT

    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(old_fast_m3, kind=c_int) /= cells) return
    if (size(trial_fast_m3, kind=c_int) /= cells) return
    if (size(reaction_loss_m3, kind=c_int) /= cells) return
    if (size(windows, kind=c_int) /= cells) return

    status = c_fusion_beam_target_burn_window_trial(channel, cells, dt_s, &
         projectile_mass_kg, target_mass_kg, c_loc(edges_J(1)), &
         c_loc(old_fast_m3(1)), target_number_m3, target_kT_J, &
         c_loc(trial_fast_m3(1)), c_loc(reaction_loss_m3(1)), &
         c_loc(windows(1)), c_loc(out))
    if (status /= PB11_STATUS_OK) then
       trial_fast_m3 = 0.0_c_double
       reaction_loss_m3 = 0.0_c_double
       call clear_windows(windows)
       call clear_target_burn(out)
    end if
  end subroutine fusion_beam_target_burn_window_trial

end module fusion_target_burn_fortran
