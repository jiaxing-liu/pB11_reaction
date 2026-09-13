module fusion_two_component_fortran
  !! ISO_C_BINDING wrappers for the shared-grid two-component kinetic trial.
  !!
  !! The component state is deliberately an explicit old/trial pair.  The
  !! wrapper validates every Fortran extent before forming a C address and
  !! passes C_NULL_PTR for zero-length bath/diffusion arrays.  It clears all
  !! returned state and ledger fields if the C routine rejects the request.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_ptr, &
       c_loc
  use fusion_kinetics_fortran, only : fusion_kinetic_ledger_v1, &
       fusion_maxwellian_bath_v1, PB11_STATUS_OK, PB11_STATUS_INVALID_ARGUMENT
  implicit none
  private

  ! Exact C layout: fusion_kinetic_ledger_v1 followed by two doubles.
  type, bind(C), public :: fusion_two_component_ledger_v1
     type(fusion_kinetic_ledger_v1) :: total
     real(c_double) :: transferred_number_m3
     real(c_double) :: transferred_energy_J_m3
  end type fusion_two_component_ledger_v1

  public :: fusion_coulomb_transfer_rate
  public :: fusion_two_component_trial

  interface
     function c_fusion_coulomb_transfer_rate(energy_J, mass_kg, &
          charge_number, ion_bath, rate_s_inv) bind(C, &
          name="fusion_c_coulomb_transfer_rate") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: energy_J
       real(c_double), value :: mass_kg
       real(c_double), value :: charge_number
       type(c_ptr), value :: ion_bath
       type(c_ptr), value :: rate_s_inv
       integer(c_int) :: status
     end function c_fusion_coulomb_transfer_rate

     function c_fusion_two_component_trial(cells, baths, dt_s, edges, &
          old_s, old_t, bath_kT, diffusion, birth_s, birth_t, escape, &
          transfer, trial_s, trial_t, heat_to_bath, ledger) bind(C, &
          name="fusion_c_two_component_trial") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: cells
       integer(c_int), value :: baths
       real(c_double), value :: dt_s
       type(c_ptr), value :: edges
       type(c_ptr), value :: old_s
       type(c_ptr), value :: old_t
       type(c_ptr), value :: bath_kT
       type(c_ptr), value :: diffusion
       type(c_ptr), value :: birth_s
       type(c_ptr), value :: birth_t
       type(c_ptr), value :: escape
       type(c_ptr), value :: transfer
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       type(c_ptr), value :: heat_to_bath
       type(c_ptr), value :: ledger
       integer(c_int) :: status
     end function c_fusion_two_component_trial
  end interface

contains

  subroutine clear_ledger(ledger)
    type(fusion_two_component_ledger_v1), intent(out) :: ledger

    ledger%total%initial_number_m3 = 0.0_c_double
    ledger%total%final_number_m3 = 0.0_c_double
    ledger%total%initial_energy_J_m3 = 0.0_c_double
    ledger%total%final_energy_J_m3 = 0.0_c_double
    ledger%total%born_number_m3 = 0.0_c_double
    ledger%total%born_energy_J_m3 = 0.0_c_double
    ledger%total%escaped_number_m3 = 0.0_c_double
    ledger%total%escaped_energy_J_m3 = 0.0_c_double
    ledger%total%thermalized_number_m3 = 0.0_c_double
    ledger%total%thermalized_energy_J_m3 = 0.0_c_double
    ledger%total%particle_balance_error_m3 = 0.0_c_double
    ledger%total%energy_balance_error_J_m3 = 0.0_c_double
    ledger%transferred_number_m3 = 0.0_c_double
    ledger%transferred_energy_J_m3 = 0.0_c_double
  end subroutine clear_ledger

  subroutine fusion_coulomb_transfer_rate(energy_J, mass_kg, charge_number, &
       ion_bath, rate_s_inv, status)
    real(c_double), intent(in) :: energy_J, mass_kg, charge_number
    type(fusion_maxwellian_bath_v1), intent(in), target :: ion_bath
    real(c_double), intent(out), target :: rate_s_inv
    integer(c_int), intent(out) :: status

    rate_s_inv = 0.0_c_double
    status = c_fusion_coulomb_transfer_rate(energy_J, mass_kg, charge_number, &
         c_loc(ion_bath), c_loc(rate_s_inv))
    if (status /= PB11_STATUS_OK) rate_s_inv = 0.0_c_double
  end subroutine fusion_coulomb_transfer_rate

  subroutine fusion_two_component_trial(cells, baths, dt_s, edges_J, &
       old_s_m3, old_t_m3, bath_kT_J, diffusion_J2_s, birth_s_m3_s, &
       birth_t_m3_s, escape_s_inv, transfer_s_inv, trial_s_m3, trial_t_m3, &
       heat_to_bath_J_m3, ledger, status)
    integer(c_int), intent(in) :: cells, baths
    real(c_double), intent(in) :: dt_s
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: old_s_m3(:)
    real(c_double), intent(in), target, contiguous :: old_t_m3(:)
    real(c_double), intent(in), target, contiguous :: bath_kT_J(:)
    real(c_double), intent(in), target, contiguous :: diffusion_J2_s(:)
    real(c_double), intent(in), target, contiguous :: birth_s_m3_s(:)
    real(c_double), intent(in), target, contiguous :: birth_t_m3_s(:)
    real(c_double), intent(in), target, contiguous :: escape_s_inv(:)
    real(c_double), intent(in), target, contiguous :: transfer_s_inv(:)
    real(c_double), intent(out), target, contiguous :: trial_s_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_t_m3(:)
    real(c_double), intent(out), target, contiguous :: heat_to_bath_J_m3(:)
    type(fusion_two_component_ledger_v1), intent(out), target :: ledger
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges, expected_diff
    type(c_ptr) :: edges_ptr, old_s_ptr, old_t_ptr, bath_ptr
    type(c_ptr) :: diffusion_ptr, birth_s_ptr, birth_t_ptr, escape_ptr
    type(c_ptr) :: transfer_ptr, trial_s_ptr, trial_t_ptr, heat_ptr
    type(c_ptr) :: ledger_ptr

    call clear_ledger(ledger)
    trial_s_m3 = 0.0_c_double
    trial_t_m3 = 0.0_c_double
    heat_to_bath_J_m3 = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Validate dimensions using c_int before any element is referenced.  This
    ! also prevents a malformed C-int dimension from wrapping an extent.
    if (cells < 1_c_int .or. baths < 0_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(old_s_m3, kind=c_int) /= cells) return
    if (size(old_t_m3, kind=c_int) /= cells) return
    if (size(birth_s_m3_s, kind=c_int) /= cells) return
    if (size(birth_t_m3_s, kind=c_int) /= cells) return
    if (size(escape_s_inv, kind=c_int) /= cells) return
    if (size(transfer_s_inv, kind=c_int) /= cells) return
    if (size(trial_s_m3, kind=c_int) /= cells) return
    if (size(trial_t_m3, kind=c_int) /= cells) return
    if (size(bath_kT_J, kind=c_int) /= baths) return
    if (size(heat_to_bath_J_m3, kind=c_int) /= baths) return

    expected_diff = 0_c_int
    if (cells > 1_c_int) then
       if (baths > huge(expected_diff) / (cells - 1_c_int)) return
       expected_diff = baths * (cells - 1_c_int)
    end if
    if (size(diffusion_J2_s, kind=c_int) /= expected_diff) return

    ! Every cell-based vector has at least one element after the checks above.
    edges_ptr = c_loc(edges_J(1))
    old_s_ptr = c_loc(old_s_m3(1))
    old_t_ptr = c_loc(old_t_m3(1))
    birth_s_ptr = c_loc(birth_s_m3_s(1))
    birth_t_ptr = c_loc(birth_t_m3_s(1))
    escape_ptr = c_loc(escape_s_inv(1))
    transfer_ptr = c_loc(transfer_s_inv(1))
    trial_s_ptr = c_loc(trial_s_m3(1))
    trial_t_ptr = c_loc(trial_t_m3(1))
    ledger_ptr = c_loc(ledger)

    ! Zero-length arrays have no element one.  C_NULL_PTR is the required
    ! representation for baths=0 and for diffusion when cells=1.
    bath_ptr = c_null_ptr
    diffusion_ptr = c_null_ptr
    heat_ptr = c_null_ptr
    if (baths > 0_c_int) then
       bath_ptr = c_loc(bath_kT_J(1))
       heat_ptr = c_loc(heat_to_bath_J_m3(1))
    end if
    if (expected_diff > 0_c_int) diffusion_ptr = c_loc(diffusion_J2_s(1))

    status = c_fusion_two_component_trial(cells, baths, dt_s, edges_ptr, &
         old_s_ptr, old_t_ptr, bath_ptr, diffusion_ptr, birth_s_ptr, &
         birth_t_ptr, escape_ptr, transfer_ptr, trial_s_ptr, trial_t_ptr, &
         heat_ptr, ledger_ptr)
    if (status /= PB11_STATUS_OK) then
       trial_s_m3 = 0.0_c_double
       trial_t_m3 = 0.0_c_double
       heat_to_bath_J_m3 = 0.0_c_double
       call clear_ledger(ledger)
    end if
  end subroutine fusion_two_component_trial

end module fusion_two_component_fortran
