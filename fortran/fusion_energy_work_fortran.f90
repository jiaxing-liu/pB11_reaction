module fusion_energy_work_fortran
  !! ISO_C_BINDING wrapper for the stateless energy-space work trial.
  !!
  !! The C layer owns the backward-Euler upwind solve.  This module exposes
  !! the C ledger with its exact field order and checks every assumed-shape
  !! extent before forming a C address.  Every output is cleared on rejection.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_size_t, &
       c_loc
  implicit none
  private

  ! Keep this binding standalone; these values are the common pb11 C ABI
  ! status values used by the other Fortran wrappers.
  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: FUSION_ENERGY_WORK_MAX_CELLS = &
       1000000_c_int

  ! Exact C layout: eleven consecutive c_double fields (88 bytes).
  type, bind(C), public :: fusion_energy_work_ledger_v1
     real(c_double) :: initial_number_m3
     real(c_double) :: final_number_m3
     real(c_double) :: initial_energy_J_m3
     real(c_double) :: final_energy_J_m3
     real(c_double) :: lower_number_m3
     real(c_double) :: lower_energy_J_m3
     real(c_double) :: upper_number_m3
     real(c_double) :: upper_energy_J_m3
     real(c_double) :: work_on_particles_J_m3
     real(c_double) :: particle_balance_error_m3
     real(c_double) :: energy_balance_error_J_m3
  end type fusion_energy_work_ledger_v1

  public :: fusion_energy_work_trial

  interface
     function c_fusion_energy_work_trial(cells, dt_s, compression_s_inv, &
          edges, old_number, trial_number, face_amount, ledger) bind(C, &
          name="fusion_c_energy_work_trial") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: cells
       real(c_double), value :: dt_s
       real(c_double), value :: compression_s_inv
       type(c_ptr), value :: edges
       type(c_ptr), value :: old_number
       type(c_ptr), value :: trial_number
       type(c_ptr), value :: face_amount
       type(c_ptr), value :: ledger
       integer(c_int) :: status
     end function c_fusion_energy_work_trial
  end interface

contains

  subroutine clear_ledger(ledger)
    type(fusion_energy_work_ledger_v1), intent(out) :: ledger

    ledger%initial_number_m3 = 0.0_c_double
    ledger%final_number_m3 = 0.0_c_double
    ledger%initial_energy_J_m3 = 0.0_c_double
    ledger%final_energy_J_m3 = 0.0_c_double
    ledger%lower_number_m3 = 0.0_c_double
    ledger%lower_energy_J_m3 = 0.0_c_double
    ledger%upper_number_m3 = 0.0_c_double
    ledger%upper_energy_J_m3 = 0.0_c_double
    ledger%work_on_particles_J_m3 = 0.0_c_double
    ledger%particle_balance_error_m3 = 0.0_c_double
    ledger%energy_balance_error_J_m3 = 0.0_c_double
  end subroutine clear_ledger

  subroutine fusion_energy_work_trial(cells, dt_s, compression_s_inv, &
       edges_J, old_number_m3, trial_number_m3, face_amount_m3, ledger, &
       status)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in) :: dt_s, compression_s_inv
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: old_number_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_number_m3(:)
    real(c_double), intent(out), target, contiguous :: face_amount_m3(:)
    type(fusion_energy_work_ledger_v1), intent(out), target :: ledger
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: cell_count, edge_count

    ! Clear every output before validating dimensions.  A malformed call
    ! returns before any C address is formed, including for zero-size arrays.
    trial_number_m3 = 0.0_c_double
    face_amount_m3 = 0.0_c_double
    call clear_ledger(ledger)
    status = PB11_STATUS_INVALID_ARGUMENT

    if (cells < 1_c_int .or. cells > FUSION_ENERGY_WORK_MAX_CELLS) return

    cell_count = int(cells, c_size_t)
    edge_count = cell_count + 1_c_size_t
    if (size(edges_J, kind=c_size_t) /= edge_count) return
    if (size(old_number_m3, kind=c_size_t) /= cell_count) return
    if (size(trial_number_m3, kind=c_size_t) /= cell_count) return
    if (size(face_amount_m3, kind=c_size_t) /= edge_count) return

    status = c_fusion_energy_work_trial(cells, dt_s, compression_s_inv, &
         c_loc(edges_J(1)), c_loc(old_number_m3(1)), &
         c_loc(trial_number_m3(1)), c_loc(face_amount_m3(1)), &
         c_loc(ledger))
    if (status /= PB11_STATUS_OK) then
       trial_number_m3 = 0.0_c_double
       face_amount_m3 = 0.0_c_double
       call clear_ledger(ledger)
    end if
  end subroutine fusion_energy_work_trial

end module fusion_energy_work_fortran
