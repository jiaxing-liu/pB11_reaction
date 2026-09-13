module fusion_handoff_fortran
  !! ISO_C_BINDING wrappers for Maxwellian energy-grid construction and the
  !! explicit kinetic-to-fluid handoff trial.
  !!
  !! Both public routines validate Fortran array extents before taking C_LOC.
  !! They clear all output arrays, scalar flags, and interoperable structures
  !! on rejected calls, so a caller can safely discard a failed trial.
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

  ! Exact C layout: four consecutive c_double fields.
  type, bind(C), public :: fusion_maxwellian_grid_v1
     real(c_double) :: below_probability
     real(c_double) :: above_probability
     real(c_double) :: represented_mean_energy_J
     real(c_double) :: probability_balance_error
  end type fusion_maxwellian_grid_v1

  ! Exact C layout: thirteen consecutive c_double fields.
  type, bind(C), public :: fusion_handoff_ledger_v1
     real(c_double) :: initial_number_m3
     real(c_double) :: initial_energy_J_m3
     real(c_double) :: remaining_number_m3
     real(c_double) :: remaining_energy_J_m3
     real(c_double) :: fluid_number_m3
     real(c_double) :: fluid_energy_J_m3
     real(c_double) :: bath_energy_correction_J_m3
     real(c_double) :: distribution_L1
     real(c_double) :: relative_mean_energy_error
     real(c_double) :: outside_grid_probability
     real(c_double) :: represented_Maxwellian_mean_energy_J
     real(c_double) :: particle_balance_error_m3
     real(c_double) :: energy_balance_error_J_m3
  end type fusion_handoff_ledger_v1

  public :: fusion_maxwellian_energy_grid
  public :: fusion_maxwellian_handoff_trial

  interface
     function c_fusion_maxwellian_energy_grid(cells, kT_J, edges, probability, &
          out) bind(C, name="fusion_c_maxwellian_energy_grid") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: cells
       real(c_double), value :: kT_J
       type(c_ptr), value :: edges
       type(c_ptr), value :: probability
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_maxwellian_energy_grid

     function c_fusion_maxwellian_handoff_trial(cells, target_kT_J, max_L1, &
          max_relative_mean_error, edges, old_number, trial_number, projected, &
          out) bind(C, name="fusion_c_maxwellian_handoff_trial") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: cells
       real(c_double), value :: target_kT_J
       real(c_double), value :: max_L1
       real(c_double), value :: max_relative_mean_error
       type(c_ptr), value :: edges
       type(c_ptr), value :: old_number
       type(c_ptr), value :: trial_number
       type(c_ptr), value :: projected
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_maxwellian_handoff_trial
  end interface

contains

  subroutine clear_grid(out)
    type(fusion_maxwellian_grid_v1), intent(out) :: out

    out%below_probability = 0.0_c_double
    out%above_probability = 0.0_c_double
    out%represented_mean_energy_J = 0.0_c_double
    out%probability_balance_error = 0.0_c_double
  end subroutine clear_grid

  subroutine clear_handoff(out)
    type(fusion_handoff_ledger_v1), intent(out) :: out

    out%initial_number_m3 = 0.0_c_double
    out%initial_energy_J_m3 = 0.0_c_double
    out%remaining_number_m3 = 0.0_c_double
    out%remaining_energy_J_m3 = 0.0_c_double
    out%fluid_number_m3 = 0.0_c_double
    out%fluid_energy_J_m3 = 0.0_c_double
    out%bath_energy_correction_J_m3 = 0.0_c_double
    out%distribution_L1 = 0.0_c_double
    out%relative_mean_energy_error = 0.0_c_double
    out%outside_grid_probability = 0.0_c_double
    out%represented_Maxwellian_mean_energy_J = 0.0_c_double
    out%particle_balance_error_m3 = 0.0_c_double
    out%energy_balance_error_J_m3 = 0.0_c_double
  end subroutine clear_handoff

  subroutine fusion_maxwellian_energy_grid(cells, kT_J, edges_J, probability, &
       out, status)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in) :: kT_J
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(out), target, contiguous :: probability(:)
    type(fusion_maxwellian_grid_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges
    type(c_ptr) :: edges_ptr, probability_ptr, out_ptr

    call clear_grid(out)
    probability = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Check the C-int arithmetic and both extents before referring to element
    ! one.  The C API has no valid zero-cell representation.
    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(probability, kind=c_int) /= cells) return

    edges_ptr = c_loc(edges_J(1))
    probability_ptr = c_loc(probability(1))
    out_ptr = c_loc(out)
    status = c_fusion_maxwellian_energy_grid(cells, kT_J, edges_ptr, &
         probability_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       probability = 0.0_c_double
       call clear_grid(out)
    end if
  end subroutine fusion_maxwellian_energy_grid

  subroutine fusion_maxwellian_handoff_trial(cells, target_kT_J, max_L1, &
       max_relative_mean_error, edges_J, old_number_m3, trial_number_m3, &
       projected, out, status)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in) :: target_kT_J, max_L1
    real(c_double), intent(in) :: max_relative_mean_error
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: old_number_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_number_m3(:)
    integer(c_int), intent(out), target :: projected
    type(fusion_handoff_ledger_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges
    type(c_ptr) :: edges_ptr, old_ptr, trial_ptr, projected_ptr, out_ptr

    call clear_handoff(out)
    trial_number_m3 = 0.0_c_double
    projected = 0_c_int
    status = PB11_STATUS_INVALID_ARGUMENT

    ! The C API requires one or more cell values.  Validate all extents before
    ! C_LOC so a malformed Fortran call cannot form an invalid address.
    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(old_number_m3, kind=c_int) /= cells) return
    if (size(trial_number_m3, kind=c_int) /= cells) return

    edges_ptr = c_loc(edges_J(1))
    old_ptr = c_loc(old_number_m3(1))
    trial_ptr = c_loc(trial_number_m3(1))
    projected_ptr = c_loc(projected)
    out_ptr = c_loc(out)
    status = c_fusion_maxwellian_handoff_trial(cells, target_kT_J, max_L1, &
         max_relative_mean_error, edges_ptr, old_ptr, trial_ptr, projected_ptr, &
         out_ptr)
    if (status /= PB11_STATUS_OK) then
       trial_number_m3 = 0.0_c_double
       projected = 0_c_int
       call clear_handoff(out)
    end if
  end subroutine fusion_maxwellian_handoff_trial

end module fusion_handoff_fortran
