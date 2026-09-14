module fusion_target_network_fortran
  !! ISO_C_BINDING wrapper for the bounded fast/target network trial.
  !!
  !! The C layer owns the numerical solve.  This module mirrors the C
  !! records exactly, infers all C dimensions from the Fortran arrays, and
  !! validates every extent before taking C_LOC.  Zero-edge calls pass a
  !! C_NULL_PTR for the optional edge and event buffers.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_ptr, &
       c_null_ptr, c_loc
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: &
       FUSION_TARGET_NETWORK_MAX_FAST = 600000_c_int
  integer(c_int), parameter, public :: &
       FUSION_TARGET_NETWORK_MAX_TARGET = 6_c_int
  integer(c_int), parameter, public :: &
       FUSION_TARGET_NETWORK_MAX_EDGE = 2000000_c_int

  ! Exact C layout: two C ints followed by two C doubles.
  type, bind(C), public :: fusion_target_network_edge_v1
     integer(c_int) :: fast_index
     integer(c_int) :: target_index
     real(c_double) :: reactivity_m3_s
     real(c_double) :: target_energy_reactivity_J_m3_s
  end type fusion_target_network_edge_v1

  ! Exact C layout: six C doubles followed by one C int.
  type, bind(C), public :: fusion_target_network_v1
     real(c_double) :: reactions_m3
     real(c_double) :: removed_fast_energy_J_m3
     real(c_double) :: removed_target_energy_J_m3
     real(c_double) :: max_fast_number_relative_residual
     real(c_double) :: max_target_number_relative_residual
     real(c_double) :: max_target_energy_relative_residual
     integer(c_int) :: iterations
  end type fusion_target_network_v1

  public :: target_network_trial

  interface
     function c_fusion_target_network_trial(nfast, ntarget, nedge, dt_s, &
          fast_energy_J, old_fast_m3, target_number_m3, &
          target_energy_J_m3, edges, trial_fast_m3, trial_target_m3, &
          trial_target_energy_J_m3, edge_reactions_m3, out) bind(C, &
          name="fusion_c_target_network_trial") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: nfast
       integer(c_int), value :: ntarget
       integer(c_int), value :: nedge
       real(c_double), value :: dt_s
       type(c_ptr), value :: fast_energy_J
       type(c_ptr), value :: old_fast_m3
       type(c_ptr), value :: target_number_m3
       type(c_ptr), value :: target_energy_J_m3
       type(c_ptr), value :: edges
       type(c_ptr), value :: trial_fast_m3
       type(c_ptr), value :: trial_target_m3
       type(c_ptr), value :: trial_target_energy_J_m3
       type(c_ptr), value :: edge_reactions_m3
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_target_network_trial
  end interface

contains

  subroutine clear_network(out)
    type(fusion_target_network_v1), intent(out) :: out

    out%reactions_m3 = 0.0_c_double
    out%removed_fast_energy_J_m3 = 0.0_c_double
    out%removed_target_energy_J_m3 = 0.0_c_double
    out%max_fast_number_relative_residual = 0.0_c_double
    out%max_target_number_relative_residual = 0.0_c_double
    out%max_target_energy_relative_residual = 0.0_c_double
    out%iterations = 0_c_int
  end subroutine clear_network

  subroutine target_network_trial(dt_s, fast_energy_J, old_fast_m3, &
       target_number_m3, target_energy_J_m3, edges, trial_fast_m3, &
       trial_target_m3, trial_target_energy_J_m3, edge_reactions_m3, out, &
       status)
    real(c_double), intent(in) :: dt_s
    real(c_double), intent(in), target, contiguous :: fast_energy_J(:)
    real(c_double), intent(in), target, contiguous :: old_fast_m3(:)
    real(c_double), intent(in), target, contiguous :: target_number_m3(:)
    real(c_double), intent(in), target, contiguous :: target_energy_J_m3(:)
    type(fusion_target_network_edge_v1), intent(in), target, contiguous :: &
         edges(:)
    real(c_double), intent(out), target, contiguous :: trial_fast_m3(:)
    real(c_double), intent(out), target, contiguous :: trial_target_m3(:)
    real(c_double), intent(out), target, contiguous :: &
         trial_target_energy_J_m3(:)
    real(c_double), intent(out), target, contiguous :: edge_reactions_m3(:)
    type(fusion_target_network_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: nfast_size, ntarget_size, nedge_size
    integer(c_int) :: nfast, ntarget, nedge
    type(c_ptr) :: fast_energy_ptr, old_fast_ptr, target_number_ptr
    type(c_ptr) :: target_energy_ptr, edges_ptr
    type(c_ptr) :: trial_fast_ptr, trial_target_ptr, trial_energy_ptr
    type(c_ptr) :: edge_reactions_ptr, out_ptr

    ! Clear every output before inspecting dimensions.  This also covers
    ! zero-length output actuals and all shape-rejection paths.
    call clear_network(out)
    trial_fast_m3 = 0.0_c_double
    trial_target_m3 = 0.0_c_double
    trial_target_energy_J_m3 = 0.0_c_double
    edge_reactions_m3 = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Infer the C dimensions and check all limits/extents before C_LOC.  The
    ! C contract permits no zero fast/target dimensions and permits zero
    ! links, whose pointers are set to C_NULL_PTR below.
    nfast_size = size(fast_energy_J, kind=c_size_t)
    ntarget_size = size(target_number_m3, kind=c_size_t)
    nedge_size = size(edges, kind=c_size_t)
    if (nfast_size < 1_c_size_t .or. nfast_size > &
         int(FUSION_TARGET_NETWORK_MAX_FAST, kind=c_size_t)) return
    if (ntarget_size < 1_c_size_t .or. ntarget_size > &
         int(FUSION_TARGET_NETWORK_MAX_TARGET, kind=c_size_t)) return
    if (nedge_size > int(FUSION_TARGET_NETWORK_MAX_EDGE, kind=c_size_t)) return
    if (size(old_fast_m3, kind=c_size_t) /= nfast_size) return
    if (size(trial_fast_m3, kind=c_size_t) /= nfast_size) return
    if (size(target_energy_J_m3, kind=c_size_t) /= ntarget_size) return
    if (size(trial_target_m3, kind=c_size_t) /= ntarget_size) return
    if (size(trial_target_energy_J_m3, kind=c_size_t) /= ntarget_size) return
    if (size(edge_reactions_m3, kind=c_size_t) /= nedge_size) return

    ! The bounds above make these narrowing conversions representable as C
    ! ints, while c_size_t protects the pre-C_LOC extent checks from wrap.
    nfast = int(nfast_size, kind=c_int)
    ntarget = int(ntarget_size, kind=c_int)
    nedge = int(nedge_size, kind=c_int)

    fast_energy_ptr = c_loc(fast_energy_J(1))
    old_fast_ptr = c_loc(old_fast_m3(1))
    target_number_ptr = c_loc(target_number_m3(1))
    target_energy_ptr = c_loc(target_energy_J_m3(1))
    trial_fast_ptr = c_loc(trial_fast_m3(1))
    trial_target_ptr = c_loc(trial_target_m3(1))
    trial_energy_ptr = c_loc(trial_target_energy_J_m3(1))
    out_ptr = c_loc(out)

    edges_ptr = c_null_ptr
    edge_reactions_ptr = c_null_ptr
    if (nedge > 0_c_int) then
       edges_ptr = c_loc(edges(1))
       edge_reactions_ptr = c_loc(edge_reactions_m3(1))
    end if

    status = c_fusion_target_network_trial(nfast, ntarget, nedge, dt_s, &
         fast_energy_ptr, old_fast_ptr, target_number_ptr, target_energy_ptr, &
         edges_ptr, trial_fast_ptr, trial_target_ptr, trial_energy_ptr, &
         edge_reactions_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       trial_fast_m3 = 0.0_c_double
       trial_target_m3 = 0.0_c_double
       trial_target_energy_J_m3 = 0.0_c_double
       edge_reactions_m3 = 0.0_c_double
       call clear_network(out)
    end if
  end subroutine target_network_trial

end module fusion_target_network_fortran
