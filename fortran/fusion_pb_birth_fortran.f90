module fusion_pb_birth_fortran
  !! ISO_C_BINDING wrappers for the p-11B CM birth-source APIs.
  !!
  !! The C API owns the physical source model.  These entry points only
  !! provide an explicit Fortran kind/layout contract and validate assumed-
  !! shape array extents before forming C pointers.
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

  ! Exact C layout: fourteen consecutive c_double fields.
  type, bind(C), public :: fusion_pb_birth_v1
     real(c_double) :: mapped_number
     real(c_double) :: mapped_energy_J
     real(c_double) :: below_number
     real(c_double) :: below_energy_J
     real(c_double) :: above_number
     real(c_double) :: above_energy_J
     real(c_double) :: number_residual
     real(c_double) :: energy_residual_J
     real(c_double) :: alpha0_fraction
     real(c_double) :: low_alpha1_fraction
     real(c_double) :: broad_alpha1_fraction
     real(c_double) :: primary_alpha0_energy_J
     real(c_double) :: secondary_min_J
     real(c_double) :: secondary_max_J
  end type fusion_pb_birth_v1

  public :: fusion_pb_alpha0_grid
  public :: fusion_pb_cm_source_grid

  interface
     function c_fusion_pb_alpha0_grid(A_J, q_J, cells, edges, birth, out) &
          bind(C, name="fusion_c_pb_alpha0_grid") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: A_J
       real(c_double), value :: q_J
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_pb_alpha0_grid

     function c_fusion_pb_cm_source_grid(A_J, q_J, f0, flow, broad_mode, &
          fsci_policy, cutoff_J, k, phase, nq, ncos, cells, edges, birth, &
          out) bind(C, name="fusion_c_pb_cm_source_grid") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: A_J
       real(c_double), value :: q_J
       real(c_double), value :: f0
       real(c_double), value :: flow
       integer(c_int), value :: broad_mode
       integer(c_int), value :: fsci_policy
       real(c_double), value :: cutoff_J
       real(c_double), value :: k
       real(c_double), value :: phase
       integer(c_int), value :: nq
       integer(c_int), value :: ncos
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_pb_cm_source_grid
  end interface

contains

  subroutine clear_birth(out)
    type(fusion_pb_birth_v1), intent(out) :: out

    out%mapped_number = 0.0_c_double
    out%mapped_energy_J = 0.0_c_double
    out%below_number = 0.0_c_double
    out%below_energy_J = 0.0_c_double
    out%above_number = 0.0_c_double
    out%above_energy_J = 0.0_c_double
    out%number_residual = 0.0_c_double
    out%energy_residual_J = 0.0_c_double
    out%alpha0_fraction = 0.0_c_double
    out%low_alpha1_fraction = 0.0_c_double
    out%broad_alpha1_fraction = 0.0_c_double
    out%primary_alpha0_energy_J = 0.0_c_double
    out%secondary_min_J = 0.0_c_double
    out%secondary_max_J = 0.0_c_double
  end subroutine clear_birth

  subroutine fusion_pb_alpha0_grid(A_J, q_J, cells, edges_J, birth, out, &
       status)
    real(c_double), intent(in) :: A_J, q_J
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(out), target, contiguous :: birth(:)
    type(fusion_pb_birth_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges
    type(c_ptr) :: edges_ptr, birth_ptr, out_ptr

    call clear_birth(out)
    birth = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Validate dimensions before C_LOC.  The C API has no zero-cell pointer
    ! form, so the first elements are formed only after these checks pass.
    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(birth, kind=c_int) /= cells) return

    edges_ptr = c_loc(edges_J(1))
    birth_ptr = c_loc(birth(1))
    out_ptr = c_loc(out)
    status = c_fusion_pb_alpha0_grid(A_J, q_J, cells, edges_ptr, birth_ptr, &
         out_ptr)
    if (status /= PB11_STATUS_OK) then
       birth = 0.0_c_double
       call clear_birth(out)
    end if
  end subroutine fusion_pb_alpha0_grid

  subroutine fusion_pb_cm_source_grid(A_J, q_J, f0, flow, broad_mode, &
       fsci_policy, cutoff_J, k, phase, nq, ncos, cells, edges_J, birth, &
       out, status)
    real(c_double), intent(in) :: A_J, q_J, f0, flow
    integer(c_int), intent(in) :: broad_mode, fsci_policy
    real(c_double), intent(in) :: cutoff_J, k, phase
    integer(c_int), intent(in) :: nq, ncos, cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(out), target, contiguous :: birth(:)
    type(fusion_pb_birth_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges
    type(c_ptr) :: edges_ptr, birth_ptr, out_ptr

    call clear_birth(out)
    birth = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Validate dimensions before C_LOC; all source controls are validated by
    ! the C implementation, including controls for zero-weight branches.
    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(birth, kind=c_int) /= cells) return

    edges_ptr = c_loc(edges_J(1))
    birth_ptr = c_loc(birth(1))
    out_ptr = c_loc(out)
    status = c_fusion_pb_cm_source_grid(A_J, q_J, f0, flow, broad_mode, &
         fsci_policy, cutoff_J, k, phase, nq, ncos, cells, edges_ptr, &
         birth_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       birth = 0.0_c_double
       call clear_birth(out)
    end if
  end subroutine fusion_pb_cm_source_grid

end module fusion_pb_birth_fortran
