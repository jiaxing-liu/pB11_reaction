module fusion_source_state_fortran
  !! ISO_C_BINDING wrappers for the caller-owned source-state context.
  !!
  !! The opaque C context is represented directly by type(c_ptr).  This
  !! module has no finalizer: the caller owns every successful context and
  !! must call fusion_source_state_destroy exactly once.  Intrinsic
  !! assignment of a c_ptr copies only the opaque handle; it does not transfer
  !! ownership, so copied handles must not be destroyed independently.
  !!
  !! The C context stores its immutable cell count.  Snapshot and stage query
  !! that count before checking array extents and before forming any C_LOC.
  !! Restart unpack likewise queries the newly-created context, so no caller
  !! supplied dimension is trusted for an opaque state.
  !! The extended stage/snapshot pair carries one signed finite inert-bath
  !! heat amount per fast species.  It is separate from the fixed seven-bath
  !! ledger fields and is preserved by version-2 restart packets.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int8_t, &
       c_int64_t, c_size_t, c_ptr, c_null_ptr, c_loc, c_associated
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: FUSION_SOURCE_STATE_SPECIES = 6_c_int
  integer(c_int), parameter, public :: FUSION_SOURCE_STATE_MAX_CELLS = &
       1000000_c_int
  integer(c_int), parameter, public :: FUSION_SOURCE_STATE_LEDGER_DOUBLES = &
       121_c_int

  ! Exact C layout: 121 consecutive c_double fields.  C arrays are traversed
  ! in their natural order, and each Fortran component is one such array.
  type, bind(C), public :: fusion_source_ledger_v1
     real(c_double) :: events_m3(5)
     real(c_double) :: nuclear_born_number_m3(6)
     real(c_double) :: nuclear_born_energy_J_m3(6)
     real(c_double) :: external_born_number_m3(6)
     real(c_double) :: external_born_energy_J_m3(6)
     real(c_double) :: thermal_consumed_number_m3(6)
     real(c_double) :: thermal_consumed_energy_J_m3(6)
     real(c_double) :: fast_consumed_number_m3(6)
     real(c_double) :: fast_consumed_energy_J_m3(6)
     real(c_double) :: escaped_number_m3(6)
     real(c_double) :: escaped_energy_J_m3(6)
     real(c_double) :: handed_off_number_m3(6)
     real(c_double) :: handed_off_energy_J_m3(6)
     real(c_double) :: heat_to_bath_J_m3(42)
     real(c_double) :: neutron_number_m3
     real(c_double) :: neutron_energy_J_m3
  end type fusion_source_ledger_v1

  public :: fusion_source_state_create
  public :: fusion_source_state_destroy
  public :: fusion_source_state_cells
  public :: fusion_source_state_snapshot
  public :: fusion_source_state_snapshot_inert
  public :: fusion_source_state_begin
  public :: fusion_source_state_stage
  public :: fusion_source_state_stage_inert
  public :: fusion_source_state_commit
  public :: fusion_source_state_discard
  public :: fusion_source_state_pack_size
  public :: fusion_source_state_pack
  public :: fusion_source_state_unpack

  interface
     function c_fusion_source_state_create(cells, edges, initial_s, initial_t, &
          initial_time, model_tag, out) bind(C, &
          name="fusion_c_source_state_create") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: initial_s
       type(c_ptr), value :: initial_t
       real(c_double), value :: initial_time
       integer(c_int64_t), value :: model_tag
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_source_state_create

     subroutine c_fusion_source_state_destroy(state) bind(C, &
          name="fusion_c_source_state_destroy")
       import :: c_ptr
       type(c_ptr), value :: state
     end subroutine c_fusion_source_state_destroy

     function c_fusion_source_state_cells(state, cells) bind(C, &
          name="fusion_c_source_state_cells") result(status)
       import :: c_int, c_ptr
       type(c_ptr), value :: state
       type(c_ptr), value :: cells
       integer(c_int) :: status
     end function c_fusion_source_state_cells

     function c_fusion_source_state_snapshot(state, accepted_s, accepted_t, &
          cumulative, accepted_time, epoch) bind(C, &
          name="fusion_c_source_state_snapshot") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       type(c_ptr), value :: accepted_s
       type(c_ptr), value :: accepted_t
       type(c_ptr), value :: cumulative
       type(c_ptr), value :: accepted_time
       type(c_ptr), value :: epoch
       integer(c_int) :: status
     end function c_fusion_source_state_snapshot

     function c_fusion_source_state_begin(state, dt_s, ticket) bind(C, &
          name="fusion_c_source_state_begin") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       real(c_double), value :: dt_s
       type(c_ptr), value :: ticket
       integer(c_int) :: status
     end function c_fusion_source_state_begin

     function c_fusion_source_state_stage(state, ticket, trial_s, trial_t, &
          step) bind(C, name="fusion_c_source_state_stage") result(status)
       import :: c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       integer(c_int64_t), value :: ticket
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       type(c_ptr), value :: step
       integer(c_int) :: status
     end function c_fusion_source_state_stage

     function c_fusion_source_state_stage_inert(state, ticket, trial_s, &
          trial_t, step, inert_heat) bind(C, &
          name="fusion_c_source_state_stage_inert") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       integer(c_int64_t), value :: ticket
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       type(c_ptr), value :: step
       type(c_ptr), value :: inert_heat
       integer(c_int) :: status
     end function c_fusion_source_state_stage_inert

     function c_fusion_source_state_snapshot_inert(state, accepted_s, &
          accepted_t, cumulative, cumulative_inert_heat, accepted_time, &
          epoch) bind(C, name="fusion_c_source_state_snapshot_inert") &
          result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       type(c_ptr), value :: accepted_s
       type(c_ptr), value :: accepted_t
       type(c_ptr), value :: cumulative
       type(c_ptr), value :: cumulative_inert_heat
       type(c_ptr), value :: accepted_time
       type(c_ptr), value :: epoch
       integer(c_int) :: status
     end function c_fusion_source_state_snapshot_inert

     function c_fusion_source_state_commit(state, ticket) bind(C, &
          name="fusion_c_source_state_commit") result(status)
       import :: c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       integer(c_int64_t), value :: ticket
       integer(c_int) :: status
     end function c_fusion_source_state_commit

     function c_fusion_source_state_discard(state, ticket) bind(C, &
          name="fusion_c_source_state_discard") result(status)
       import :: c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       integer(c_int64_t), value :: ticket
       integer(c_int) :: status
     end function c_fusion_source_state_discard

     function c_fusion_source_state_pack_size(state, required) bind(C, &
          name="fusion_c_source_state_pack_size") result(status)
       import :: c_int, c_ptr, c_size_t
       type(c_ptr), value :: state
       type(c_ptr), value :: required
       integer(c_int) :: status
     end function c_fusion_source_state_pack_size

     function c_fusion_source_state_pack(state, buffer, capacity, written) &
          bind(C, name="fusion_c_source_state_pack") result(status)
       import :: c_int, c_ptr, c_size_t
       type(c_ptr), value :: state
       type(c_ptr), value :: buffer
       integer(c_size_t), value :: capacity
       type(c_ptr), value :: written
       integer(c_int) :: status
     end function c_fusion_source_state_pack

     function c_fusion_source_state_unpack(buffer, length, expected_model_tag, &
          out) bind(C, name="fusion_c_source_state_unpack") result(status)
       import :: c_int, c_int64_t, c_ptr, c_size_t
       type(c_ptr), value :: buffer
       integer(c_size_t), value :: length
       integer(c_int64_t), value :: expected_model_tag
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_source_state_unpack
  end interface

contains

  subroutine clear_ledger(ledger)
    type(fusion_source_ledger_v1), intent(out) :: ledger

    ledger%events_m3 = 0.0_c_double
    ledger%nuclear_born_number_m3 = 0.0_c_double
    ledger%nuclear_born_energy_J_m3 = 0.0_c_double
    ledger%external_born_number_m3 = 0.0_c_double
    ledger%external_born_energy_J_m3 = 0.0_c_double
    ledger%thermal_consumed_number_m3 = 0.0_c_double
    ledger%thermal_consumed_energy_J_m3 = 0.0_c_double
    ledger%fast_consumed_number_m3 = 0.0_c_double
    ledger%fast_consumed_energy_J_m3 = 0.0_c_double
    ledger%escaped_number_m3 = 0.0_c_double
    ledger%escaped_energy_J_m3 = 0.0_c_double
    ledger%handed_off_number_m3 = 0.0_c_double
    ledger%handed_off_energy_J_m3 = 0.0_c_double
    ledger%heat_to_bath_J_m3 = 0.0_c_double
    ledger%neutron_number_m3 = 0.0_c_double
    ledger%neutron_energy_J_m3 = 0.0_c_double
  end subroutine clear_ledger

  logical function valid_cells(cells)
    integer(c_int), intent(in) :: cells

    valid_cells = cells >= 1_c_int .and. &
         cells <= FUSION_SOURCE_STATE_MAX_CELLS
  end function valid_cells

  integer(c_size_t) function state_value_count(cells)
    integer(c_int), intent(in) :: cells

    state_value_count = int(cells, c_size_t) * &
         int(FUSION_SOURCE_STATE_SPECIES, c_size_t)
  end function state_value_count

  subroutine query_cells(state, cells, status)
    type(c_ptr), intent(in) :: state
    integer(c_int), intent(out), target :: cells
    integer(c_int), intent(out) :: status

    type(c_ptr) :: cells_ptr

    cells = 0_c_int
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return

    cells_ptr = c_loc(cells)
    status = c_fusion_source_state_cells(state, cells_ptr)
    if (status /= PB11_STATUS_OK) then
       cells = 0_c_int
       return
    end if
    if (.not. valid_cells(cells)) then
       cells = 0_c_int
       status = PB11_STATUS_INVALID_ARGUMENT
    end if
  end subroutine query_cells

  subroutine invalidate_staged_candidate(state, ticket)
    type(c_ptr), intent(in) :: state
    integer(c_int64_t), intent(in) :: ticket
    integer(c_int) :: ignored_status

    ! c_fusion_source_state_stage clears its staged flag before validating
    ! inputs.  Supplying three null pointers therefore invalidates a prior
    ! candidate without forming an address for malformed Fortran arrays.
    ignored_status = c_fusion_source_state_stage(state, ticket, c_null_ptr, &
         c_null_ptr, c_null_ptr)
  end subroutine invalidate_staged_candidate

  subroutine invalidate_staged_candidate_inert(state, ticket)
    type(c_ptr), intent(in) :: state
    integer(c_int64_t), intent(in) :: ticket
    integer(c_int) :: ignored_status

    ! The extended C stage clears its staged flag before checking pointers.
    ! Passing all null pointers invalidates an earlier candidate without ever
    ! taking C_LOC of malformed Fortran arrays.
    ignored_status = c_fusion_source_state_stage_inert(state, ticket, &
         c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr)
  end subroutine invalidate_staged_candidate_inert

  logical function state_arrays_match(cells, first_size, second_size)
    integer(c_int), intent(in) :: cells
    integer(c_size_t), intent(in) :: first_size, second_size
    integer(c_size_t) :: expected

    state_arrays_match = .false.
    if (.not. valid_cells(cells)) return
    expected = state_value_count(cells)
    state_arrays_match = first_size == expected .and. second_size == expected
  end function state_arrays_match

  subroutine fusion_source_state_create(cells, edges_J, initial_s_m3, &
       initial_t_m3, initial_time_s, model_tag, state, status)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: initial_s_m3(:)
    real(c_double), intent(in), target, contiguous :: initial_t_m3(:)
    real(c_double), intent(in) :: initial_time_s
    integer(c_int64_t), intent(in) :: model_tag
    type(c_ptr), intent(out), target :: state
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: expected_values
    type(c_ptr), target :: created
    integer(c_int) :: actual_cells, query_status

    state = c_null_ptr
    created = c_null_ptr
    status = PB11_STATUS_INVALID_ARGUMENT

    if (.not. valid_cells(cells)) return
    if (size(edges_J, kind=c_size_t) /= int(cells, c_size_t) + &
         1_c_size_t) return
    expected_values = state_value_count(cells)
    if (size(initial_s_m3, kind=c_size_t) /= expected_values) return
    if (size(initial_t_m3, kind=c_size_t) /= expected_values) return

    status = c_fusion_source_state_create(cells, c_loc(edges_J(1)), &
         c_loc(initial_s_m3(1)), c_loc(initial_t_m3(1)), initial_time_s, &
         model_tag, c_loc(created))
    if (status /= PB11_STATUS_OK) then
       state = c_null_ptr
       return
    end if

    ! Verify the immutable C-side dimension before publishing the handle.
    call query_cells(created, actual_cells, query_status)
    if (query_status /= PB11_STATUS_OK .or. actual_cells /= cells) then
       call c_fusion_source_state_destroy(created)
       state = c_null_ptr
       if (query_status == PB11_STATUS_OK) then
          status = PB11_STATUS_INVALID_ARGUMENT
       else
          status = query_status
       end if
       return
    end if
    state = created
  end subroutine fusion_source_state_create

  subroutine fusion_source_state_destroy(state)
    type(c_ptr), intent(inout) :: state

    if (c_associated(state)) call c_fusion_source_state_destroy(state)
    state = c_null_ptr
  end subroutine fusion_source_state_destroy

  subroutine fusion_source_state_cells(state, cells, status)
    type(c_ptr), intent(in) :: state
    integer(c_int), intent(out), target :: cells
    integer(c_int), intent(out) :: status

    call query_cells(state, cells, status)
  end subroutine fusion_source_state_cells

  subroutine fusion_source_state_snapshot(state, accepted_s_m3, &
       accepted_t_m3, cumulative, accepted_time_s, epoch, status)
    type(c_ptr), intent(in) :: state
    real(c_double), intent(out), target, contiguous :: accepted_s_m3(:)
    real(c_double), intent(out), target, contiguous :: accepted_t_m3(:)
    type(fusion_source_ledger_v1), intent(out), target :: cumulative
    real(c_double), intent(out), target :: accepted_time_s
    integer(c_int64_t), intent(out), target :: epoch
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, query_status
    type(c_ptr) :: state_ptr

    accepted_s_m3 = 0.0_c_double
    accepted_t_m3 = 0.0_c_double
    call clear_ledger(cumulative)
    accepted_time_s = 0.0_c_double
    epoch = 0_c_int64_t
    status = PB11_STATUS_INVALID_ARGUMENT

    call query_cells(state, cells, query_status)
    if (query_status /= PB11_STATUS_OK) then
       status = query_status
       return
    end if
    if (.not. state_arrays_match(cells, size(accepted_s_m3, kind=c_size_t), &
         size(accepted_t_m3, kind=c_size_t))) return

    state_ptr = state
    status = c_fusion_source_state_snapshot(state_ptr, &
         c_loc(accepted_s_m3(1)), c_loc(accepted_t_m3(1)), &
         c_loc(cumulative), c_loc(accepted_time_s), c_loc(epoch))
    if (status /= PB11_STATUS_OK) then
       accepted_s_m3 = 0.0_c_double
       accepted_t_m3 = 0.0_c_double
       call clear_ledger(cumulative)
       accepted_time_s = 0.0_c_double
       epoch = 0_c_int64_t
    end if
  end subroutine fusion_source_state_snapshot

  subroutine fusion_source_state_snapshot_inert(state, accepted_s_m3, &
       accepted_t_m3, cumulative, cumulative_inert_heat_J_m3, &
       accepted_time_s, epoch, status)
    type(c_ptr), intent(in) :: state
    real(c_double), intent(out), target, contiguous :: accepted_s_m3(:)
    real(c_double), intent(out), target, contiguous :: accepted_t_m3(:)
    type(fusion_source_ledger_v1), intent(out), target :: cumulative
    real(c_double), intent(out), target, contiguous :: &
         cumulative_inert_heat_J_m3(:)
    real(c_double), intent(out), target :: accepted_time_s
    integer(c_int64_t), intent(out), target :: epoch
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, query_status
    integer(c_size_t) :: inert_count
    type(c_ptr) :: state_ptr

    ! Match the original Fortran snapshot policy: clear every output before
    ! validation and again if the C call rejects the request.  This is stricter
    ! than the C snapshot_inert contract, which leaves kinetic arrays alone on
    ! failure, and gives Fortran callers deterministic output state.
    accepted_s_m3 = 0.0_c_double
    accepted_t_m3 = 0.0_c_double
    call clear_ledger(cumulative)
    cumulative_inert_heat_J_m3 = 0.0_c_double
    accepted_time_s = 0.0_c_double
    epoch = 0_c_int64_t
    status = PB11_STATUS_INVALID_ARGUMENT

    call query_cells(state, cells, query_status)
    if (query_status /= PB11_STATUS_OK) then
       status = query_status
       return
    end if
    inert_count = int(FUSION_SOURCE_STATE_SPECIES, c_size_t)
    if (.not. state_arrays_match(cells, size(accepted_s_m3, kind=c_size_t), &
         size(accepted_t_m3, kind=c_size_t))) return
    if (size(cumulative_inert_heat_J_m3, kind=c_size_t) /= inert_count) return

    state_ptr = state
    status = c_fusion_source_state_snapshot_inert(state_ptr, &
         c_loc(accepted_s_m3(1)), c_loc(accepted_t_m3(1)), &
         c_loc(cumulative), c_loc(cumulative_inert_heat_J_m3(1)), &
         c_loc(accepted_time_s), c_loc(epoch))
    if (status /= PB11_STATUS_OK) then
       accepted_s_m3 = 0.0_c_double
       accepted_t_m3 = 0.0_c_double
       call clear_ledger(cumulative)
       cumulative_inert_heat_J_m3 = 0.0_c_double
       accepted_time_s = 0.0_c_double
       epoch = 0_c_int64_t
    end if
  end subroutine fusion_source_state_snapshot_inert

  subroutine fusion_source_state_begin(state, dt_s, ticket, status)
    type(c_ptr), intent(in) :: state
    real(c_double), intent(in) :: dt_s
    integer(c_int64_t), intent(out), target :: ticket
    integer(c_int), intent(out) :: status

    ticket = 0_c_int64_t
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return

    status = c_fusion_source_state_begin(state, dt_s, c_loc(ticket))
    if (status /= PB11_STATUS_OK) ticket = 0_c_int64_t
  end subroutine fusion_source_state_begin

  subroutine fusion_source_state_stage(state, ticket, trial_s_m3, &
       trial_t_m3, step, status)
    type(c_ptr), intent(in) :: state
    integer(c_int64_t), intent(in) :: ticket
    real(c_double), intent(in), target, contiguous :: trial_s_m3(:)
    real(c_double), intent(in), target, contiguous :: trial_t_m3(:)
    type(fusion_source_ledger_v1), intent(in), target :: step
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, query_status
    type(c_ptr) :: state_ptr

    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return

    call query_cells(state, cells, query_status)
    if (query_status /= PB11_STATUS_OK) then
       ! A failed C-side stage invalidates a previously staged candidate.  Do
       ! the same for a query failure while the opaque handle is still known.
       call invalidate_staged_candidate(state, ticket)
       status = query_status
       return
    end if
    if (.not. state_arrays_match(cells, size(trial_s_m3, kind=c_size_t), &
         size(trial_t_m3, kind=c_size_t))) then
       ! Do not form C_LOC for malformed extents.  A null stage is the C ABI's
       ! validation path and clears any earlier stage for this ticket while
       ! retaining the ticket for a corrected re-stage.
       call invalidate_staged_candidate(state, ticket)
       status = PB11_STATUS_INVALID_ARGUMENT
       return
    end if

    state_ptr = state
    status = c_fusion_source_state_stage(state_ptr, ticket, &
         c_loc(trial_s_m3(1)), c_loc(trial_t_m3(1)), c_loc(step))
  end subroutine fusion_source_state_stage

  subroutine fusion_source_state_stage_inert(state, ticket, trial_s_m3, &
       trial_t_m3, step, inert_heat_J_m3, status)
    type(c_ptr), intent(in) :: state
    integer(c_int64_t), intent(in) :: ticket
    real(c_double), intent(in), target, contiguous :: trial_s_m3(:)
    real(c_double), intent(in), target, contiguous :: trial_t_m3(:)
    type(fusion_source_ledger_v1), intent(in), target :: step
    real(c_double), intent(in), target, contiguous :: inert_heat_J_m3(:)
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, query_status
    integer(c_size_t) :: inert_count
    type(c_ptr) :: state_ptr

    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return

    call query_cells(state, cells, query_status)
    if (query_status /= PB11_STATUS_OK) then
       call invalidate_staged_candidate_inert(state, ticket)
       status = query_status
       return
    end if
    if (.not. state_arrays_match(cells, size(trial_s_m3, kind=c_size_t), &
         size(trial_t_m3, kind=c_size_t))) then
       call invalidate_staged_candidate_inert(state, ticket)
       status = PB11_STATUS_INVALID_ARGUMENT
       return
    end if
    inert_count = int(FUSION_SOURCE_STATE_SPECIES, c_size_t)
    if (size(inert_heat_J_m3, kind=c_size_t) /= inert_count) then
       call invalidate_staged_candidate_inert(state, ticket)
       status = PB11_STATUS_INVALID_ARGUMENT
       return
    end if

    state_ptr = state
    status = c_fusion_source_state_stage_inert(state_ptr, ticket, &
         c_loc(trial_s_m3(1)), c_loc(trial_t_m3(1)), c_loc(step), &
         c_loc(inert_heat_J_m3(1)))
  end subroutine fusion_source_state_stage_inert

  subroutine fusion_source_state_commit(state, ticket, status)
    type(c_ptr), intent(in) :: state
    integer(c_int64_t), intent(in) :: ticket
    integer(c_int), intent(out) :: status

    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return
    status = c_fusion_source_state_commit(state, ticket)
  end subroutine fusion_source_state_commit

  subroutine fusion_source_state_discard(state, ticket, status)
    type(c_ptr), intent(in) :: state
    integer(c_int64_t), intent(in) :: ticket
    integer(c_int), intent(out) :: status

    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return
    status = c_fusion_source_state_discard(state, ticket)
  end subroutine fusion_source_state_discard

  subroutine fusion_source_state_pack_size(state, required, status)
    type(c_ptr), intent(in) :: state
    integer(c_size_t), intent(out), target :: required
    integer(c_int), intent(out) :: status

    required = 0_c_size_t
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return

    status = c_fusion_source_state_pack_size(state, c_loc(required))
    if (status /= PB11_STATUS_OK) required = 0_c_size_t
  end subroutine fusion_source_state_pack_size

  subroutine fusion_source_state_pack(state, buffer, capacity, &
       bytes_written, status)
    type(c_ptr), intent(in) :: state
    integer(c_int8_t), intent(out), target, contiguous :: buffer(:)
    integer(c_size_t), intent(in) :: capacity
    integer(c_size_t), intent(out), target :: bytes_written
    integer(c_int), intent(out) :: status

    bytes_written = 0_c_size_t
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return
    ! The C routine requires a nonnull output buffer.  Check capacity and the
    ! actual Fortran extent before referring to buffer(1), including for a
    ! valid zero-length Fortran actual.
    if (capacity <= 0_c_size_t) return
    if (size(buffer, kind=c_size_t) < capacity) return

    status = c_fusion_source_state_pack(state, c_loc(buffer(1)), capacity, &
         c_loc(bytes_written))
    if (status /= PB11_STATUS_OK) bytes_written = 0_c_size_t
  end subroutine fusion_source_state_pack

  subroutine fusion_source_state_unpack(buffer, length, expected_model_tag, &
       state, status)
    integer(c_int8_t), intent(in), target, contiguous :: buffer(:)
    integer(c_size_t), intent(in) :: length
    integer(c_int64_t), intent(in) :: expected_model_tag
    type(c_ptr), intent(out), target :: state
    integer(c_int), intent(out) :: status

    type(c_ptr), target :: unpacked
    integer(c_int) :: cells, query_status

    state = c_null_ptr
    unpacked = c_null_ptr
    status = PB11_STATUS_INVALID_ARGUMENT
    ! Do not form C_LOC for a zero-length actual or for a requested length
    ! beyond the actual buffer.  The C ABI receives exactly the requested
    ! prefix when a larger Fortran buffer is supplied.
    if (length <= 0_c_size_t) return
    if (size(buffer, kind=c_size_t) < length) return

    status = c_fusion_source_state_unpack(c_loc(buffer(1)), length, &
         expected_model_tag, c_loc(unpacked))
    if (status /= PB11_STATUS_OK) then
       state = c_null_ptr
       return
    end if

    call query_cells(unpacked, cells, query_status)
    if (query_status /= PB11_STATUS_OK) then
       call c_fusion_source_state_destroy(unpacked)
       state = c_null_ptr
       status = query_status
       return
    end if
    state = unpacked
  end subroutine fusion_source_state_unpack

end module fusion_source_state_fortran
