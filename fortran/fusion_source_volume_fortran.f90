module fusion_source_volume_fortran
  !! ISO_C_BINDING wrappers for the volume-aware source-state context.
  !!
  !! The opaque context remains a type(c_ptr) owned by the caller.  The base
  !! source-state module owns destroy/begin/commit/discard/restart operations;
  !! this module adds only the volume-aware create, stage, and snapshot calls.
  !! Densities use the Fortran shape (cells,6).  Since the first dimension is
  !! contiguous, that shape has the C species-major [species][cell] layout.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t, &
       c_size_t, c_ptr, c_null_ptr, c_loc, c_associated
  use fusion_source_state_fortran, only : fusion_source_ledger_v1, &
       FUSION_SOURCE_STATE_SPECIES, FUSION_SOURCE_STATE_MAX_CELLS, &
       FUSION_SOURCE_STATE_LEDGER_DOUBLES, PB11_STATUS_OK, &
       PB11_STATUS_NULL_OUTPUT, PB11_STATUS_INVALID_ARGUMENT, &
       PB11_STATUS_OUT_OF_RANGE, PB11_STATUS_NUMERICAL_FAILURE, &
       PB11_STATUS_EXCEPTION, PB11_STATUS_UNKNOWN_METHOD, &
       fusion_source_state_cells, fusion_source_state_destroy
  implicit none
  private

  public :: PB11_STATUS_OK
  public :: PB11_STATUS_NULL_OUTPUT
  public :: PB11_STATUS_INVALID_ARGUMENT
  public :: PB11_STATUS_OUT_OF_RANGE
  public :: PB11_STATUS_NUMERICAL_FAILURE
  public :: PB11_STATUS_EXCEPTION
  public :: PB11_STATUS_UNKNOWN_METHOD
  public :: FUSION_SOURCE_STATE_SPECIES
  public :: FUSION_SOURCE_STATE_MAX_CELLS
  public :: FUSION_SOURCE_STATE_LEDGER_DOUBLES
  public :: fusion_source_ledger_v1
  public :: fusion_source_state_create_volume
  public :: fusion_source_state_stage_volume
  public :: fusion_source_state_snapshot_volume

  ! Exact C layout: seven consecutive six-element c_double arrays.
  type, bind(C), public :: fusion_transport_ledger_v1
     real(c_double) :: spatial_number(6)
     real(c_double) :: spatial_energy_J(6)
     real(c_double) :: work_J(6)
     real(c_double) :: lower_number(6)
     real(c_double) :: lower_energy_J(6)
     real(c_double) :: upper_number(6)
     real(c_double) :: upper_energy_J(6)
  end type fusion_transport_ledger_v1

  interface
     function c_fusion_source_state_create_volume(cells, edges, initial_s, &
          initial_t, initial_volume, initial_time, model_tag, out) bind(C, &
          name="fusion_c_source_state_create_volume") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: initial_s
       type(c_ptr), value :: initial_t
       real(c_double), value :: initial_volume
       real(c_double), value :: initial_time
       integer(c_int64_t), value :: model_tag
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_source_state_create_volume

     function c_fusion_source_state_stage_volume(state, ticket, trial_s, &
          trial_t, trial_volume, source_volume, source, inert, transport) &
          bind(C, name="fusion_c_source_state_stage_volume") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       integer(c_int64_t), value :: ticket
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       real(c_double), value :: trial_volume
       real(c_double), value :: source_volume
       type(c_ptr), value :: source
       type(c_ptr), value :: inert
       type(c_ptr), value :: transport
       integer(c_int) :: status
     end function c_fusion_source_state_stage_volume

     function c_fusion_source_state_snapshot_volume(state, accepted_s, &
          accepted_t, cumulative_reference, inert_reference, &
          cumulative_transport, reference_volume, accepted_volume, &
          accepted_time, epoch) bind(C, &
          name="fusion_c_source_state_snapshot_volume") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       type(c_ptr), value :: accepted_s
       type(c_ptr), value :: accepted_t
       type(c_ptr), value :: cumulative_reference
       type(c_ptr), value :: inert_reference
       type(c_ptr), value :: cumulative_transport
       type(c_ptr), value :: reference_volume
       type(c_ptr), value :: accepted_volume
       type(c_ptr), value :: accepted_time
       type(c_ptr), value :: epoch
       integer(c_int) :: status
     end function c_fusion_source_state_snapshot_volume
  end interface

contains

  subroutine clear_source_ledger(ledger)
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
  end subroutine clear_source_ledger

  subroutine clear_transport_ledger(ledger)
    type(fusion_transport_ledger_v1), intent(out) :: ledger

    ledger%spatial_number = 0.0_c_double
    ledger%spatial_energy_J = 0.0_c_double
    ledger%work_J = 0.0_c_double
    ledger%lower_number = 0.0_c_double
    ledger%lower_energy_J = 0.0_c_double
    ledger%upper_number = 0.0_c_double
    ledger%upper_energy_J = 0.0_c_double
  end subroutine clear_transport_ledger

  logical function valid_cells(cells)
    integer(c_int), intent(in) :: cells

    valid_cells = cells >= 1_c_int .and. &
         cells <= FUSION_SOURCE_STATE_MAX_CELLS
  end function valid_cells

  logical function state_shapes_match(cells, first_rows, first_cols, &
       second_rows, second_cols)
    integer(c_int), intent(in) :: cells
    integer(c_size_t), intent(in) :: first_rows, first_cols
    integer(c_size_t), intent(in) :: second_rows, second_cols
    integer(c_size_t) :: expected_rows, expected_cols

    state_shapes_match = .false.
    if (.not. valid_cells(cells)) return
    expected_rows = int(cells, c_size_t)
    expected_cols = int(FUSION_SOURCE_STATE_SPECIES, c_size_t)
    state_shapes_match = first_rows == expected_rows .and. &
         first_cols == expected_cols .and. second_rows == expected_rows .and. &
         second_cols == expected_cols
  end function state_shapes_match

  subroutine query_cells(state, cells, status)
    type(c_ptr), intent(in) :: state
    integer(c_int), intent(out), target :: cells
    integer(c_int), intent(out) :: status

    cells = 0_c_int
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return

    call fusion_source_state_cells(state, cells, status)
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

    ! The C stage clears its staged flag before validating pointers.  Passing
    ! null pointers therefore invalidates a stale candidate without forming a
    ! C address for malformed Fortran arrays.  The returned status is ignored.
    ignored_status = c_fusion_source_state_stage_volume(state, ticket, &
         c_null_ptr, c_null_ptr, 0.0_c_double, 0.0_c_double, c_null_ptr, &
         c_null_ptr, c_null_ptr)
  end subroutine invalidate_staged_candidate

  subroutine fusion_source_state_create_volume(cells, edges_J, initial_s_m3, &
       initial_t_m3, initial_volume_m3, initial_time_s, model_tag, state, &
       status)
    integer(c_int), intent(in) :: cells
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: initial_s_m3(:,:)
    real(c_double), intent(in), target, contiguous :: initial_t_m3(:,:)
    real(c_double), intent(in) :: initial_volume_m3
    real(c_double), intent(in) :: initial_time_s
    integer(c_int64_t), intent(in) :: model_tag
    type(c_ptr), intent(inout), target :: state
    integer(c_int), intent(out) :: status

    type(c_ptr), target :: created
    integer(c_int) :: actual_cells, query_status
    integer(c_size_t) :: expected_rows, expected_cols

    status = PB11_STATUS_INVALID_ARGUMENT
    ! A live handle belongs to the caller.  Requiring an initialized c_ptr and
    ! rejecting it here prevents a successful create from leaking/overwriting
    ! an existing context.
    if (c_associated(state)) return
    state = c_null_ptr
    created = c_null_ptr

    if (.not. valid_cells(cells)) return
    expected_rows = int(cells, c_size_t)
    expected_cols = int(FUSION_SOURCE_STATE_SPECIES, c_size_t)
    if (size(edges_J, kind=c_size_t) /= expected_rows + 1_c_size_t) return
    if (size(initial_s_m3, 1, kind=c_size_t) /= expected_rows) return
    if (size(initial_s_m3, 2, kind=c_size_t) /= expected_cols) return
    if (size(initial_t_m3, 1, kind=c_size_t) /= expected_rows) return
    if (size(initial_t_m3, 2, kind=c_size_t) /= expected_cols) return

    status = c_fusion_source_state_create_volume(cells, c_loc(edges_J(1)), &
         c_loc(initial_s_m3(1,1)), c_loc(initial_t_m3(1,1)), &
         initial_volume_m3, initial_time_s, model_tag, c_loc(created))
    if (status /= PB11_STATUS_OK) then
       state = c_null_ptr
       return
    end if

    ! Verify the immutable C-side dimension before publishing ownership to the
    ! caller, matching the base source-state create wrapper.
    call query_cells(created, actual_cells, query_status)
    if (query_status /= PB11_STATUS_OK .or. actual_cells /= cells) then
       call fusion_source_state_destroy(created)
       state = c_null_ptr
       if (query_status == PB11_STATUS_OK) then
          status = PB11_STATUS_INVALID_ARGUMENT
       else
          status = query_status
       end if
       return
    end if
    state = created
  end subroutine fusion_source_state_create_volume

  subroutine fusion_source_state_stage_volume(state, ticket, trial_s_m3, &
       trial_t_m3, trial_volume_m3, source_volume_m3, source_step, &
       inert_heat_J_m3, transport_step, status)
    type(c_ptr), intent(in) :: state
    integer(c_int64_t), intent(in) :: ticket
    real(c_double), intent(in), target, contiguous :: trial_s_m3(:,:)
    real(c_double), intent(in), target, contiguous :: trial_t_m3(:,:)
    real(c_double), intent(in) :: trial_volume_m3
    real(c_double), intent(in) :: source_volume_m3
    type(fusion_source_ledger_v1), intent(in), target :: source_step
    real(c_double), intent(in), target, contiguous :: inert_heat_J_m3(:)
    type(fusion_transport_ledger_v1), intent(in), target :: transport_step
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, query_status
    integer(c_size_t) :: inert_count
    type(c_ptr) :: state_ptr

    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(state)) return

    call query_cells(state, cells, query_status)
    if (query_status /= PB11_STATUS_OK) then
       call invalidate_staged_candidate(state, ticket)
       status = query_status
       return
    end if
    if (.not. state_shapes_match(cells, &
         size(trial_s_m3, 1, kind=c_size_t), &
         size(trial_s_m3, 2, kind=c_size_t), &
         size(trial_t_m3, 1, kind=c_size_t), &
         size(trial_t_m3, 2, kind=c_size_t))) then
       call invalidate_staged_candidate(state, ticket)
       return
    end if
    inert_count = int(FUSION_SOURCE_STATE_SPECIES, c_size_t)
    if (size(inert_heat_J_m3, kind=c_size_t) /= inert_count) then
       call invalidate_staged_candidate(state, ticket)
       return
    end if

    state_ptr = state
    status = c_fusion_source_state_stage_volume(state_ptr, ticket, &
         c_loc(trial_s_m3(1,1)), c_loc(trial_t_m3(1,1)), trial_volume_m3, &
         source_volume_m3, c_loc(source_step), c_loc(inert_heat_J_m3(1)), &
         c_loc(transport_step))
    if (status /= PB11_STATUS_OK) then
       ! C clears its staged flag on all validation failures.  Repeat the
       ! null-pointer invalidation so this wrapper keeps that guarantee if a
       ! future C implementation rejects before entering its normal stage path.
       call invalidate_staged_candidate(state, ticket)
    end if
  end subroutine fusion_source_state_stage_volume

  subroutine fusion_source_state_snapshot_volume(state, accepted_s_m3, &
       accepted_t_m3, cumulative_reference, inert_heat_reference_J_m3, &
       cumulative_transport, reference_volume_m3, accepted_volume_m3, &
       accepted_time_s, epoch, status)
    type(c_ptr), intent(in) :: state
    real(c_double), intent(out), target, contiguous :: accepted_s_m3(:,:)
    real(c_double), intent(out), target, contiguous :: accepted_t_m3(:,:)
    type(fusion_source_ledger_v1), intent(out), target :: cumulative_reference
    real(c_double), intent(out), target, contiguous :: &
         inert_heat_reference_J_m3(:)
    type(fusion_transport_ledger_v1), intent(out), target :: &
         cumulative_transport
    real(c_double), intent(out), target :: reference_volume_m3
    real(c_double), intent(out), target :: accepted_volume_m3
    real(c_double), intent(out), target :: accepted_time_s
    integer(c_int64_t), intent(out), target :: epoch
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, query_status
    integer(c_size_t) :: inert_count
    type(c_ptr) :: state_ptr

    ! Clear every output before validation.  This also covers null handles,
    ! malformed extents, and a C-side failure after the address checks.
    accepted_s_m3 = 0.0_c_double
    accepted_t_m3 = 0.0_c_double
    call clear_source_ledger(cumulative_reference)
    inert_heat_reference_J_m3 = 0.0_c_double
    call clear_transport_ledger(cumulative_transport)
    reference_volume_m3 = 0.0_c_double
    accepted_volume_m3 = 0.0_c_double
    accepted_time_s = 0.0_c_double
    epoch = 0_c_int64_t
    status = PB11_STATUS_INVALID_ARGUMENT

    call query_cells(state, cells, query_status)
    if (query_status /= PB11_STATUS_OK) then
       status = query_status
       return
    end if
    if (.not. state_shapes_match(cells, &
         size(accepted_s_m3, 1, kind=c_size_t), &
         size(accepted_s_m3, 2, kind=c_size_t), &
         size(accepted_t_m3, 1, kind=c_size_t), &
         size(accepted_t_m3, 2, kind=c_size_t))) return
    inert_count = int(FUSION_SOURCE_STATE_SPECIES, c_size_t)
    if (size(inert_heat_reference_J_m3, kind=c_size_t) /= inert_count) return

    state_ptr = state
    status = c_fusion_source_state_snapshot_volume(state_ptr, &
         c_loc(accepted_s_m3(1,1)), c_loc(accepted_t_m3(1,1)), &
         c_loc(cumulative_reference), c_loc(inert_heat_reference_J_m3(1)), &
         c_loc(cumulative_transport), c_loc(reference_volume_m3), &
         c_loc(accepted_volume_m3), c_loc(accepted_time_s), c_loc(epoch))
    if (status /= PB11_STATUS_OK) then
       accepted_s_m3 = 0.0_c_double
       accepted_t_m3 = 0.0_c_double
       call clear_source_ledger(cumulative_reference)
       inert_heat_reference_J_m3 = 0.0_c_double
       call clear_transport_ledger(cumulative_transport)
       reference_volume_m3 = 0.0_c_double
       accepted_volume_m3 = 0.0_c_double
       accepted_time_s = 0.0_c_double
       epoch = 0_c_int64_t
    end if
  end subroutine fusion_source_state_snapshot_volume

end module fusion_source_volume_fortran
