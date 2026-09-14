program test_fusion_source_volume_fortran
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_int8_t, c_int64_t, c_null_ptr, c_ptr, c_size_t, c_sizeof
  use fusion_source_state_fortran, only : PB11_STATUS_OK, &
       PB11_STATUS_INVALID_ARGUMENT, fusion_source_state_begin, &
       fusion_source_state_commit, fusion_source_state_destroy, &
       fusion_source_state_pack, fusion_source_state_pack_size, &
       fusion_source_state_unpack
  use fusion_source_volume_fortran, only : fusion_source_ledger_v1, &
       fusion_transport_ledger_v1, fusion_source_state_create_volume, &
       fusion_source_state_stage_volume, fusion_source_state_snapshot_volume
  implicit none

  integer(c_int), parameter :: test_cells = 2_c_int
  integer(c_int64_t), parameter :: test_tag = 9988776655443322_c_int64_t
  integer :: failures

  failures = 0
  call test_layout_and_live_handle()
  call test_volume_dilution_and_source_volume()
  call test_bad_stage_invalidates_candidate()
  call test_bad_snapshot_clears_outputs()
  call test_volume_v3_restart()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, &
          ' Fortran source-volume binding test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran source-volume binding tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  subroutine check_status(status, message)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: message

    call check(status == PB11_STATUS_OK, message)
  end subroutine check_status

  subroutine zero_source(ledger)
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
  end subroutine zero_source

  subroutine zero_transport(ledger)
    type(fusion_transport_ledger_v1), intent(out) :: ledger

    ledger%spatial_number = 0.0_c_double
    ledger%spatial_energy_J = 0.0_c_double
    ledger%work_J = 0.0_c_double
    ledger%lower_number = 0.0_c_double
    ledger%lower_energy_J = 0.0_c_double
    ledger%upper_number = 0.0_c_double
    ledger%upper_energy_J = 0.0_c_double
  end subroutine zero_transport

  logical function source_is_zero(ledger)
    type(fusion_source_ledger_v1), intent(in) :: ledger

    source_is_zero = all(ledger%events_m3 == 0.0_c_double) .and. &
         all(ledger%nuclear_born_number_m3 == 0.0_c_double) .and. &
         all(ledger%nuclear_born_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%external_born_number_m3 == 0.0_c_double) .and. &
         all(ledger%external_born_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%thermal_consumed_number_m3 == 0.0_c_double) .and. &
         all(ledger%thermal_consumed_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%fast_consumed_number_m3 == 0.0_c_double) .and. &
         all(ledger%fast_consumed_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%escaped_number_m3 == 0.0_c_double) .and. &
         all(ledger%escaped_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%handed_off_number_m3 == 0.0_c_double) .and. &
         all(ledger%handed_off_energy_J_m3 == 0.0_c_double) .and. &
         all(ledger%heat_to_bath_J_m3 == 0.0_c_double) .and. &
         ledger%neutron_number_m3 == 0.0_c_double .and. &
         ledger%neutron_energy_J_m3 == 0.0_c_double
  end function source_is_zero

  logical function transport_is_zero(ledger)
    type(fusion_transport_ledger_v1), intent(in) :: ledger

    transport_is_zero = all(ledger%spatial_number == 0.0_c_double) .and. &
         all(ledger%spatial_energy_J == 0.0_c_double) .and. &
         all(ledger%work_J == 0.0_c_double) .and. &
         all(ledger%lower_number == 0.0_c_double) .and. &
         all(ledger%lower_energy_J == 0.0_c_double) .and. &
         all(ledger%upper_number == 0.0_c_double) .and. &
         all(ledger%upper_energy_J == 0.0_c_double)
  end function transport_is_zero

  subroutine make_initial(edges, initial_s, initial_t)
    real(c_double), intent(out), target :: edges(:)
    real(c_double), intent(out), target :: initial_s(:,:)
    real(c_double), intent(out), target :: initial_t(:,:)

    edges = [0.0_c_double, 1.0_c_double, 2.0_c_double]
    initial_s = 0.0_c_double
    initial_t = 0.0_c_double
    initial_s(1,1) = 4.0_c_double
  end subroutine make_initial

  subroutine test_layout_and_live_handle()
    real(c_double), target :: edges(3), initial_s(2,6), initial_t(2,6)
    type(c_ptr) :: state
    integer(c_int) :: status
    type(fusion_transport_ledger_v1), target :: transport

    call make_initial(edges, initial_s, initial_t)
    call check(c_sizeof(transport) == 42_c_size_t * 8_c_size_t, &
         'transport ledger layout is exactly 42 C doubles')
    state = c_null_ptr
    call fusion_source_state_create_volume(test_cells, edges, initial_s, &
         initial_t, 2.0_c_double, 0.0_c_double, test_tag, state, status)
    call check_status(status, 'create volume context')
    call check(c_associated(state), 'create publishes a context')

    call fusion_source_state_create_volume(test_cells, edges, initial_s, &
         initial_t, 2.0_c_double, 0.0_c_double, test_tag, state, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'create refuses to overwrite a live context')
    call check(c_associated(state), 'live context remains owned after rejection')
    call fusion_source_state_destroy(state)
  end subroutine test_layout_and_live_handle

  subroutine test_volume_dilution_and_source_volume()
    real(c_double), target :: edges(3), initial_s(2,6), initial_t(2,6)
    real(c_double), target :: candidate_s(2,6), candidate_t(2,6)
    real(c_double), target :: inert(6)
    type(fusion_source_ledger_v1), target :: source, cumulative
    type(fusion_transport_ledger_v1), target :: transport, cumulative_transport
    real(c_double), target :: accepted_s(2,6), accepted_t(2,6)
    real(c_double), target :: reference_volume, accepted_volume, accepted_time
    integer(c_int64_t), target :: ticket, epoch
    integer(c_int) :: status
    type(c_ptr) :: state

    call make_initial(edges, initial_s, initial_t)
    candidate_t = 0.0_c_double
    inert = 0.0_c_double
    call zero_source(source)
    call zero_transport(transport)
    state = c_null_ptr
    call fusion_source_state_create_volume(test_cells, edges, initial_s, &
         initial_t, 2.0_c_double, 0.0_c_double, test_tag + 1_c_int64_t, &
         state, status)
    call check_status(status, 'create dilution context')

    candidate_s = initial_s
    candidate_s(1,1) = 8.0_c_double
    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin volume dilution')
    call fusion_source_state_stage_volume(state, ticket, candidate_s, &
         candidate_t, 1.0_c_double, 1.0_c_double, source, inert, transport, &
         status)
    call check_status(status, 'stage V=2 to V=1 dilution')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit V=2 to V=1 dilution')

    call fusion_source_state_snapshot_volume(state, accepted_s, accepted_t, &
         cumulative, inert, cumulative_transport, reference_volume, &
         accepted_volume, accepted_time, epoch, status)
    call check_status(status, 'snapshot diluted volume')
    call check(accepted_s(1,1) == 8.0_c_double .and. &
         reference_volume == 2.0_c_double .and. accepted_volume == 1.0_c_double, &
         'physical V=2 to V=1 dilution doubles the density')

    ! The source density is integrated over source_volume.  With a reference
    ! volume of 2 and source_volume of 4, 0.5 m^-3 contributes one reference
    ! amount. Bin 1 has midpoint 0.5 J, hence the energy density is
    ! 0.25 J/m3 at the source volume, or 0.5 J/m3 at the reference volume.
    call zero_source(source)
    source%external_born_number_m3(1) = 0.5_c_double
    source%external_born_energy_J_m3(1) = 0.25_c_double
    candidate_s(1,1) = 10.0_c_double
    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin explicit source-volume step')
    call fusion_source_state_stage_volume(state, ticket, candidate_s, &
         candidate_t, 1.0_c_double, 4.0_c_double, source, inert, transport, &
         status)
    call check_status(status, 'stage explicit source-volume normalization')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit explicit source-volume normalization')
    call fusion_source_state_snapshot_volume(state, accepted_s, accepted_t, &
         cumulative, inert, cumulative_transport, reference_volume, &
         accepted_volume, accepted_time, epoch, status)
    call check_status(status, 'snapshot source-volume step')
    call check(accepted_s(1,1) == 10.0_c_double .and. &
         cumulative%external_born_number_m3(1) == 1.0_c_double .and. &
         cumulative%external_born_energy_J_m3(1) == 0.5_c_double, &
         'source_volume scales source density into reference ledger units')
    call fusion_source_state_destroy(state)
  end subroutine test_volume_dilution_and_source_volume

  subroutine test_bad_stage_invalidates_candidate()
    real(c_double), target :: edges(3), initial_s(2,6), initial_t(2,6)
    real(c_double), target :: candidate_s(2,6), candidate_t(2,6)
    real(c_double), target :: bad_s(1,6), inert(6)
    type(fusion_source_ledger_v1), target :: source
    type(fusion_transport_ledger_v1), target :: transport
    integer(c_int64_t), target :: ticket
    integer(c_int) :: status
    type(c_ptr) :: state

    call make_initial(edges, initial_s, initial_t)
    candidate_s = initial_s
    candidate_s(1,1) = 8.0_c_double
    candidate_t = 0.0_c_double
    bad_s = 0.0_c_double
    inert = 0.0_c_double
    call zero_source(source)
    call zero_transport(transport)
    state = c_null_ptr
    call fusion_source_state_create_volume(test_cells, edges, initial_s, &
         initial_t, 2.0_c_double, 0.0_c_double, test_tag + 2_c_int64_t, &
         state, status)
    call check_status(status, 'create invalid-stage context')
    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin invalid-stage context')
    call fusion_source_state_stage_volume(state, ticket, candidate_s, &
         candidate_t, 1.0_c_double, 1.0_c_double, source, inert, transport, &
         status)
    call check_status(status, 'stage valid candidate before bad shape')

    call fusion_source_state_stage_volume(state, ticket, bad_s, candidate_t, &
         1.0_c_double, 1.0_c_double, source, inert, transport, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'stage rejects wrong cell extent')
    call fusion_source_state_commit(state, ticket, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong extent invalidates the previous staged candidate')

    call fusion_source_state_stage_volume(state, ticket, candidate_s, &
         candidate_t, 1.0_c_double, 1.0_c_double, source, inert, transport, &
         status)
    call check_status(status, 'restage after invalid extent')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit corrected restage')
    call fusion_source_state_destroy(state)
  end subroutine test_bad_stage_invalidates_candidate

  subroutine test_bad_snapshot_clears_outputs()
    real(c_double), target :: edges(3), initial_s(2,6), initial_t(2,6)
    real(c_double), target :: candidate_s(2,6), candidate_t(2,6)
    real(c_double), target :: bad_accepted_s(1,6), accepted_t(2,6)
    real(c_double), target :: inert(6)
    type(fusion_source_ledger_v1), target :: source, cumulative
    type(fusion_transport_ledger_v1), target :: transport, cumulative_transport
    real(c_double), target :: reference_volume, accepted_volume, accepted_time
    integer(c_int64_t), target :: ticket, epoch
    integer(c_int) :: status
    type(c_ptr) :: state

    call make_initial(edges, initial_s, initial_t)
    candidate_s = initial_s
    candidate_s(1,1) = 8.0_c_double
    candidate_t = 0.0_c_double
    inert = 0.0_c_double
    call zero_source(source)
    call zero_transport(transport)
    state = c_null_ptr
    call fusion_source_state_create_volume(test_cells, edges, initial_s, &
         initial_t, 2.0_c_double, 0.0_c_double, test_tag + 3_c_int64_t, &
         state, status)
    call check_status(status, 'create snapshot-guard context')
    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin snapshot-guard context')
    call fusion_source_state_stage_volume(state, ticket, candidate_s, &
         candidate_t, 1.0_c_double, 1.0_c_double, source, inert, transport, &
         status)
    call check_status(status, 'stage snapshot-guard context')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit snapshot-guard context')

    bad_accepted_s = -1.0_c_double
    accepted_t = -1.0_c_double
    call zero_source(cumulative)
    cumulative%events_m3 = 1.0_c_double
    inert = 1.0_c_double
    call zero_transport(cumulative_transport)
    cumulative_transport%work_J = 1.0_c_double
    reference_volume = -1.0_c_double
    accepted_volume = -1.0_c_double
    accepted_time = -1.0_c_double
    epoch = -1_c_int64_t
    call fusion_source_state_snapshot_volume(state, bad_accepted_s, accepted_t, &
         cumulative, inert, cumulative_transport, reference_volume, &
         accepted_volume, accepted_time, epoch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'snapshot rejects wrong cell extent')
    call check(all(bad_accepted_s == 0.0_c_double) .and. &
         all(accepted_t == 0.0_c_double) .and. source_is_zero(cumulative) .and. &
         all(inert == 0.0_c_double) .and. &
         transport_is_zero(cumulative_transport) .and. &
         reference_volume == 0.0_c_double .and. accepted_volume == 0.0_c_double .and. &
         accepted_time == 0.0_c_double .and. epoch == 0_c_int64_t, &
         'invalid snapshot clears every output')
    call fusion_source_state_destroy(state)
  end subroutine test_bad_snapshot_clears_outputs

  subroutine test_volume_v3_restart()
    real(c_double), target :: edges(3), initial_s(2,6), initial_t(2,6)
    real(c_double), target :: candidate_s(2,6), candidate_t(2,6)
    real(c_double), target :: accepted_s(2,6), accepted_t(2,6)
    real(c_double), target :: inert(6)
    type(fusion_source_ledger_v1), target :: source, cumulative
    type(fusion_transport_ledger_v1), target :: transport, cumulative_transport
    real(c_double), target :: reference_volume, accepted_volume, accepted_time
    integer(c_int64_t), target :: ticket, epoch
    integer(c_int) :: status
    integer(c_size_t), target :: required, written
    integer(c_int8_t), allocatable, target :: buffer(:)
    type(c_ptr) :: state, restored

    call make_initial(edges, initial_s, initial_t)
    candidate_s = initial_s
    candidate_s(1,1) = 8.0_c_double
    candidate_t = 0.0_c_double
    inert = 0.0_c_double
    call zero_source(source)
    call zero_transport(transport)
    state = c_null_ptr
    call fusion_source_state_create_volume(test_cells, edges, initial_s, &
         initial_t, 2.0_c_double, 0.0_c_double, test_tag + 4_c_int64_t, &
         state, status)
    call check_status(status, 'create restart context')
    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin restart context')
    call fusion_source_state_stage_volume(state, ticket, candidate_s, &
         candidate_t, 1.0_c_double, 1.0_c_double, source, inert, transport, &
         status)
    call check_status(status, 'stage restart context')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit restart context')

    call fusion_source_state_pack_size(state, required, status)
    call check_status(status, 'query v3 packet size')
    allocate(buffer(required))
    call fusion_source_state_pack(state, buffer, required, written, status)
    call check_status(status, 'pack v3 volume context')
    call check(written == required .and. int(buffer(9), c_int) == 3_c_int, &
         'volume context uses restart version 3')

    restored = c_null_ptr
    call fusion_source_state_unpack(buffer, required, test_tag + 4_c_int64_t, &
         restored, status)
    call check_status(status, 'unpack v3 volume context')
    call check(c_associated(restored), 'v3 unpack returns a context')
    call fusion_source_state_snapshot_volume(restored, accepted_s, accepted_t, &
         cumulative, inert, cumulative_transport, reference_volume, &
         accepted_volume, accepted_time, epoch, status)
    call check_status(status, 'snapshot unpacked v3 context')
    call check(all(accepted_s == candidate_s) .and. all(accepted_t == candidate_t) &
         .and. reference_volume == 2.0_c_double .and. &
         accepted_volume == 1.0_c_double .and. epoch == 1_c_int64_t, &
         'unpacked v3 context preserves populations and volume state')
    call fusion_source_state_destroy(restored)
    call fusion_source_state_destroy(state)
    deallocate(buffer)
  end subroutine test_volume_v3_restart

end program test_fusion_source_volume_fortran
