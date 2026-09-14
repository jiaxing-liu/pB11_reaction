program test_fusion_source_state_fortran
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_int8_t, c_int64_t, c_null_ptr, c_ptr, c_size_t, c_sizeof
  use fusion_source_state_fortran
  implicit none

  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  integer(c_int), parameter :: test_cells = 2_c_int
  ! A negative c_int64_t carries a high-bit uint64_t model tag unchanged.
  integer(c_int64_t), parameter :: test_tag = -1122334455667788_c_int64_t
  integer :: failures

  failures = 0
  call test_initial_snapshot()
  call test_repeated_stage_replaces()
  call test_bad_extent_invalidates_stage()
  call test_discard()
  call test_pack_roundtrip_and_tag()
  call test_inert_stage_and_v2_restart()
  call test_invalid_extent_guards()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, &
          ' Fortran source-state binding test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran source-state binding tests passed'

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

  subroutine zero_ledger(ledger)
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
  end subroutine zero_ledger

  logical function ledger_is_zero(ledger)
    type(fusion_source_ledger_v1), intent(in) :: ledger

    ledger_is_zero = all(ledger%events_m3 == 0.0_c_double) .and. &
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
  end function ledger_is_zero

  subroutine make_initial(edges_J, initial_s_m3, initial_t_m3)
    real(c_double), intent(out), target :: edges_J(:)
    real(c_double), intent(out), target :: initial_s_m3(:)
    real(c_double), intent(out), target :: initial_t_m3(:)

    edges_J = [0.0_c_double, 2.0_c_double * mev_j, &
         4.0_c_double * mev_j]
    initial_s_m3 = 0.0_c_double
    initial_t_m3 = 0.0_c_double
    ! Species-major [species][cell], with He4 species ID 4 and cell 1.
    initial_s_m3(9) = 1.0e12_c_double
  end subroutine make_initial

  subroutine test_initial_snapshot()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: accepted_s(12), accepted_t(12)
    real(c_double), target :: accepted_inert(6)
    real(c_double), target :: accepted_time
    integer(c_int64_t), target :: epoch
    integer(c_int) :: cells, status
    type(fusion_source_ledger_v1), target :: cumulative
    type(c_ptr) :: state

    call make_initial(edges, initial_s, initial_t)
    state = c_null_ptr
    call zero_ledger(cumulative)
    call check(c_sizeof(cumulative) == &
         FUSION_SOURCE_STATE_LEDGER_DOUBLES * 8_c_size_t, &
         'Fortran ledger layout is exactly 121 C doubles')
    call fusion_source_state_create(test_cells, edges, initial_s, initial_t, &
         3.25_c_double, test_tag, state, status)
    call check_status(status, 'create source state')
    call check(c_associated(state), 'create returns an opaque context')

    cells = 0_c_int
    call fusion_source_state_cells(state, cells, status)
    call check_status(status, 'query immutable source-state cells')
    call check(cells == test_cells, 'cell query returns the configured count')

    call fusion_source_state_snapshot(state, accepted_s, accepted_t, &
         cumulative, accepted_time, epoch, status)
    call check_status(status, 'snapshot initial accepted state')
    call check(all(accepted_s == initial_s), &
         'initial snapshot returns S populations')
    call check(all(accepted_t == initial_t), &
         'initial snapshot returns T populations')
    call check(ledger_is_zero(cumulative), 'initial cumulative ledger is zero')
    call check(accepted_time == 3.25_c_double, &
         'initial snapshot returns time')
    call check(epoch == 0_c_int64_t, 'initial accepted epoch is zero')

    call fusion_source_state_snapshot_inert(state, accepted_s, accepted_t, &
         cumulative, accepted_inert, accepted_time, epoch, status)
    call check_status(status, 'extended snapshot of legacy context')
    call check(all(accepted_inert == 0.0_c_double), &
         'extended snapshot reports zero inert heat for legacy state')

    call fusion_source_state_destroy(state)
    call check(.not. c_associated(state), 'destroy clears the opaque handle')
  end subroutine test_initial_snapshot

  subroutine test_repeated_stage_replaces()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: first_s(12), second_s(12), trial_t(12)
    real(c_double), target :: accepted_s(12), accepted_t(12)
    real(c_double), target :: accepted_time
    integer(c_int64_t), target :: ticket, epoch
    integer(c_int) :: status
    type(fusion_source_ledger_v1), target :: first_step, second_step
    type(fusion_source_ledger_v1), target :: cumulative
    type(c_ptr) :: state

    call make_initial(edges, initial_s, initial_t)
    state = c_null_ptr
    call fusion_source_state_create(test_cells, edges, initial_s, initial_t, &
         3.25_c_double, test_tag + 1_c_int64_t, state, status)
    call check_status(status, 'create replacement-stage state')

    first_s = initial_s
    first_s(9) = initial_s(9) - 1.0e11_c_double
    trial_t = initial_t
    call zero_ledger(first_step)
    first_step%escaped_number_m3(5) = 1.0e11_c_double
    first_step%escaped_energy_J_m3(5) = 1.0e11_c_double * mev_j

    second_s = initial_s
    second_s(9) = initial_s(9) - 2.0e11_c_double
    call zero_ledger(second_step)
    second_step%escaped_number_m3(5) = 2.0e11_c_double
    second_step%escaped_energy_J_m3(5) = 2.0e11_c_double * mev_j

    call fusion_source_state_begin(state, 0.1_c_double, ticket, status)
    call check_status(status, 'begin replacement-stage trial')
    call fusion_source_state_stage(state, ticket, first_s, trial_t, &
         first_step, status)
    call check_status(status, 'stage first nonlinear candidate')
    call fusion_source_state_stage(state, ticket, second_s, trial_t, &
         second_step, status)
    call check_status(status, 'stage replacement nonlinear candidate')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit replacement candidate')

    call fusion_source_state_snapshot(state, accepted_s, accepted_t, &
         cumulative, accepted_time, epoch, status)
    call check_status(status, 'snapshot replacement result')
    call check(all(accepted_s == second_s), &
         'repeated stage commits the replacement S state')
    call check(all(accepted_t == trial_t), &
         'repeated stage commits the replacement T state')
    call check(cumulative%escaped_number_m3(5) == 2.0e11_c_double, &
         'repeated stage does not accumulate the first candidate')
    call check(cumulative%escaped_energy_J_m3(5) == &
         2.0e11_c_double * mev_j, &
         'replacement ledger keeps only the second energy amount')
    call check(accepted_time == 3.25_c_double + 0.1_c_double, &
         'replacement commit advances time once')
    call check(epoch == 1_c_int64_t, &
         'replacement commit increments epoch once')

    call fusion_source_state_destroy(state)
  end subroutine test_repeated_stage_replaces

  subroutine test_bad_extent_invalidates_stage()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: candidate_s(12), candidate_t(12)
    real(c_double), target :: bad_s(11)
    integer(c_int64_t), target :: ticket
    integer(c_int) :: status
    type(fusion_source_ledger_v1), target :: step
    type(c_ptr) :: state

    call make_initial(edges, initial_s, initial_t)
    state = c_null_ptr
    call fusion_source_state_create(test_cells, edges, initial_s, initial_t, &
         0.0_c_double, test_tag + 2_c_int64_t, state, status)
    call check_status(status, 'create wrong-extent regression state')
    candidate_s = initial_s
    candidate_s(9) = initial_s(9) - 1.0e11_c_double
    candidate_t = initial_t
    call zero_ledger(step)
    step%escaped_number_m3(5) = 1.0e11_c_double
    step%escaped_energy_J_m3(5) = 1.0e11_c_double * mev_j

    call fusion_source_state_begin(state, 0.2_c_double, ticket, status)
    call check_status(status, 'begin wrong-extent regression trial')
    call fusion_source_state_stage(state, ticket, candidate_s, candidate_t, &
         step, status)
    call check_status(status, 'stage valid candidate before bad extent')

    bad_s = 0.0_c_double
    call fusion_source_state_stage(state, ticket, bad_s, candidate_t, step, &
         status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong stage extent is rejected before C_LOC')
    call fusion_source_state_commit(state, ticket, status)
    call check(status /= PB11_STATUS_OK, &
         'wrong stage extent invalidates the prior staged candidate')
    call fusion_source_state_discard(state, ticket, status)
    call check_status(status, 'discard ticket after bad stage')

    call fusion_source_state_destroy(state)
  end subroutine test_bad_extent_invalidates_stage

  subroutine test_discard()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: candidate_s(12), accepted_s(12), accepted_t(12)
    real(c_double), target :: accepted_time
    integer(c_int64_t), target :: ticket, epoch
    integer(c_int) :: status
    type(fusion_source_ledger_v1), target :: step, cumulative
    type(c_ptr) :: state

    call make_initial(edges, initial_s, initial_t)
    state = c_null_ptr
    call fusion_source_state_create(test_cells, edges, initial_s, initial_t, &
         2.0_c_double, test_tag + 3_c_int64_t, state, status)
    call check_status(status, 'create discard state')
    candidate_s = initial_s
    candidate_s(9) = initial_s(9) - 1.0e11_c_double
    call zero_ledger(step)
    step%escaped_number_m3(5) = 1.0e11_c_double
    step%escaped_energy_J_m3(5) = 1.0e11_c_double * mev_j

    call fusion_source_state_begin(state, 0.5_c_double, ticket, status)
    call check_status(status, 'begin discard trial')
    call fusion_source_state_stage(state, ticket, candidate_s, initial_t, &
         step, status)
    call check_status(status, 'stage candidate to discard')
    call fusion_source_state_discard(state, ticket, status)
    call check_status(status, 'discard pending candidate')
    call fusion_source_state_commit(state, ticket, status)
    call check(status /= PB11_STATUS_OK, 'discarded ticket cannot commit')

    call fusion_source_state_snapshot(state, accepted_s, accepted_t, &
         cumulative, accepted_time, epoch, status)
    call check_status(status, 'snapshot after discard')
    call check(all(accepted_s == initial_s) .and. all(accepted_t == initial_t), &
         'discard leaves accepted populations unchanged')
    call check(ledger_is_zero(cumulative), &
         'discard leaves cumulative ledger unchanged')
    call check(accepted_time == 2.0_c_double .and. epoch == 0_c_int64_t, &
         'discard leaves accepted time and epoch unchanged')

    call fusion_source_state_destroy(state)
  end subroutine test_discard

  subroutine test_pack_roundtrip_and_tag()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: candidate_s(12), original_s(12), original_t(12)
    real(c_double), target :: restored_s(12), restored_t(12)
    real(c_double), target :: original_time, restored_time
    integer(c_int64_t), target :: ticket, original_epoch, restored_epoch
    integer(c_int) :: status
    integer(c_size_t), target :: required, written
    integer(c_int8_t), allocatable, target :: buffer(:)
    type(fusion_source_ledger_v1), target :: step
    type(fusion_source_ledger_v1), target :: original_ledger, restored_ledger
    type(c_ptr) :: state, restored, mismatch

    call make_initial(edges, initial_s, initial_t)
    state = c_null_ptr
    call fusion_source_state_create(test_cells, edges, initial_s, initial_t, &
         4.0_c_double, test_tag + 4_c_int64_t, state, status)
    call check_status(status, 'create pack state')
    candidate_s = initial_s
    candidate_s(9) = initial_s(9) - 2.0e11_c_double
    call zero_ledger(step)
    step%escaped_number_m3(5) = 2.0e11_c_double
    step%escaped_energy_J_m3(5) = 2.0e11_c_double * mev_j

    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin state before packing')
    call fusion_source_state_stage(state, ticket, candidate_s, initial_t, &
         step, status)
    call check_status(status, 'stage state before packing')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit state before packing')

    call fusion_source_state_snapshot(state, original_s, original_t, &
         original_ledger, original_time, original_epoch, status)
    call check_status(status, 'snapshot state before packing')
    call fusion_source_state_pack_size(state, required, status)
    call check_status(status, 'query accepted restart packet size')
    call check(required > 0_c_size_t, 'restart packet size is positive')
    allocate(buffer(required))
    call fusion_source_state_pack(state, buffer, required, written, status)
    call check_status(status, 'pack accepted restart state')
    call check(written == required, 'pack writes the required byte count')

    restored = c_null_ptr
    call fusion_source_state_unpack(buffer, required, test_tag + 4_c_int64_t, &
         restored, status)
    call check_status(status, 'unpack accepted restart state')
    call check(c_associated(restored), 'unpack returns a new context')
    call fusion_source_state_snapshot(restored, restored_s, restored_t, &
         restored_ledger, restored_time, restored_epoch, status)
    call check_status(status, 'snapshot unpacked restart state')
    call check(all(restored_s == original_s) .and. all(restored_t == original_t), &
         'unpacked populations equal the accepted original')
    call check(ledger_equal(restored_ledger, original_ledger) .and. &
         restored_time == original_time .and. restored_epoch == original_epoch, &
         'unpacked ledger, time, and epoch equal the original')

    ! The unpacked context remains usable for a fresh accepted step.
    call zero_ledger(step)
    call fusion_source_state_begin(restored, 0.125_c_double, ticket, status)
    call check_status(status, 'begin step after unpack')
    call fusion_source_state_stage(restored, ticket, restored_s, restored_t, &
         step, status)
    call check_status(status, 'stage step after unpack')
    call fusion_source_state_commit(restored, ticket, status)
    call check_status(status, 'commit step after unpack')

    mismatch = c_null_ptr
    call fusion_source_state_unpack(buffer, required, test_tag + 5_c_int64_t, &
         mismatch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'unpack rejects a model-tag mismatch')
    call check(.not. c_associated(mismatch), &
         'tag mismatch does not publish a context')

    call fusion_source_state_destroy(restored)
    call fusion_source_state_destroy(state)
    deallocate(buffer)
  end subroutine test_pack_roundtrip_and_tag

  subroutine test_inert_stage_and_v2_restart()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: candidate_s(12), candidate_t(12)
    real(c_double), target :: accepted_s(12), accepted_t(12)
    real(c_double), target :: restored_s(12), restored_t(12)
    real(c_double), target :: bad_s(11), bad_inert(5)
    real(c_double), target :: inert_step(6), zero_inert(6)
    real(c_double), target :: accepted_inert(6), restored_inert(6)
    real(c_double), target :: accepted_time, restored_time
    integer(c_int64_t), target :: ticket, accepted_epoch, restored_epoch
    integer(c_int) :: status
    integer(c_size_t), target :: required, written
    integer(c_int8_t), allocatable, target :: buffer(:)
    type(fusion_source_ledger_v1), target :: step, cumulative
    type(fusion_source_ledger_v1), target :: restored_cumulative
    type(c_ptr) :: state, restored, mismatch

    call make_initial(edges, initial_s, initial_t)
    ! Two alpha energy cells lose 1e12 MeV to the inert bath.  A proton is
    ! shifted from the first to the second cell, requiring signed negative
    ! inert heat for that species.  Particle numbers stay unchanged.
    initial_s = 0.0_c_double
    initial_t = 0.0_c_double
    initial_s(1) = 1.0e11_c_double
    initial_s(9) = 1.0e12_c_double
    initial_s(10) = 1.0e12_c_double
    candidate_s = initial_s
    candidate_s(1) = 0.0_c_double
    candidate_s(2) = 1.0e11_c_double
    candidate_s(9) = 1.5e12_c_double
    candidate_s(10) = 0.5e12_c_double
    candidate_t = initial_t
    call zero_ledger(step)
    zero_inert = 0.0_c_double
    inert_step = 0.0_c_double
    inert_step(1) = -2.0e11_c_double * mev_j
    inert_step(5) = 1.0e12_c_double * mev_j

    state = c_null_ptr
    call fusion_source_state_create(test_cells, edges, initial_s, initial_t, &
         5.0_c_double, test_tag + 8_c_int64_t, state, status)
    call check_status(status, 'create inert-heat state')
    call fusion_source_state_begin(state, 0.5_c_double, ticket, status)
    call check_status(status, 'begin inert-heat trial')
    call fusion_source_state_stage_inert(state, ticket, candidate_s, &
         candidate_t, step, inert_step, status)
    call check_status(status, 'stage signed inert-heat candidate')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit inert-heat candidate')

    call fusion_source_state_snapshot_inert(state, accepted_s, accepted_t, &
         cumulative, accepted_inert, accepted_time, accepted_epoch, status)
    call check_status(status, 'snapshot extended accepted state')
    call check(all(accepted_s == candidate_s) .and. &
         all(accepted_t == candidate_t), &
         'extended snapshot returns the two-cell accepted populations')
    call check(accepted_inert(1) == inert_step(1) .and. &
         accepted_inert(5) == inert_step(5), &
         'extended snapshot returns signed inert-bath heat')
    call check(accepted_time == 5.0_c_double + 0.5_c_double .and. &
         accepted_epoch == 1_c_int64_t, &
         'extended commit advances time and epoch once')

    bad_inert = -1.0_c_double
    accepted_s = -1.0_c_double
    accepted_t = -1.0_c_double
    cumulative%events_m3 = 1.0_c_double
    accepted_time = -1.0_c_double
    accepted_epoch = -1_c_int64_t
    call fusion_source_state_snapshot_inert(state, accepted_s, accepted_t, &
         cumulative, bad_inert, accepted_time, accepted_epoch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'extended snapshot rejects a wrong inert-heat extent')
    call check(all(accepted_s == 0.0_c_double) .and. &
         all(accepted_t == 0.0_c_double) .and. ledger_is_zero(cumulative) .and. &
         all(bad_inert == 0.0_c_double) .and. &
         accepted_time == 0.0_c_double .and. accepted_epoch == 0_c_int64_t, &
         'extended snapshot extent failure clears every Fortran output')

    ! The v1 snapshot cannot silently omit inert-bath accounting after an
    ! extended commit.  The Fortran wrapper clears every output on failure.
    accepted_s = -1.0_c_double
    accepted_t = -1.0_c_double
    cumulative%events_m3 = 1.0_c_double
    accepted_time = -1.0_c_double
    accepted_epoch = -1_c_int64_t
    call fusion_source_state_snapshot(state, accepted_s, accepted_t, &
         cumulative, accepted_time, accepted_epoch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'legacy snapshot rejects an extended context')
    call check(all(accepted_s == 0.0_c_double) .and. &
         all(accepted_t == 0.0_c_double) .and. ledger_is_zero(cumulative) .and. &
         accepted_time == 0.0_c_double .and. accepted_epoch == 0_c_int64_t, &
         'legacy snapshot failure clears all Fortran outputs')

    ! A malformed extended stage must cancel a previously valid candidate even
    ! though its bad extent cannot be passed to C_LOC.
    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin inert wrong-extent trial')
    call fusion_source_state_stage_inert(state, ticket, candidate_s, candidate_t, &
         step, zero_inert, status)
    call check_status(status, 'stage candidate before inert extent failure')
    bad_inert = -1.0_c_double
    call fusion_source_state_stage_inert(state, ticket, candidate_s, candidate_t, &
         step, bad_inert, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong inert-heat extent is rejected before C_LOC')
    call fusion_source_state_commit(state, ticket, status)
    call check(status /= PB11_STATUS_OK, &
         'wrong inert-heat extent invalidates the prior candidate')
    call fusion_source_state_discard(state, ticket, status)
    call check_status(status, 'discard ticket after inert-heat extent failure')

    call fusion_source_state_begin(state, 0.25_c_double, ticket, status)
    call check_status(status, 'begin inert state-array extent trial')
    call fusion_source_state_stage_inert(state, ticket, candidate_s, candidate_t, &
         step, zero_inert, status)
    call check_status(status, 'stage before inert state-array failure')
    bad_s = -1.0_c_double
    call fusion_source_state_stage_inert(state, ticket, bad_s, candidate_t, &
         step, zero_inert, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong inert-stage state extent is rejected before C_LOC')
    call fusion_source_state_commit(state, ticket, status)
    call check(status /= PB11_STATUS_OK, &
         'wrong inert-stage extent invalidates the prior candidate')
    call fusion_source_state_discard(state, ticket, status)
    call check_status(status, 'discard inert ticket after extent failure')

    ! Packing an extended context produces v2 and unpacking preserves the
    ! separate inert account.  The restored context must remain stageable.
    call fusion_source_state_pack_size(state, required, status)
    call check_status(status, 'query extended restart packet size')
    call check(required > 0_c_size_t, 'extended restart packet size is positive')
    allocate(buffer(required))
    call fusion_source_state_pack(state, buffer, required, written, status)
    call check_status(status, 'pack extended restart state')
    call check(written == required, 'extended pack writes the required bytes')

    restored = c_null_ptr
    call fusion_source_state_unpack(buffer, required, test_tag + 8_c_int64_t, &
         restored, status)
    call check_status(status, 'unpack extended v2 restart state')
    call check(c_associated(restored), 'v2 unpack returns a context')
    call fusion_source_state_snapshot_inert(restored, restored_s, restored_t, &
         restored_cumulative, restored_inert, restored_time, restored_epoch, &
         status)
    call check_status(status, 'snapshot unpacked extended state')
    call check(all(restored_s == candidate_s) .and. &
         all(restored_t == candidate_t) .and. &
         ledger_equal(restored_cumulative, cumulative) .and. &
         all(restored_inert == inert_step) .and. &
         restored_time == 5.0_c_double + 0.5_c_double .and. &
         restored_epoch == 1_c_int64_t, &
         'v2 round trip preserves populations, ledger, inert heat, and epoch')

    zero_inert = 0.0_c_double
    call zero_ledger(step)
    call fusion_source_state_begin(restored, 0.125_c_double, ticket, status)
    call check_status(status, 'begin continuation after v2 unpack')
    call fusion_source_state_stage_inert(restored, ticket, restored_s, &
         restored_t, step, zero_inert, status)
    call check_status(status, 'stage continuation after v2 unpack')
    call fusion_source_state_commit(restored, ticket, status)
    call check_status(status, 'commit continuation after v2 unpack')
    call fusion_source_state_snapshot_inert(restored, restored_s, restored_t, &
         restored_cumulative, restored_inert, restored_time, restored_epoch, &
         status)
    call check_status(status, 'snapshot continuation after v2 unpack')
    call check(restored_time == 5.0_c_double + 0.5_c_double + 0.125_c_double .and. &
         restored_epoch == 2_c_int64_t .and. all(restored_inert == inert_step), &
         'v2 context continues with cumulative inert heat')

    mismatch = c_null_ptr
    call fusion_source_state_unpack(buffer, required, test_tag + 9_c_int64_t, &
         mismatch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         .not. c_associated(mismatch), &
         'v2 unpack rejects a model-tag mismatch')

    call fusion_source_state_destroy(restored)
    call fusion_source_state_destroy(state)
    deallocate(buffer)
  end subroutine test_inert_stage_and_v2_restart

  subroutine test_invalid_extent_guards()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: wrong_s(11), accepted_t(12)
    real(c_double), target :: accepted_time
    integer(c_int64_t), target :: epoch
    integer(c_int) :: status
    integer(c_int8_t), allocatable, target :: empty_buffer(:)
    integer(c_size_t), target :: written
    type(fusion_source_ledger_v1), target :: cumulative
    type(c_ptr) :: state, rejected

    call make_initial(edges, initial_s, initial_t)
    state = c_null_ptr
    call fusion_source_state_create(test_cells, edges, wrong_s, initial_t, &
         0.0_c_double, test_tag + 6_c_int64_t, state, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'create rejects a wrong state-array extent before C_LOC')
    call check(.not. c_associated(state), &
         'failed create leaves its context output null')

    state = c_null_ptr
    call fusion_source_state_create(test_cells, edges, initial_s, initial_t, &
         0.0_c_double, test_tag + 7_c_int64_t, state, status)
    call check_status(status, 'create state for snapshot extent guard')
    wrong_s = -1.0_c_double
    accepted_t = -1.0_c_double
    cumulative%events_m3 = 1.0_c_double
    accepted_time = -1.0_c_double
    epoch = -1_c_int64_t
    call fusion_source_state_snapshot(state, wrong_s, accepted_t, cumulative, &
         accepted_time, epoch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'snapshot rejects wrong extents before C_LOC')
    call check(all(wrong_s == 0.0_c_double) .and. &
         all(accepted_t == 0.0_c_double) .and. ledger_is_zero(cumulative) .and. &
         accepted_time == 0.0_c_double .and. epoch == 0_c_int64_t, &
         'invalid snapshot clears all outputs')

    allocate(empty_buffer(0))
    written = -1_c_size_t
    call fusion_source_state_pack(state, empty_buffer, 0_c_size_t, written, &
         status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         written == 0_c_size_t, &
         'pack rejects zero capacity without forming C_LOC')
    rejected = c_null_ptr
    call fusion_source_state_unpack(empty_buffer, 0_c_size_t, &
         test_tag + 7_c_int64_t, rejected, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         .not. c_associated(rejected), &
         'unpack rejects a zero-length buffer without forming C_LOC')

    call fusion_source_state_destroy(state)
    deallocate(empty_buffer)
  end subroutine test_invalid_extent_guards

  logical function ledger_equal(left, right)
    type(fusion_source_ledger_v1), intent(in) :: left, right

    ledger_equal = all(left%events_m3 == right%events_m3) .and. &
         all(left%nuclear_born_number_m3 == right%nuclear_born_number_m3) .and. &
         all(left%nuclear_born_energy_J_m3 == right%nuclear_born_energy_J_m3) .and. &
         all(left%external_born_number_m3 == right%external_born_number_m3) .and. &
         all(left%external_born_energy_J_m3 == right%external_born_energy_J_m3) .and. &
         all(left%thermal_consumed_number_m3 == right%thermal_consumed_number_m3) .and. &
         all(left%thermal_consumed_energy_J_m3 == right%thermal_consumed_energy_J_m3) .and. &
         all(left%fast_consumed_number_m3 == right%fast_consumed_number_m3) .and. &
         all(left%fast_consumed_energy_J_m3 == right%fast_consumed_energy_J_m3) .and. &
         all(left%escaped_number_m3 == right%escaped_number_m3) .and. &
         all(left%escaped_energy_J_m3 == right%escaped_energy_J_m3) .and. &
         all(left%handed_off_number_m3 == right%handed_off_number_m3) .and. &
         all(left%handed_off_energy_J_m3 == right%handed_off_energy_J_m3) .and. &
         all(left%heat_to_bath_J_m3 == right%heat_to_bath_J_m3) .and. &
         left%neutron_number_m3 == right%neutron_number_m3 .and. &
         left%neutron_energy_J_m3 == right%neutron_energy_J_m3
  end function ledger_equal

end program test_fusion_source_state_fortran
