program test_fusion_source_state_numerical_fortran
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_int8_t, c_int64_t, c_loc, c_null_ptr, c_ptr, c_size_t
  use fusion_source_state_fortran, only : &
       PB11_STATUS_OK, PB11_STATUS_INVALID_ARGUMENT, &
       fusion_source_ledger_v1, fusion_source_state_begin, &
       fusion_source_state_commit, fusion_source_state_create, &
       fusion_source_state_destroy, fusion_source_state_discard, &
       fusion_source_state_pack, fusion_source_state_pack_size, &
       fusion_source_state_snapshot_inert, fusion_source_state_stage_inert, &
       fusion_source_state_unpack
  use fusion_source_volume_fortran, only : &
       fusion_transport_ledger_v1, fusion_source_state_create_volume, &
       fusion_source_state_snapshot_volume
  use fusion_source_state_numerical_fortran, only : &
       c_fusion_source_state_stage_numerical, &
       c_fusion_source_state_stage_volume_numerical, &
       c_fusion_source_state_snapshot_numerical, &
       c_fusion_source_state_snapshot_volume_numerical
  implicit none

  integer(c_int), parameter :: cells = 2_c_int
  integer(c_int64_t), parameter :: tag = 107_c_int64_t
  integer :: failures

  failures = 0
  call test_stationary()
  call test_moving()

  if (failures /= 0) then
     write (*, '(I0,A)') failures, &
          ' Fortran numerical binding test(s) failed'
     error stop 1
  end if
  write (*, '(A)') &
       'PASS numerical Fortran binding: stationary/moving remap, old API rejection, restart'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write (*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  subroutine check_status(status, message)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: message

    call check(status == PB11_STATUS_OK, message)
  end subroutine check_status

  subroutine zero_source(source)
    type(fusion_source_ledger_v1), intent(out) :: source

    source%events_m3 = 0.0_c_double
    source%nuclear_born_number_m3 = 0.0_c_double
    source%nuclear_born_energy_J_m3 = 0.0_c_double
    source%external_born_number_m3 = 0.0_c_double
    source%external_born_energy_J_m3 = 0.0_c_double
    source%thermal_consumed_number_m3 = 0.0_c_double
    source%thermal_consumed_energy_J_m3 = 0.0_c_double
    source%fast_consumed_number_m3 = 0.0_c_double
    source%fast_consumed_energy_J_m3 = 0.0_c_double
    source%escaped_number_m3 = 0.0_c_double
    source%escaped_energy_J_m3 = 0.0_c_double
    source%handed_off_number_m3 = 0.0_c_double
    source%handed_off_energy_J_m3 = 0.0_c_double
    source%heat_to_bath_J_m3 = 0.0_c_double
    source%neutron_number_m3 = 0.0_c_double
    source%neutron_energy_J_m3 = 0.0_c_double
  end subroutine zero_source

  subroutine zero_transport(transport)
    type(fusion_transport_ledger_v1), intent(out) :: transport

    transport%spatial_number = 0.0_c_double
    transport%spatial_energy_J = 0.0_c_double
    transport%work_J = 0.0_c_double
    transport%lower_number = 0.0_c_double
    transport%lower_energy_J = 0.0_c_double
    transport%upper_number = 0.0_c_double
    transport%upper_energy_J = 0.0_c_double
  end subroutine zero_transport

  subroutine test_stationary()
    real(c_double), target :: edges(3), initial_s(12), initial_t(12)
    real(c_double), target :: trial_s(12), trial_t(12)
    real(c_double), target :: accepted_s(12), accepted_t(12)
    real(c_double), target :: inert(6), numerical(6), out_inert(6)
    real(c_double), target :: accepted_time
    integer(c_int64_t), target :: ticket, epoch
    integer(c_int) :: status
    type(fusion_source_ledger_v1), target :: source, cumulative
    type(c_ptr) :: state, restored
    integer(c_size_t), target :: required, written
    integer(c_int8_t), allocatable, target :: packet(:)

    edges = [0.0_c_double, 2.0_c_double, 4.0_c_double]
    initial_s = 0.0_c_double
    initial_t = 0.0_c_double
    ! Arithmetic cell centres are 1 and 3 J.  Splitting one particle 50/50
    ! across those cells remaps N=1 at E=1 to N=1 at E=2.
    initial_s(1) = 1.0_c_double
    trial_s = 0.0_c_double
    trial_s(1) = 0.5_c_double
    trial_s(2) = 0.5_c_double
    trial_t = 0.0_c_double
    inert = 0.0_c_double
    numerical = 0.0_c_double
    numerical(1) = -1.0_c_double
    call zero_source(source)

    state = c_null_ptr
    call fusion_source_state_create(cells, edges, initial_s, initial_t, &
         0.0_c_double, tag, state, status)
    call check_status(status, 'create stationary numerical owner')
    call check(c_associated(state), 'stationary owner is associated')

    call fusion_source_state_begin(state, 1.0_c_double, ticket, status)
    call check_status(status, 'begin stationary remap')
    status = c_fusion_source_state_stage_numerical(state, ticket, &
         c_loc(trial_s(1)), c_loc(trial_t(1)), c_loc(source), &
         c_loc(inert(1)), c_loc(numerical(1)))
    call check_status(status, 'stage stationary signed numerical remap')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit stationary signed numerical remap')

    accepted_s = -7.0_c_double
    accepted_t = -7.0_c_double
    call zero_source(cumulative)
    out_inert = -7.0_c_double
    accepted_time = -7.0_c_double
    epoch = -7_c_int64_t
    status = c_fusion_source_state_snapshot_numerical(state, &
         c_loc(accepted_s(1)), c_loc(accepted_t(1)), c_loc(cumulative), &
         c_loc(out_inert(1)), c_loc(numerical(1)), c_loc(accepted_time), &
         c_loc(epoch))
    call check_status(status, 'snapshot stationary numerical account')
    call check(all(accepted_s == trial_s) .and. all(accepted_t == trial_t), &
         'stationary snapshot returns the remapped populations')
    call check(numerical(1) == -1.0_c_double .and. &
         all(numerical(2:6) == 0.0_c_double), &
         'stationary snapshot returns delta E=-1 numerical account')
    call check(accepted_time == 1.0_c_double .and. epoch == 1_c_int64_t, &
         'stationary snapshot returns time and uint64 epoch')

    ! A numerical commit promotes the owner.  Existing inert stage/snapshot
    ! calls must reject it instead of silently dropping the new account.
    call fusion_source_state_snapshot_inert(state, accepted_s, accepted_t, &
         cumulative, out_inert, accepted_time, epoch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'old stationary snapshot rejects promoted numerical owner')
    call fusion_source_state_begin(state, 1.0_c_double, ticket, status)
    call check_status(status, 'begin old-stage rejection probe')
    call fusion_source_state_stage_inert(state, ticket, trial_s, trial_t, &
         source, inert, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'old stationary stage rejects promoted numerical owner')
    call fusion_source_state_discard(state, ticket, status)
    call check_status(status, 'discard rejected old-stage ticket')

    ! Pack/unpack through the existing lifecycle wrappers.  The restored
    ! numerical account is read with the new direct snapshot declaration.
    call fusion_source_state_pack_size(state, required, status)
    call check_status(status, 'query promoted stationary restart size')
    allocate(packet(required))
    call fusion_source_state_pack(state, packet, required, written, status)
    call check_status(status, 'pack promoted stationary owner')
    call check(written == required, 'stationary restart writes all bytes')
    restored = c_null_ptr
    call fusion_source_state_unpack(packet, required, tag, restored, status)
    call check_status(status, 'unpack promoted stationary owner')
    call check(c_associated(restored), 'unpack returns stationary owner')
    numerical = 0.0_c_double
    status = c_fusion_source_state_snapshot_numerical(restored, &
         c_loc(accepted_s(1)), c_loc(accepted_t(1)), c_loc(cumulative), &
         c_loc(out_inert(1)), c_loc(numerical(1)), c_loc(accepted_time), &
         c_loc(epoch))
    call check_status(status, 'snapshot restored stationary numerical account')
    call check(numerical(1) == -1.0_c_double .and. epoch == 1_c_int64_t, &
         'stationary numerical account survives restart')

    call fusion_source_state_destroy(restored)
    call fusion_source_state_destroy(state)
    deallocate(packet)
  end subroutine test_stationary

  subroutine test_moving()
    real(c_double), target :: edges(3)
    real(c_double), target :: initial_s(2,6), initial_t(2,6)
    real(c_double), target :: trial_s(2,6), trial_t(2,6)
    real(c_double), target :: accepted_s(2,6), accepted_t(2,6)
    real(c_double), target :: inert(6), numerical(6), out_inert(6)
    real(c_double), target :: reference_volume, accepted_volume, accepted_time
    integer(c_int64_t), target :: ticket, epoch
    integer(c_int) :: status
    type(fusion_source_ledger_v1), target :: source, cumulative
    type(fusion_transport_ledger_v1), target :: transport, cumulative_transport
    type(c_ptr) :: state

    edges = [0.0_c_double, 2.0_c_double, 4.0_c_double]
    initial_s = 0.0_c_double
    initial_t = 0.0_c_double
    initial_s(1,1) = 1.0_c_double
    trial_s = 0.0_c_double
    trial_s(1,1) = 0.5_c_double
    trial_s(2,1) = 0.5_c_double
    trial_t = 0.0_c_double
    inert = 0.0_c_double
    numerical = 0.0_c_double
    numerical(1) = -1.0_c_double
    call zero_source(source)
    call zero_transport(transport)

    state = c_null_ptr
    call fusion_source_state_create_volume(cells, edges, initial_s, initial_t, &
         2.0_c_double, 0.0_c_double, tag + 1_c_int64_t, state, status)
    call check_status(status, 'create moving numerical owner')
    call check(c_associated(state), 'moving owner is associated')

    call fusion_source_state_begin(state, 1.0_c_double, ticket, status)
    call check_status(status, 'begin moving remap')
    status = c_fusion_source_state_stage_volume_numerical(state, ticket, &
         c_loc(trial_s(1,1)), c_loc(trial_t(1,1)), 2.0_c_double, &
         2.0_c_double, c_loc(source), c_loc(inert(1)), c_loc(transport), &
         c_loc(numerical(1)))
    call check_status(status, 'stage moving signed numerical remap')
    call fusion_source_state_commit(state, ticket, status)
    call check_status(status, 'commit moving signed numerical remap')

    accepted_s = -7.0_c_double
    accepted_t = -7.0_c_double
    call zero_source(cumulative)
    call zero_transport(cumulative_transport)
    out_inert = -7.0_c_double
    reference_volume = -7.0_c_double
    accepted_volume = -7.0_c_double
    accepted_time = -7.0_c_double
    epoch = -7_c_int64_t
    status = c_fusion_source_state_snapshot_volume_numerical(state, &
         c_loc(accepted_s(1,1)), c_loc(accepted_t(1,1)), c_loc(cumulative), &
         c_loc(out_inert(1)), c_loc(cumulative_transport), &
         c_loc(reference_volume), c_loc(accepted_volume), c_loc(accepted_time), &
         c_loc(epoch), c_loc(numerical(1)))
    call check_status(status, 'snapshot moving numerical account')
    call check(all(accepted_s == trial_s) .and. all(accepted_t == trial_t), &
         'moving snapshot returns the remapped populations')
    call check(numerical(1) == -1.0_c_double .and. &
         reference_volume == 2.0_c_double .and. accepted_volume == 2.0_c_double, &
         'moving snapshot preserves scaled numerical account and volumes')
    call check(accepted_time == 1.0_c_double .and. epoch == 1_c_int64_t, &
         'moving snapshot returns time and uint64 epoch')

    call fusion_source_state_snapshot_volume(state, accepted_s, accepted_t, &
         cumulative, out_inert, cumulative_transport, reference_volume, &
         accepted_volume, accepted_time, epoch, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'old moving snapshot rejects promoted numerical owner')

    call fusion_source_state_destroy(state)
  end subroutine test_moving

end program test_fusion_source_state_numerical_fortran
