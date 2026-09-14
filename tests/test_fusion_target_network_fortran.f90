program test_fusion_target_network_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_sizeof
  use fusion_target_network_fortran
  implicit none

  integer :: failures

  failures = 0
  call test_layout()
  call test_golden_root()
  call test_zero_edges()
  call test_shape_rejection()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, &
          ' Fortran target-network test assertion(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran target-network tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  subroutine check_close(actual, expected, message)
    real(c_double), intent(in) :: actual, expected
    character(len=*), intent(in) :: message
    real(c_double) :: scale

    scale = max(abs(actual), abs(expected))
    call check(abs(actual - expected) <= 5.0e-12_c_double * scale + &
         5.0e-13_c_double, message)
  end subroutine check_close

  subroutine check_network_zero(out, message)
    type(fusion_target_network_v1), intent(in) :: out
    character(len=*), intent(in) :: message

    call check_close(out%reactions_m3, 0.0_c_double, message // &
         ': reactions cleared')
    call check_close(out%removed_fast_energy_J_m3, 0.0_c_double, message // &
         ': fast-energy ledger cleared')
    call check_close(out%removed_target_energy_J_m3, 0.0_c_double, message // &
         ': target-energy ledger cleared')
    call check_close(out%max_fast_number_relative_residual, 0.0_c_double, &
         message // ': fast residual cleared')
    call check_close(out%max_target_number_relative_residual, &
         0.0_c_double, message // ': target residual cleared')
    call check_close(out%max_target_energy_relative_residual, &
         0.0_c_double, message // ': target-energy residual cleared')
    call check(out%iterations == 0_c_int, message // &
         ': iteration count cleared')
  end subroutine check_network_zero

  subroutine test_layout()
    type(fusion_target_network_edge_v1) :: edge
    type(fusion_target_network_v1) :: out

    call check(c_sizeof(edge) == 24, 'edge bind(C) size is 24 bytes')
    call check(c_sizeof(out) == 56, 'network bind(C) size is 56 bytes')
  end subroutine test_layout

  subroutine test_golden_root()
    real(c_double), target :: fast_energy(1), old_fast(1)
    real(c_double), target :: target_number(1), target_energy(1)
    real(c_double), target :: trial_fast(1), trial_target(1)
    real(c_double), target :: trial_target_energy(1), edge_events(1)
    type(fusion_target_network_edge_v1), target :: edges(1)
    type(fusion_target_network_v1), target :: out
    integer(c_int) :: status
    real(c_double) :: root, reaction

    fast_energy = [5.0_c_double]
    old_fast = [1.0_c_double]
    target_number = [1.0_c_double]
    target_energy = [10.0_c_double]
    edges(1)%fast_index = 0_c_int
    edges(1)%target_index = 0_c_int
    edges(1)%reactivity_m3_s = 1.0_c_double
    edges(1)%target_energy_reactivity_J_m3_s = 0.0_c_double
    trial_fast = -1.0_c_double
    trial_target = -1.0_c_double
    trial_target_energy = -1.0_c_double
    edge_events = -1.0_c_double

    call target_network_trial(1.0_c_double, fast_energy, old_fast, &
         target_number, target_energy, edges, trial_fast, trial_target, &
         trial_target_energy, edge_events, out, status)

    call check(status == PB11_STATUS_OK, 'golden root call succeeds')
    root = (sqrt(5.0_c_double) - 1.0_c_double) / 2.0_c_double
    reaction = 1.0_c_double - root
    call check_close(trial_fast(1), root, 'golden fast inventory')
    call check_close(trial_target(1), root, 'golden target inventory')
    call check_close(trial_target_energy(1), 10.0_c_double, &
         'golden zero-M target energy')
    call check_close(edge_events(1), reaction, 'golden edge reaction')
    call check_close(out%reactions_m3, reaction, 'golden total reaction')
    call check_close(out%removed_fast_energy_J_m3, 5.0_c_double * reaction, &
         'golden removed fast energy')
    call check_close(out%removed_target_energy_J_m3, 0.0_c_double, &
         'golden removed target energy')
  end subroutine test_golden_root

  subroutine test_zero_edges()
    real(c_double), target :: fast_energy(2), old_fast(2)
    real(c_double), target :: target_number(2), target_energy(2)
    real(c_double), target :: trial_fast(2), trial_target(2)
    real(c_double), target :: trial_target_energy(2)
    real(c_double), allocatable, target :: edge_events(:)
    type(fusion_target_network_edge_v1), allocatable, target :: edges(:)
    type(fusion_target_network_v1), target :: out
    integer(c_int) :: status

    allocate(edges(0), edge_events(0))
    fast_energy = [2.0_c_double, 3.0_c_double]
    old_fast = [4.0_c_double, 5.0_c_double]
    target_number = [6.0_c_double, 7.0_c_double]
    target_energy = [8.0_c_double, 9.0_c_double]
    trial_fast = -1.0_c_double
    trial_target = -1.0_c_double
    trial_target_energy = -1.0_c_double

    call target_network_trial(1.0_c_double, fast_energy, old_fast, &
         target_number, target_energy, edges, trial_fast, trial_target, &
         trial_target_energy, edge_events, out, status)

    call check(status == PB11_STATUS_OK, 'zero-edge call succeeds')
    call check_close(trial_fast(1), old_fast(1), 'zero-edge fast[1]')
    call check_close(trial_fast(2), old_fast(2), 'zero-edge fast[2]')
    call check_close(trial_target(1), target_number(1), 'zero-edge target[1]')
    call check_close(trial_target(2), target_number(2), 'zero-edge target[2]')
    call check_close(trial_target_energy(1), target_energy(1), &
         'zero-edge target energy[1]')
    call check_close(trial_target_energy(2), target_energy(2), &
         'zero-edge target energy[2]')
    call check_network_zero(out, 'zero-edge output')
    deallocate(edges, edge_events)
  end subroutine test_zero_edges

  subroutine test_shape_rejection()
    real(c_double), target :: fast_energy(2), old_fast(1)
    real(c_double), target :: target_number(1), target_energy(1)
    real(c_double), target :: trial_fast(1), trial_target(1)
    real(c_double), target :: trial_target_energy(1), edge_events(1)
    type(fusion_target_network_edge_v1), target :: edges(1)
    type(fusion_target_network_v1), target :: out
    integer(c_int) :: status

    fast_energy = [5.0_c_double, 6.0_c_double]
    old_fast = [1.0_c_double]
    target_number = [1.0_c_double]
    target_energy = [1.0_c_double]
    edges(1)%fast_index = 0_c_int
    edges(1)%target_index = 0_c_int
    edges(1)%reactivity_m3_s = 1.0_c_double
    edges(1)%target_energy_reactivity_J_m3_s = 0.0_c_double
    trial_fast = 77.0_c_double
    trial_target = 77.0_c_double
    trial_target_energy = 77.0_c_double
    edge_events = 77.0_c_double
    out%reactions_m3 = 77.0_c_double
    out%removed_fast_energy_J_m3 = 77.0_c_double
    out%removed_target_energy_J_m3 = 77.0_c_double
    out%max_fast_number_relative_residual = 77.0_c_double
    out%max_target_number_relative_residual = 77.0_c_double
    out%max_target_energy_relative_residual = 77.0_c_double
    out%iterations = 77_c_int

    ! The fast-energy extent (2) disagrees with the old/trial extent (1).
    ! The wrapper must reject this before taking any C_LOC and clear outputs.
    call target_network_trial(1.0_c_double, fast_energy, old_fast, &
         target_number, target_energy, edges, trial_fast, trial_target, &
         trial_target_energy, edge_events, out, status)

    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'shape mismatch returns invalid-argument status')
    call check_close(trial_fast(1), 0.0_c_double, &
         'shape mismatch clears fast output')
    call check_close(trial_target(1), 0.0_c_double, &
         'shape mismatch clears target output')
    call check_close(trial_target_energy(1), 0.0_c_double, &
         'shape mismatch clears target-energy output')
    call check_close(edge_events(1), 0.0_c_double, &
         'shape mismatch clears edge output')
    call check_network_zero(out, 'shape mismatch output')
  end subroutine test_shape_rejection

end program test_fusion_target_network_fortran
