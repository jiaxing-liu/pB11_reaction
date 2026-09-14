program test_fusion_radial_transport_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use fusion_radial_transport_fortran, only : &
       fusion_radial_ledger_v1, fusion_radial_transport_trial, &
       PB11_STATUS_OK, PB11_STATUS_INVALID_ARGUMENT
  implicit none

  call test_ledger_layout()
  call test_two_component_diffusion()
  call test_malformed_extents()
  call test_invalid_final_component_conductance()

  print *, 'PASS radial transport Fortran ABI, diffusion, signed faces, inventory and atomic failures'

contains

  subroutine test_ledger_layout()
    type(fusion_radial_ledger_v1) :: ledger(1)

    if (c_sizeof(ledger(1)) /= 40_c_size_t) stop 1
  end subroutine test_ledger_layout

  subroutine test_two_component_diffusion()
    integer(c_int), parameter :: zones = 2_c_int
    integer(c_int), parameter :: components = 2_c_int
    real(c_double) :: volume_old(zones), volume_new(zones)
    real(c_double) :: advection(zones + 1_c_int)
    real(c_double) :: conductance(zones + 1_c_int, components)
    real(c_double) :: old_density(zones, components)
    real(c_double) :: boundary_density(2, components)
    real(c_double) :: trial_density(zones, components)
    real(c_double) :: face_amount(zones + 1_c_int, components)
    type(fusion_radial_ledger_v1) :: ledger(components)
    integer(c_int) :: status

    volume_old = 1._c_double
    volume_new = 1._c_double
    advection = 0._c_double
    conductance = 0._c_double
    conductance(2, 1) = 2._c_double
    conductance(2, 2) = 1._c_double
    old_density(:, 1) = [10._c_double, 0._c_double]
    old_density(:, 2) = [0._c_double, 20._c_double]
    boundary_density = 0._c_double
    trial_density = -1._c_double
    face_amount = -1._c_double
    call set_ledger_sentinel(ledger, -1._c_double)

    call fusion_radial_transport_trial(zones, components, 0.1_c_double, &
         volume_old, volume_new, advection, conductance, old_density, &
         boundary_density, trial_density, face_amount, ledger, status)
    if (status /= PB11_STATUS_OK) stop 2

    call assert_close(trial_density(1, 1), 60._c_double / 7._c_double)
    call assert_close(trial_density(2, 1), 10._c_double / 7._c_double)
    call assert_close(trial_density(1, 2), 5._c_double / 3._c_double)
    call assert_close(trial_density(2, 2), 55._c_double / 3._c_double)

    ! Face amounts are signed outward amounts: component 1 flows outward,
    ! while component 2 flows inward across the same interior face.
    call assert_close(face_amount(1, 1), 0._c_double)
    call assert_close(face_amount(2, 1), 10._c_double / 7._c_double)
    call assert_close(face_amount(3, 1), 0._c_double)
    call assert_close(face_amount(1, 2), 0._c_double)
    call assert_close(face_amount(2, 2), -5._c_double / 3._c_double)
    call assert_close(face_amount(3, 2), 0._c_double)

    ! Closed boundaries preserve each component inventory independently.
    call assert_close(sum(volume_old * old_density(:, 1)), &
         ledger(1)%initial_number)
    call assert_close(sum(volume_new * trial_density(:, 1)), &
         ledger(1)%final_number)
    call assert_close(ledger(1)%initial_number, 10._c_double)
    call assert_close(ledger(1)%final_number, 10._c_double)
    call assert_close(ledger(1)%inner_inward_number, 0._c_double)
    call assert_close(ledger(1)%outer_outward_number, 0._c_double)
    call assert_close(ledger(1)%balance_error, 0._c_double)

    call assert_close(sum(volume_old * old_density(:, 2)), &
         ledger(2)%initial_number)
    call assert_close(sum(volume_new * trial_density(:, 2)), &
         ledger(2)%final_number)
    call assert_close(ledger(2)%initial_number, 20._c_double)
    call assert_close(ledger(2)%final_number, 20._c_double)
    call assert_close(ledger(2)%inner_inward_number, 0._c_double)
    call assert_close(ledger(2)%outer_outward_number, 0._c_double)
    call assert_close(ledger(2)%balance_error, 0._c_double)
  end subroutine test_two_component_diffusion

  subroutine test_malformed_extents()
    integer(c_int), parameter :: zones = 2_c_int
    integer(c_int), parameter :: components = 2_c_int
    real(c_double) :: volume_old(zones), volume_old_bad(1)
    real(c_double) :: volume_new(zones)
    real(c_double) :: advection(zones + 1_c_int)
    real(c_double) :: conductance(zones + 1_c_int, components)
    real(c_double) :: old_density(zones, components)
    real(c_double) :: boundary_density(2, components)
    real(c_double) :: trial_density(zones, components)
    real(c_double) :: trial_density_bad(1, components)
    real(c_double) :: face_amount(zones + 1_c_int, components)
    type(fusion_radial_ledger_v1) :: ledger(components)
    integer(c_int) :: status

    volume_old = 1._c_double
    volume_old_bad = 1._c_double
    volume_new = 1._c_double
    advection = 0._c_double
    conductance = 0._c_double
    old_density = 1._c_double
    boundary_density = 0._c_double

    ! A malformed input extent must return before any C address is formed,
    ! while every output remains deterministically cleared.
    trial_density = 17._c_double
    face_amount = 18._c_double
    call set_ledger_sentinel(ledger, 19._c_double)
    call fusion_radial_transport_trial(zones, components, 0.1_c_double, &
         volume_old_bad, volume_new, advection, conductance, old_density, &
         boundary_density, trial_density, face_amount, ledger, status)
    if (status /= PB11_STATUS_INVALID_ARGUMENT) stop 27
    call assert_outputs_clear(trial_density, face_amount, ledger)

    ! A malformed output extent must clear the supplied short output and all
    ! other output objects before returning.
    trial_density_bad = 27._c_double
    face_amount = 28._c_double
    call set_ledger_sentinel(ledger, 29._c_double)
    call fusion_radial_transport_trial(zones, components, 0.1_c_double, &
         volume_old, volume_new, advection, conductance, old_density, &
         boundary_density, trial_density_bad, face_amount, ledger, status)
    if (status /= PB11_STATUS_INVALID_ARGUMENT) stop 30
    call assert_outputs_clear(trial_density_bad, face_amount, ledger)
  end subroutine test_malformed_extents

  subroutine test_invalid_final_component_conductance()
    integer(c_int), parameter :: zones = 2_c_int
    integer(c_int), parameter :: components = 2_c_int
    real(c_double) :: volume_old(zones), volume_new(zones)
    real(c_double) :: advection(zones + 1_c_int)
    real(c_double) :: conductance(zones + 1_c_int, components)
    real(c_double) :: old_density(zones, components)
    real(c_double) :: boundary_density(2, components)
    real(c_double) :: trial_density(zones, components)
    real(c_double) :: face_amount(zones + 1_c_int, components)
    type(fusion_radial_ledger_v1) :: ledger(components)
    integer(c_int) :: status

    volume_old = 1._c_double
    volume_new = 1._c_double
    advection = 0._c_double
    conductance = 0._c_double
    conductance(2, 1) = 2._c_double
    conductance(2, 2) = -1._c_double
    old_density = 1._c_double
    boundary_density = 0._c_double
    trial_density = 37._c_double
    face_amount = 38._c_double
    call set_ledger_sentinel(ledger, 39._c_double)

    ! The bad conductance is in the final component.  No earlier component
    ! may be published when the shared call fails validation.
    call fusion_radial_transport_trial(zones, components, 0.1_c_double, &
         volume_old, volume_new, advection, conductance, old_density, &
         boundary_density, trial_density, face_amount, ledger, status)
    if (status /= PB11_STATUS_INVALID_ARGUMENT) stop 32
    call assert_outputs_clear(trial_density, face_amount, ledger)
  end subroutine test_invalid_final_component_conductance

  subroutine assert_close(actual, expected)
    real(c_double), intent(in) :: actual, expected
    real(c_double), parameter :: tolerance = 2.e-13_c_double
    real(c_double) :: scale

    scale = max(1._c_double, abs(actual), abs(expected))
    if (abs(actual - expected) > tolerance * scale) error stop 1
  end subroutine assert_close

  subroutine set_ledger_sentinel(ledger, value)
    type(fusion_radial_ledger_v1), intent(out) :: ledger(:)
    real(c_double), intent(in) :: value
    integer(c_int) :: i, count

    count = int(size(ledger), kind=c_int)
    do i = 1_c_int, count
       ledger(i)%initial_number = value
       ledger(i)%final_number = value
       ledger(i)%inner_inward_number = value
       ledger(i)%outer_outward_number = value
       ledger(i)%balance_error = value
    end do
  end subroutine set_ledger_sentinel

  subroutine assert_outputs_clear(trial_density, face_amount, ledger)
    real(c_double), intent(in) :: trial_density(:,:), face_amount(:,:)
    type(fusion_radial_ledger_v1), intent(in) :: ledger(:)
    integer(c_int) :: i, count

    if (any(trial_density /= 0._c_double)) error stop 1
    if (any(face_amount /= 0._c_double)) error stop 1
    count = int(size(ledger), kind=c_int)
    do i = 1_c_int, count
       if (ledger(i)%initial_number /= 0._c_double) error stop 1
       if (ledger(i)%final_number /= 0._c_double) error stop 1
       if (ledger(i)%inner_inward_number /= 0._c_double) error stop 1
       if (ledger(i)%outer_outward_number /= 0._c_double) error stop 1
       if (ledger(i)%balance_error /= 0._c_double) error stop 1
    end do
  end subroutine assert_outputs_clear

end program test_fusion_radial_transport_fortran
