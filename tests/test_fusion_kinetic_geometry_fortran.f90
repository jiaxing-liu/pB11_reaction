program test_fusion_kinetic_geometry_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use fusion_kinetic_geometry_fortran, only : &
       fusion_kinetic_geometry_trial, fusion_transport_ledger_v1, &
       PB11_STATUS_OK, PB11_STATUS_INVALID_ARGUMENT
  implicit none

  integer :: failures

  failures = 0
  call test_shape_and_layout()
  call test_identity_with_nonfirst_species()
  call test_volume_dilution_with_s_and_t()
  call test_invalid_shape_clears_outputs()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, &
          ' Fortran kinetic-geometry binding test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran kinetic-geometry binding tests passed'

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

  subroutine clear_transport(ledger)
    type(fusion_transport_ledger_v1), intent(out) :: ledger(:)
    integer :: i

    do i = 1, size(ledger)
       ledger(i)%spatial_number = 0.0_c_double
       ledger(i)%spatial_energy_J = 0.0_c_double
       ledger(i)%work_J = 0.0_c_double
       ledger(i)%lower_number = 0.0_c_double
       ledger(i)%lower_energy_J = 0.0_c_double
       ledger(i)%upper_number = 0.0_c_double
       ledger(i)%upper_energy_J = 0.0_c_double
    end do
  end subroutine clear_transport

  logical function transport_is_zero(ledger)
    type(fusion_transport_ledger_v1), intent(in) :: ledger(:)
    integer :: i

    transport_is_zero = .true.
    do i = 1, size(ledger)
       transport_is_zero = transport_is_zero .and. &
            all(ledger(i)%spatial_number == 0.0_c_double) .and. &
            all(ledger(i)%spatial_energy_J == 0.0_c_double) .and. &
            all(ledger(i)%work_J == 0.0_c_double) .and. &
            all(ledger(i)%lower_number == 0.0_c_double) .and. &
            all(ledger(i)%lower_energy_J == 0.0_c_double) .and. &
            all(ledger(i)%upper_number == 0.0_c_double) .and. &
            all(ledger(i)%upper_energy_J == 0.0_c_double)
    end do
  end function transport_is_zero

  subroutine set_common_inputs(edges, volume_old, volume_new, advection, &
       compression, conductance, boundary_s, boundary_t)
    real(c_double), intent(out) :: edges(:), volume_old(:), volume_new(:)
    real(c_double), intent(out) :: advection(:), compression(:)
    real(c_double), intent(out) :: conductance(:,:,:)
    real(c_double), intent(out) :: boundary_s(:,:,:), boundary_t(:,:,:)

    edges = [0.0_c_double, 1.0_c_double, 2.0_c_double]
    volume_old = 1.0_c_double
    volume_new = 1.0_c_double
    advection = 0.0_c_double
    compression = 0.0_c_double
    conductance = 0.0_c_double
    boundary_s = 0.0_c_double
    boundary_t = 0.0_c_double
  end subroutine set_common_inputs

  subroutine test_shape_and_layout()
    integer(c_int), parameter :: zones = 2_c_int, cells = 2_c_int
    real(c_double) :: edges(cells + 1_c_int)
    real(c_double) :: volume_old(zones), volume_new(zones)
    real(c_double) :: advection(zones + 1_c_int), compression(zones)
    real(c_double) :: conductance(zones + 1_c_int, cells, 6)
    real(c_double) :: boundary_s(2, cells, 6), boundary_t(2, cells, 6)
    real(c_double) :: old_s(cells, 6, zones), old_t(cells, 6, zones)
    real(c_double) :: trial_s(cells, 6, zones), trial_t(cells, 6, zones)
    type(fusion_transport_ledger_v1) :: ledger(zones)
    integer(c_int) :: status

    call check(c_sizeof(ledger(1)) == 42_c_size_t * 8_c_size_t, &
         'transport ledger has the C 42-double layout')
    call check(size(conductance) == (zones + 1_c_int) * cells * 6_c_int, &
         'conductance uses face, cell, species shape')
    call check(size(old_s) == cells * 6_c_int * zones, &
         'density uses cell, species, zone shape')

    call set_common_inputs(edges, volume_old, volume_new, advection, &
         compression, conductance, boundary_s, boundary_t)
    old_s = 0.0_c_double
    old_t = 0.0_c_double
    old_s(2, 4, 2) = 7.25_c_double
    old_t(1, 6, 1) = 3.5_c_double
    trial_s = -1.0_c_double
    trial_t = -1.0_c_double
    call clear_transport(ledger)

    call fusion_kinetic_geometry_trial(zones, cells, 0.25_c_double, edges, &
         volume_old, volume_new, advection, compression, conductance, &
         boundary_s, boundary_t, old_s, old_t, trial_s, trial_t, ledger, &
         status)
    call check_status(status, 'shape/layout trial succeeds')
    call check(all(trial_s == old_s) .and. all(trial_t == old_t), &
         'shape/layout trial preserves identity populations')
  end subroutine test_shape_and_layout

  subroutine test_identity_with_nonfirst_species()
    integer(c_int), parameter :: zones = 2_c_int, cells = 2_c_int
    real(c_double) :: edges(cells + 1_c_int)
    real(c_double) :: volume_old(zones), volume_new(zones)
    real(c_double) :: advection(zones + 1_c_int), compression(zones)
    real(c_double) :: conductance(zones + 1_c_int, cells, 6)
    real(c_double) :: boundary_s(2, cells, 6), boundary_t(2, cells, 6)
    real(c_double) :: old_s(cells, 6, zones), old_t(cells, 6, zones)
    real(c_double) :: trial_s(cells, 6, zones), trial_t(cells, 6, zones)
    type(fusion_transport_ledger_v1) :: ledger(zones)
    integer(c_int) :: status

    call set_common_inputs(edges, volume_old, volume_new, advection, &
         compression, conductance, boundary_s, boundary_t)
    old_s = 0.0_c_double
    old_t = 0.0_c_double
    old_s(1, 2, 1) = 2.0_c_double
    old_s(2, 5, 2) = 11.0_c_double
    old_t(2, 3, 1) = 4.0_c_double
    old_t(1, 6, 2) = 13.0_c_double
    trial_s = -9.0_c_double
    trial_t = -9.0_c_double
    call clear_transport(ledger)

    call fusion_kinetic_geometry_trial(zones, cells, 0.5_c_double, edges, &
         volume_old, volume_new, advection, compression, conductance, &
         boundary_s, boundary_t, old_s, old_t, trial_s, trial_t, ledger, &
         status)
    call check_status(status, 'identity trial succeeds')
    call check(all(trial_s == old_s) .and. all(trial_t == old_t), &
         'zero transport/work preserves all S/T species and cells')
    call check(transport_is_zero(ledger), 'identity trial has zero transport ledger')
  end subroutine test_identity_with_nonfirst_species

  subroutine test_volume_dilution_with_s_and_t()
    integer(c_int), parameter :: zones = 2_c_int, cells = 2_c_int
    real(c_double) :: edges(cells + 1_c_int)
    real(c_double) :: volume_old(zones), volume_new(zones)
    real(c_double) :: advection(zones + 1_c_int), compression(zones)
    real(c_double) :: conductance(zones + 1_c_int, cells, 6)
    real(c_double) :: boundary_s(2, cells, 6), boundary_t(2, cells, 6)
    real(c_double) :: old_s(cells, 6, zones), old_t(cells, 6, zones)
    real(c_double) :: trial_s(cells, 6, zones), trial_t(cells, 6, zones)
    type(fusion_transport_ledger_v1) :: ledger(zones)
    integer(c_int) :: status

    call set_common_inputs(edges, volume_old, volume_new, advection, &
         compression, conductance, boundary_s, boundary_t)
    volume_old = 2.0_c_double
    volume_new = 1.0_c_double
    old_s = 0.0_c_double
    old_t = 0.0_c_double
    old_s(1, 3, 1) = 1.5_c_double
    old_s(2, 6, 2) = 2.25_c_double
    old_t(2, 2, 1) = 3.75_c_double
    old_t(1, 5, 2) = 4.5_c_double
    trial_s = -1.0_c_double
    trial_t = -1.0_c_double
    call clear_transport(ledger)

    call fusion_kinetic_geometry_trial(zones, cells, 0.25_c_double, edges, &
         volume_old, volume_new, advection, compression, conductance, &
         boundary_s, boundary_t, old_s, old_t, trial_s, trial_t, ledger, &
         status)
    call check_status(status, 'volume-dilution trial succeeds')
    call check(all(trial_s == 2.0_c_double * old_s) .and. &
         all(trial_t == 2.0_c_double * old_t), &
         'volume dilution preserves S/T inventory while doubling density')
    call check(transport_is_zero(ledger), &
         'volume dilution has no spatial or energy-domain ledger')
  end subroutine test_volume_dilution_with_s_and_t

  subroutine test_invalid_shape_clears_outputs()
    integer(c_int), parameter :: zones = 2_c_int, cells = 2_c_int
    real(c_double) :: edges(cells + 1_c_int)
    real(c_double) :: volume_old(zones), volume_new(zones)
    real(c_double) :: advection(zones + 1_c_int), compression(zones)
    real(c_double) :: conductance_bad(zones, cells, 6)
    real(c_double) :: conductance(zones + 1_c_int, cells, 6)
    real(c_double) :: boundary_s(2, cells, 6), boundary_t(2, cells, 6)
    real(c_double) :: old_s(cells, 6, zones), old_t(cells, 6, zones)
    real(c_double) :: trial_s(cells, 6, zones), trial_t(cells, 6, zones)
    type(fusion_transport_ledger_v1) :: ledger(zones)
    integer(c_int) :: status

    call set_common_inputs(edges, volume_old, volume_new, advection, &
         compression, conductance, boundary_s, boundary_t)
    conductance_bad = 1.0_c_double
    old_s = 1.0_c_double
    old_t = 2.0_c_double
    trial_s = 9.0_c_double
    trial_t = 9.0_c_double
    call clear_transport(ledger)
    ledger(1)%work_J = 7.0_c_double
    ledger(2)%work_J = 7.0_c_double

    call fusion_kinetic_geometry_trial(zones, cells, 0.25_c_double, edges, &
         volume_old, volume_new, advection, compression, conductance_bad, &
         boundary_s, boundary_t, old_s, old_t, trial_s, trial_t, ledger, &
         status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'wrong conductance extent is rejected before C_LOC')
    call check(all(trial_s == 0.0_c_double) .and. &
         all(trial_t == 0.0_c_double) .and. transport_is_zero(ledger), &
         'invalid shape clears all outputs and ledger fields')
  end subroutine test_invalid_shape_clears_outputs

end program test_fusion_kinetic_geometry_fortran
