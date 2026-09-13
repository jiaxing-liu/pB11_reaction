program test_fusion_laboratory_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, &
       ieee_value
  use fusion_laboratory_fortran
  implicit none

  integer :: failures

  failures = 0
  call verify_layout()
  call verify_two_body_cm()
  call verify_boost_particles()
  call verify_noncontiguous_actuals()
  call verify_errors()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran fusion-laboratory test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran fusion-laboratory tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function close_relative(actual, expected, tolerance)
    real(c_double), intent(in) :: actual, expected, tolerance
    real(c_double) :: scale

    if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
       close_relative = .false.
       return
    end if
    if (actual == expected) then
       close_relative = .true.
       return
    end if
    scale = max(abs(actual), abs(expected))
    close_relative = scale > 0.0_c_double .and. &
         abs(actual - expected) <= tolerance * scale
  end function close_relative

  logical function particles_are_zero(particles)
    type(fusion_particle_four_vector_v1), intent(in) :: particles(:)
    integer :: i

    particles_are_zero = .true.
    do i = 1, size(particles)
       particles_are_zero = particles_are_zero .and. &
            particles(i)%mass_kg == 0.0_c_double .and. &
            particles(i)%kinetic_energy_J == 0.0_c_double .and. &
            all(particles(i)%momentum_kg_m_s == 0.0_c_double)
    end do
  end function particles_are_zero

  logical function ledger_is_zero(ledger)
    type(fusion_boost_ledger_v1), intent(in) :: ledger

    ledger_is_zero = ledger%expected_kinetic_energy_J == 0.0_c_double .and. &
         ledger%output_kinetic_energy_J == 0.0_c_double .and. &
         ledger%energy_residual_J == 0.0_c_double .and. &
         ledger%momentum_residual_kg_m_s == 0.0_c_double
  end function ledger_is_zero

  subroutine fill_particle(particle, mass, kinetic, px, py, pz)
    type(fusion_particle_four_vector_v1), intent(out) :: particle
    real(c_double), intent(in) :: mass, kinetic, px, py, pz

    particle%mass_kg = mass
    particle%kinetic_energy_J = kinetic
    particle%momentum_kg_m_s = [px, py, pz]
  end subroutine fill_particle

  subroutine fill_particles(particles, value)
    type(fusion_particle_four_vector_v1), intent(out) :: particles(:)
    real(c_double), intent(in) :: value
    integer :: i

    do i = 1, size(particles)
       call fill_particle(particles(i), value, value, value, value, value)
    end do
  end subroutine fill_particles

  logical function particles_match(left, right, tolerance)
    type(fusion_particle_four_vector_v1), intent(in) :: left, right
    real(c_double), intent(in) :: tolerance

    particles_match = close_relative(left%mass_kg, right%mass_kg, tolerance) &
         .and. close_relative(left%kinetic_energy_J, right%kinetic_energy_J, &
         tolerance) .and. all(abs(left%momentum_kg_m_s - &
         right%momentum_kg_m_s) <= tolerance * &
         max(1.0e-30_c_double, maxval(abs(left%momentum_kg_m_s)), &
         maxval(abs(right%momentum_kg_m_s))))
  end function particles_match

  subroutine verify_layout()
    type(fusion_particle_four_vector_v1) :: particle
    type(fusion_boost_ledger_v1) :: ledger
    real(c_double) :: double_size

    double_size = real(c_sizeof(0.0_c_double), c_double)
    call check(c_sizeof(particle) == 5 * c_sizeof(0.0_c_double), &
         'particle C layout contains five doubles')
    call check(c_sizeof(ledger) == 4 * c_sizeof(0.0_c_double), &
         'boost ledger C layout contains four doubles')
    call check(double_size > 0.0_c_double, 'C double size is positive')
  end subroutine verify_layout

  subroutine verify_two_body_cm()
    real(c_double), parameter :: mass_a = 2.0e-27_c_double
    real(c_double), parameter :: mass_b = 3.0e-27_c_double
    real(c_double), parameter :: available = 1.0e-13_c_double
    real(c_double) :: direction(3), momentum_scale
    type(fusion_particle_four_vector_v1) :: output(2)
    integer(c_int) :: status

    direction = [0.6_c_double, 0.8_c_double, 0.0_c_double]
    call fill_particle(output(1), -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double)
    call fill_particle(output(2), -1.0_c_double, -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double)

    call fusion_two_body_cm(mass_a, mass_b, available, direction, output, &
         status)
    call check(status == PB11_STATUS_OK, 'two-body CM call returns OK')
    call check(output(1)%mass_kg == mass_a .and. &
         output(2)%mass_kg == mass_b, 'two-body CM preserves product masses')
    call check(all(ieee_is_finite(output(1)%momentum_kg_m_s)) .and. &
         all(ieee_is_finite(output(2)%momentum_kg_m_s)) .and. &
         all(ieee_is_finite(output%kinetic_energy_J)), &
         'two-body CM output is finite')
    call check(all(output%kinetic_energy_J >= 0.0_c_double), &
         'two-body CM kinetic energies are nonnegative')
    call check(close_relative(sum(output%kinetic_energy_J), available, &
         1.0e-10_c_double), 'two-body CM kinetic energy is conserved')

    momentum_scale = sum(abs(output(1)%momentum_kg_m_s)) + &
         sum(abs(output(2)%momentum_kg_m_s))
    call check(all(abs(output(1)%momentum_kg_m_s + &
         output(2)%momentum_kg_m_s) <= &
         1.0e-12_c_double * max(1.0e-30_c_double, momentum_scale)), &
         'two-body CM momentum is opposite')
    call check(sum(output(1)%momentum_kg_m_s * direction) > 0.0_c_double, &
         'two-body CM product 0 follows supplied direction')
  end subroutine verify_two_body_cm

  subroutine verify_boost_particles()
    real(c_double), parameter :: mass_a = 2.0e-27_c_double
    real(c_double), parameter :: mass_b = 3.0e-27_c_double
    real(c_double) :: velocity(3), input_energy, output_energy
    type(fusion_particle_four_vector_v1) :: input(2), output(2)
    type(fusion_boost_ledger_v1) :: ledger
    integer(c_int) :: status

    call fill_particle(input(1), mass_a, 0.0_c_double, 0.0_c_double, &
         0.0_c_double, 0.0_c_double)
    call fill_particle(input(2), mass_b, 0.0_c_double, 0.0_c_double, &
         0.0_c_double, 0.0_c_double)
    output = input
    velocity = [1.0e6_c_double, 2.0e5_c_double, 0.0_c_double]
    ledger = fusion_boost_ledger_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double)

    call fusion_boost_particles(2_c_int, velocity, input, output, ledger, &
         status)
    call check(status == PB11_STATUS_OK, 'Lorentz boost call returns OK')
    call check(input(1)%kinetic_energy_J == 0.0_c_double .and. &
         input(2)%kinetic_energy_J == 0.0_c_double .and. &
         all(input(1)%momentum_kg_m_s == 0.0_c_double) .and. &
         all(input(2)%momentum_kg_m_s == 0.0_c_double), &
         'Lorentz boost leaves input particles unchanged')
    call check(output(1)%mass_kg == mass_a .and. &
         output(2)%mass_kg == mass_b, 'Lorentz boost preserves masses')
    call check(all(output%kinetic_energy_J > 0.0_c_double) .and. &
         all(ieee_is_finite(output%kinetic_energy_J)) .and. &
         all(ieee_is_finite(output(1)%momentum_kg_m_s)) .and. &
         all(ieee_is_finite(output(2)%momentum_kg_m_s)), &
         'boosted rest particles have finite positive kinetic energy')
    call check(output(1)%momentum_kg_m_s(1) > 0.0_c_double .and. &
         output(1)%momentum_kg_m_s(2) > 0.0_c_double, &
         'boosted particles move with the supplied velocity')

    input_energy = sum(input%kinetic_energy_J)
    output_energy = sum(output%kinetic_energy_J)
    call check(close_relative(ledger%expected_kinetic_energy_J, output_energy, &
         1.0e-12_c_double), 'boost ledger expected energy matches output')
    call check(close_relative(ledger%output_kinetic_energy_J, output_energy, &
         1.0e-12_c_double), 'boost ledger reports output energy')
    call check(abs(ledger%energy_residual_J) <= &
         1.0e-10_c_double * max(1.0e-30_c_double, output_energy) .and. &
         abs(ledger%momentum_residual_kg_m_s) <= 1.0e-10_c_double * &
         max(1.0e-30_c_double, sum(abs(output(1)%momentum_kg_m_s)) + &
         sum(abs(output(2)%momentum_kg_m_s))), &
         'boost ledger residuals are small')
    call check(input_energy == 0.0_c_double, &
         'boost input energy remains zero after call')
  end subroutine verify_boost_particles

  subroutine verify_noncontiguous_actuals()
    real(c_double) :: direction_store(5), direction_ref(3)
    type(fusion_particle_four_vector_v1) :: output_store(4), output_ref(2)
    integer(c_int) :: status, status_ref

    direction_store = [0.6_c_double, -1.0_c_double, 0.8_c_double, &
         -1.0_c_double, 0.0_c_double]
    direction_ref = [0.6_c_double, 0.8_c_double, 0.0_c_double]
    call fill_particles(output_store, -1.0_c_double)
    call fill_particle(output_ref(1), -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double)
    call fill_particle(output_ref(2), -1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double, -1.0_c_double)

    call fusion_two_body_cm(2.0e-27_c_double, 3.0e-27_c_double, &
         1.0e-13_c_double, direction_store(1:5:2), output_store(1:4:2), &
         status)
    call fusion_two_body_cm(2.0e-27_c_double, 3.0e-27_c_double, &
         1.0e-13_c_double, direction_ref, output_ref, status_ref)
    call check(status == PB11_STATUS_OK .and. status_ref == PB11_STATUS_OK, &
         'noncontiguous actuals return OK')
    call check(particles_match(output_store(1), output_ref(1), &
         1.0e-12_c_double) .and. &
         particles_match(output_store(3), output_ref(2), 1.0e-12_c_double), &
         'noncontiguous direction and output match contiguous reference')
  end subroutine verify_noncontiguous_actuals

  subroutine verify_errors()
    real(c_double) :: direction(3), bad_direction(2), nan_direction(3)
    real(c_double) :: velocity(3), superluminal(3)
    type(fusion_particle_four_vector_v1) :: output(2), input(2)
    type(fusion_particle_four_vector_v1) :: short_input(1)
    type(fusion_particle_four_vector_v1) :: empty_input(0), empty_output(0)
    type(fusion_boost_ledger_v1) :: ledger
    integer(c_int) :: status

    direction = [1.0_c_double, 0.0_c_double, 0.0_c_double]
    bad_direction = [1.0_c_double, 0.0_c_double]
    nan_direction = direction
    nan_direction(2) = ieee_value(0.0_c_double, ieee_quiet_nan)
    call fill_particles(output, -1.0_c_double)

    call fusion_two_body_cm(2.0e-27_c_double, 3.0e-27_c_double, &
         1.0e-13_c_double, bad_direction, output, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         particles_are_zero(output), 'direction extent error clears output')

    call fill_particles(output, -1.0_c_double)
    call fusion_two_body_cm(-2.0e-27_c_double, 3.0e-27_c_double, &
         1.0e-13_c_double, direction, output, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         particles_are_zero(output), 'negative mass clears two-body output')

    call fill_particles(output, -1.0_c_double)
    call fusion_two_body_cm(2.0e-27_c_double, 3.0e-27_c_double, &
         1.0e-13_c_double, nan_direction, output, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         particles_are_zero(output), 'NaN direction clears two-body output')

    call fill_particle(input(1), 2.0e-27_c_double, 0.0_c_double, 0.0_c_double, &
         0.0_c_double, 0.0_c_double)
    call fill_particle(input(2), 3.0e-27_c_double, 0.0_c_double, 0.0_c_double, &
         0.0_c_double, 0.0_c_double)
    velocity = [1.0e6_c_double, 0.0_c_double, 0.0_c_double]
    superluminal = [3.0e8_c_double, 0.0_c_double, 0.0_c_double]
    call fill_particles(output, -1.0_c_double)
    ledger = fusion_boost_ledger_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double)
    call fusion_boost_particles(2_c_int, superluminal, input, output, ledger, &
         status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         particles_are_zero(output) .and. ledger_is_zero(ledger), &
         'superluminal boost clears output and ledger')

    call fill_particles(output, -1.0_c_double)
    ledger = fusion_boost_ledger_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double)
    call fusion_boost_particles(2_c_int, velocity, short_input, output, &
         ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         particles_are_zero(output) .and. ledger_is_zero(ledger), &
         'boost input extent error clears output and ledger')

    ledger = fusion_boost_ledger_v1(-1.0_c_double, -1.0_c_double, &
         -1.0_c_double, -1.0_c_double)
    call fusion_boost_particles(0_c_int, velocity, empty_input, empty_output, &
         ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         particles_are_zero(empty_output) .and. ledger_is_zero(ledger), &
         'zero boost count is rejected and cleared')
  end subroutine verify_errors

end program test_fusion_laboratory_fortran
