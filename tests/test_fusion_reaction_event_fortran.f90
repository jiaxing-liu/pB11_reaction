program test_fusion_reaction_event_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_laboratory_fortran, only : fusion_particle_four_vector_v1, &
       fusion_boost_ledger_v1
  use fusion_reaction_event_fortran
  implicit none

  integer :: failures

  failures = 0
  call verify_layout()
  call verify_dt_zero_momentum()
  call verify_pb_lab_event()
  call verify_errors_and_clearing()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran reaction-event test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran reaction-event tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function near(actual, expected, relative_tolerance, absolute_tolerance)
    real(c_double), intent(in) :: actual, expected
    real(c_double), intent(in) :: relative_tolerance, absolute_tolerance

    near = ieee_is_finite(actual) .and. ieee_is_finite(expected) .and. &
         abs(actual - expected) <= absolute_tolerance + relative_tolerance * &
         max(abs(actual), abs(expected))
  end function near

  subroutine poison_parent(parent)
    type(fusion_reaction_parent_v1), intent(out) :: parent

    parent%reactant_kinetic_J = -1.0_c_double
    parent%classical_kinetic_J = -1.0_c_double
    parent%on_shell_kinetic_J = -1.0_c_double
    parent%relative_classical_energy_J = -1.0_c_double
    parent%available_cm_energy_J = -1.0_c_double
    parent%boost_velocity_m_s = -1.0_c_double
    parent%expected_product_lab_kinetic_J = -1.0_c_double
    parent%classical_minus_on_shell_J = -1.0_c_double
    parent%momentum_sum_kg_m_s = -1.0_c_double
    parent%product_count = -1_c_int
    parent%convention = -1_c_int
  end subroutine poison_parent

  subroutine poison_particles(output)
    type(fusion_particle_four_vector_v1), intent(out) :: output(:)
    integer :: i

    do i = 1, size(output)
       output(i)%mass_kg = -1.0_c_double
       output(i)%kinetic_energy_J = -1.0_c_double
       output(i)%momentum_kg_m_s = -1.0_c_double
    end do
  end subroutine poison_particles

  subroutine poison_ledger(ledger)
    type(fusion_boost_ledger_v1), intent(out) :: ledger

    ledger%expected_kinetic_energy_J = -1.0_c_double
    ledger%output_kinetic_energy_J = -1.0_c_double
    ledger%energy_residual_J = -1.0_c_double
    ledger%momentum_residual_kg_m_s = -1.0_c_double
  end subroutine poison_ledger

  logical function parent_is_zero(parent)
    type(fusion_reaction_parent_v1), intent(in) :: parent

    parent_is_zero = all(parent%reactant_kinetic_J == 0.0_c_double) .and. &
         all(parent%classical_kinetic_J == 0.0_c_double) .and. &
         all(parent%on_shell_kinetic_J == 0.0_c_double) .and. &
         parent%relative_classical_energy_J == 0.0_c_double .and. &
         parent%available_cm_energy_J == 0.0_c_double .and. &
         all(parent%boost_velocity_m_s == 0.0_c_double) .and. &
         parent%expected_product_lab_kinetic_J == 0.0_c_double .and. &
         parent%classical_minus_on_shell_J == 0.0_c_double .and. &
         all(parent%momentum_sum_kg_m_s == 0.0_c_double) .and. &
         parent%product_count == 0_c_int .and. parent%convention == 0_c_int
  end function parent_is_zero

  logical function particles_are_zero(output)
    type(fusion_particle_four_vector_v1), intent(in) :: output(:)
    integer :: i

    particles_are_zero = .true.
    do i = 1, size(output)
       particles_are_zero = particles_are_zero .and. &
            output(i)%mass_kg == 0.0_c_double .and. &
            output(i)%kinetic_energy_J == 0.0_c_double .and. &
            all(output(i)%momentum_kg_m_s == 0.0_c_double)
    end do
  end function particles_are_zero

  logical function ledger_is_zero(ledger)
    type(fusion_boost_ledger_v1), intent(in) :: ledger

    ledger_is_zero = ledger%expected_kinetic_energy_J == 0.0_c_double .and. &
         ledger%output_kinetic_energy_J == 0.0_c_double .and. &
         ledger%energy_residual_J == 0.0_c_double .and. &
         ledger%momentum_residual_kg_m_s == 0.0_c_double
  end function ledger_is_zero

  subroutine sum_momenta(particles, count, result)
    type(fusion_particle_four_vector_v1), intent(in) :: particles(:)
    integer, intent(in) :: count
    real(c_double), intent(out) :: result(3)
    integer :: i

    result = 0.0_c_double
    do i = 1, count
       result = result + particles(i)%momentum_kg_m_s
    end do
  end subroutine sum_momenta

  subroutine verify_layout()
    type(fusion_reaction_parent_v1) :: parent
    type(fusion_particle_four_vector_v1) :: particle
    type(fusion_boost_ledger_v1) :: ledger

    call check(c_sizeof(parent) == 16_c_size_t * c_sizeof(0.0_c_double) + &
         2_c_size_t * c_sizeof(0_c_int), &
         'parent layout is sixteen doubles plus two ints')
    call check(c_sizeof(particle) == 5_c_size_t * c_sizeof(0.0_c_double), &
         'particle layout is reused from laboratory module')
    call check(c_sizeof(ledger) == 4_c_size_t * c_sizeof(0.0_c_double), &
         'boost ledger layout is reused from laboratory module')
  end subroutine verify_layout

  subroutine verify_dt_zero_momentum()
    real(c_double) :: pa(3), pb(3), direction(3), product_energy
    real(c_double) :: product_momentum(3), momentum_scale, energy_scale
    type(fusion_reaction_parent_v1) :: parent
    type(fusion_particle_four_vector_v1) :: output(3)
    type(fusion_boost_ledger_v1) :: ledger
    integer(c_int) :: status

    pa = 0.0_c_double
    pb = 0.0_c_double
    direction = [0.0_c_double, 0.0_c_double, 1.0_c_double]

    call fusion_reaction_parent(FUSION_DT_ALPHAN, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, parent, status)
    call check(status == PB11_STATUS_OK, 'DT zero-momentum parent returns OK')
    call check(parent%product_count == 2_c_int .and. &
         parent%convention == FUSION_REACTANT_CLASSICAL_BUDGET, &
         'DT parent records channel size and classical option')
    call check(all(parent%reactant_kinetic_J == 0.0_c_double) .and. &
         all(parent%classical_kinetic_J == 0.0_c_double) .and. &
         all(parent%on_shell_kinetic_J == 0.0_c_double) .and. &
         parent%relative_classical_energy_J == 0.0_c_double .and. &
         all(parent%momentum_sum_kg_m_s == 0.0_c_double), &
         'DT zero-momentum kinetic and momentum fields are zero')
    call check(parent%available_cm_energy_J > 0.0_c_double .and. &
         parent%expected_product_lab_kinetic_J > 0.0_c_double .and. &
         all(ieee_is_finite(parent%boost_velocity_m_s)) .and. &
         all(parent%boost_velocity_m_s == 0.0_c_double), &
         'DT zero-momentum parent has finite Q energy and zero boost')

    call poison_particles(output)
    call poison_parent(parent)
    call poison_ledger(ledger)
    call fusion_reaction_lab_event(FUSION_DT_ALPHAN, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, direction, &
         0.0_c_double, 0.0_c_double, 0.0_c_double, output, parent, ledger, &
         status)
    call check(status == PB11_STATUS_OK, 'DT zero-momentum event returns OK')
    call check(parent%product_count == 2_c_int .and. &
         parent%convention == FUSION_REACTANT_CLASSICAL_BUDGET, &
         'DT event preserves parent metadata')
    call check(output(1)%mass_kg > 0.0_c_double .and. &
         output(2)%mass_kg > 0.0_c_double .and. output(3)%mass_kg == 0.0_c_double .and. &
         ieee_is_finite(output(1)%kinetic_energy_J) .and. &
         ieee_is_finite(output(2)%kinetic_energy_J), &
         'DT event writes two finite products and clears the third slot')
    product_energy = output(1)%kinetic_energy_J + output(2)%kinetic_energy_J
    call sum_momenta(output, 2, product_momentum)
    momentum_scale = max(1.0e-30_c_double, sum(abs(output(1)%momentum_kg_m_s)) + &
         sum(abs(output(2)%momentum_kg_m_s)))
    energy_scale = max(1.0e-30_c_double, product_energy)
    call check(near(product_energy, parent%expected_product_lab_kinetic_J, &
         5.0e-10_c_double, 1.0e-30_c_double), &
         'DT event product kinetic energy matches parent budget')
    call check(maxval(abs(product_momentum)) <= 5.0e-12_c_double * momentum_scale, &
         'DT zero-momentum products have zero total momentum')
    call check(near(ledger%expected_kinetic_energy_J, product_energy, &
         5.0e-10_c_double, 1.0e-30_c_double) .and. &
         near(ledger%output_kinetic_energy_J, product_energy, &
         5.0e-10_c_double, 1.0e-30_c_double) .and. &
         abs(ledger%energy_residual_J) <= 5.0e-10_c_double * energy_scale .and. &
         abs(ledger%momentum_residual_kg_m_s) <= 5.0e-10_c_double * momentum_scale, &
         'DT event ledger closes energy and momentum')
  end subroutine verify_dt_zero_momentum

  subroutine verify_pb_lab_event()
    real(c_double), parameter :: kev = 1.602176634e-16_c_double
    real(c_double), parameter :: speed_of_light = 299792458.0_c_double
    real(c_double) :: pa(3), pb(3), direction(3), q, cosine, azimuth
    real(c_double) :: product_energy, product_momentum(3), input_momentum(3)
    real(c_double) :: momentum_scale, energy_scale, boost_speed
    type(fusion_reaction_parent_v1) :: parent
    type(fusion_particle_four_vector_v1) :: output(3)
    type(fusion_boost_ledger_v1) :: ledger
    integer(c_int) :: status
    integer :: i

    ! Small but nonzero p/B11 momenta exercise the parent CM boost while
    ! retaining the requested 91.84-keV sequential pB event parameter.
    pa = [2.0e-20_c_double, -1.0e-20_c_double, 0.5e-20_c_double]
    pb = [-0.5e-20_c_double, 0.75e-20_c_double, -1.25e-20_c_double]
    direction = [0.6_c_double, 0.8_c_double, 0.0_c_double]
    q = 91.84_c_double * kev
    cosine = 0.37_c_double
    azimuth = 0.83_c_double

    call fusion_reaction_parent(FUSION_PB11_3ALPHA, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, parent, status)
    call check(status == PB11_STATUS_OK .and. parent%product_count == 3_c_int, &
         'pB nonzero-momentum parent returns a three-product event')
    call check(all(ieee_is_finite(parent%reactant_kinetic_J)) .and. &
         all(parent%reactant_kinetic_J >= 0.0_c_double) .and. &
         all(parent%classical_kinetic_J >= parent%on_shell_kinetic_J) .and. &
         parent%relative_classical_energy_J >= 0.0_c_double .and. &
         parent%available_cm_energy_J > 0.0_c_double .and. &
         parent%expected_product_lab_kinetic_J > 0.0_c_double, &
         'pB parent reports finite nonnegative energy budgets')
    call check(parent%classical_minus_on_shell_J >= 0.0_c_double, &
         'pB parent reports a nonnegative convention difference')
    boost_speed = sqrt(sum(parent%boost_velocity_m_s * &
         parent%boost_velocity_m_s))
    call check(ieee_is_finite(boost_speed) .and. boost_speed < speed_of_light, &
         'pB parent boost is subluminal')
    input_momentum = pa + pb
    call check(maxval(abs(parent%momentum_sum_kg_m_s - input_momentum)) <= &
         1.0e-12_c_double * max(1.0e-30_c_double, sum(abs(input_momentum))), &
         'pB parent records the input momentum sum')

    call poison_particles(output)
    call poison_parent(parent)
    call poison_ledger(ledger)
    call fusion_reaction_lab_event(FUSION_PB11_3ALPHA, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, direction, q, cosine, &
         azimuth, output, parent, ledger, status)
    call check(status == PB11_STATUS_OK, 'pB 91.84-keV event returns OK')
    call check(parent%product_count == 3_c_int .and. &
         parent%convention == FUSION_REACTANT_CLASSICAL_BUDGET, &
         'pB event preserves three-product classical metadata')
    do i = 1, 3
       call check(output(i)%mass_kg > 0.0_c_double .and. &
            output(i)%kinetic_energy_J >= 0.0_c_double .and. &
            ieee_is_finite(output(i)%mass_kg) .and. &
            ieee_is_finite(output(i)%kinetic_energy_J) .and. &
            all(ieee_is_finite(output(i)%momentum_kg_m_s)), &
            'pB product is finite and has a positive mass shell')
    end do
    product_energy = sum(output%kinetic_energy_J)
    call sum_momenta(output, 3, product_momentum)
    input_momentum = pa + pb
    momentum_scale = max(1.0e-30_c_double, sum(abs(product_momentum)) + &
         sum(abs(input_momentum)))
    energy_scale = max(1.0e-30_c_double, product_energy + &
         parent%expected_product_lab_kinetic_J)
    call check(near(product_energy, parent%expected_product_lab_kinetic_J, &
         5.0e-10_c_double, 1.0e-30_c_double), &
         'pB event product energy matches selected reactant budget plus Q')
    call check(maxval(abs(product_momentum - input_momentum)) <= &
         5.0e-11_c_double * momentum_scale, &
         'pB event product momentum matches the incoming pair')
    call check(near(ledger%expected_kinetic_energy_J, product_energy, &
         5.0e-10_c_double, 1.0e-30_c_double) .and. &
         near(ledger%output_kinetic_energy_J, product_energy, &
         5.0e-10_c_double, 1.0e-30_c_double) .and. &
         abs(ledger%energy_residual_J) <= 5.0e-10_c_double * energy_scale .and. &
         abs(ledger%momentum_residual_kg_m_s) <= 5.0e-10_c_double * momentum_scale, &
         'pB event ledger closes energy and momentum')
  end subroutine verify_pb_lab_event

  subroutine verify_errors_and_clearing()
    real(c_double) :: pa(3), pb(3), short_pa(2), direction(3), bad_direction(3)
    real(c_double) :: short_direction(2)
    type(fusion_reaction_parent_v1) :: parent
    type(fusion_particle_four_vector_v1) :: output(3), short_output(2)
    type(fusion_boost_ledger_v1) :: ledger
    integer(c_int) :: status

    pa = 0.0_c_double
    pb = 0.0_c_double
    short_pa = 0.0_c_double
    direction = [0.0_c_double, 0.0_c_double, 1.0_c_double]
    bad_direction = [1.0_c_double, 1.0_c_double, 1.0_c_double]
    short_direction = 0.0_c_double

    call poison_parent(parent)
    call fusion_reaction_parent(FUSION_DT_ALPHAN, &
         FUSION_REACTANT_CLASSICAL_BUDGET, short_pa, pb, parent, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. parent_is_zero(parent), &
         'parent rejects a short reactant vector and clears output')

    call poison_parent(parent)
    call fusion_reaction_parent(FUSION_DT_ALPHAN, 9_c_int, pa, pb, parent, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. parent_is_zero(parent), &
         'parent rejects an invalid energy convention and clears output')

    call poison_particles(output)
    call poison_parent(parent)
    call poison_ledger(ledger)
    call fusion_reaction_lab_event(FUSION_DT_ALPHAN, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, short_direction, &
         0.0_c_double, 0.0_c_double, 0.0_c_double, output, parent, ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         particles_are_zero(output) .and. parent_is_zero(parent) .and. &
         ledger_is_zero(ledger), &
         'event rejects a short direction before C_LOC and clears all outputs')

    call poison_particles(short_output)
    call poison_parent(parent)
    call poison_ledger(ledger)
    call fusion_reaction_lab_event(FUSION_DT_ALPHAN, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, direction, &
         0.0_c_double, 0.0_c_double, 0.0_c_double, short_output, parent, ledger, &
         status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         particles_are_zero(short_output) .and. parent_is_zero(parent) .and. &
         ledger_is_zero(ledger), &
         'event rejects a short output array before C_LOC and clears all outputs')

    call poison_particles(output)
    call poison_parent(parent)
    call poison_ledger(ledger)
    call fusion_reaction_lab_event(FUSION_DT_ALPHAN, 9_c_int, pa, pb, direction, &
         0.0_c_double, 0.0_c_double, 0.0_c_double, output, parent, ledger, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         particles_are_zero(output) .and. parent_is_zero(parent) .and. &
         ledger_is_zero(ledger), &
         'event rejects an invalid convention and clears all outputs')

    call poison_particles(output)
    call poison_parent(parent)
    call poison_ledger(ledger)
    call fusion_reaction_lab_event(FUSION_DT_ALPHAN, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, bad_direction, &
         0.0_c_double, 0.0_c_double, 0.0_c_double, output, parent, ledger, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         particles_are_zero(output) .and. parent_is_zero(parent) .and. &
         ledger_is_zero(ledger), &
         'event rejects a nonunit direction and clears all outputs')

    call poison_particles(output)
    call poison_parent(parent)
    call poison_ledger(ledger)
    call fusion_reaction_lab_event(FUSION_DT_ALPHAN, &
         FUSION_REACTANT_CLASSICAL_BUDGET, pa, pb, direction, &
         1.0_c_double, 0.0_c_double, 0.0_c_double, output, parent, ledger, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         particles_are_zero(output) .and. parent_is_zero(parent) .and. &
         ledger_is_zero(ledger), &
         'two-product event rejects a nonzero q and clears all outputs')
  end subroutine verify_errors_and_clearing

end program test_fusion_reaction_event_fortran
