program test_fusion_nuclear_data_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use fusion_laboratory_fortran, only : fusion_particle_four_vector_v1
  use fusion_nuclear_data_fortran
  implicit none
  real(c_double), parameter :: rel_tol = 5.0e-12_c_double
  real(c_double), parameter :: abs_tol = 1.0e-30_c_double
  integer :: failures = 0

  call data_checks()
  call two_body_checks()
  call extent_checks()
  if (failures /= 0) error stop 1
  write(*, '(A)') 'All Fortran nuclear-data tests passed'

contains
  subroutine check(ok, text)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: text
    if (.not. ok) then
       write(*, '(A)') 'FAIL: ' // text
       failures = failures + 1
    end if
  end subroutine check

  logical function near(a, b)
    real(c_double), intent(in) :: a, b
    near = abs(a - b) <= abs_tol + rel_tol * max(abs(a), abs(b), 1.0e-300_c_double)
  end function near

  logical function zero_particles(p)
    type(fusion_particle_four_vector_v1), intent(in) :: p(:)
    integer :: i
    zero_particles = .true.
    do i = 1, size(p)
       zero_particles = zero_particles .and. &
            abs(p(i)%mass_kg) <= abs_tol .and. &
            abs(p(i)%kinetic_energy_J) <= abs_tol .and. &
            maxval(abs(p(i)%momentum_kg_m_s)) <= abs_tol
    end do
  end function zero_particles

  subroutine data_checks()
    type(fusion_nuclear_mass_v1) :: p, d, t, a, b, n
    type(fusion_nuclear_channel_v1) :: pb, dt
    real(c_double) :: q_pb, q_dt
    integer(c_int) :: status
    call check(c_sizeof(p) == 4_c_size_t * c_sizeof(q_pb) + &
         2_c_size_t * c_sizeof(status), 'mass layout is 4 doubles plus 2 ints')
    call check(c_sizeof(pb) == 2_c_size_t * c_sizeof(q_pb) + &
         6_c_size_t * c_sizeof(status), 'channel layout is 2 doubles plus 6 ints')
    call fusion_nuclear_mass(FUSION_PROTON, p, status); call check(status == 0, 'proton mass query')
    call fusion_nuclear_mass(FUSION_BORON11, b, status); call check(status == 0, 'boron mass query')
    call fusion_nuclear_mass(FUSION_MASS_NEUTRON, n, status); call check(status == 0, 'neutron mass query')
    call check(p%nuclear_charge == 1 .and. p%mass_number == 1 .and. &
         b%nuclear_charge == 5 .and. b%mass_number == 11 .and. &
         n%nuclear_charge == 0 .and. n%mass_number == 1, 'p/B/neutron IDs')
    call fusion_nuclear_mass(FUSION_DEUTERON, d, status)
    call fusion_nuclear_mass(FUSION_TRITON, t, status)
    call fusion_nuclear_mass(FUSION_HELIUM4, a, status)
    call fusion_nuclear_channel(FUSION_PB11_3ALPHA, pb, status)
    call check(status == 0 .and. pb%product_count == 3 .and. &
         all(pb%reactant_ids == [FUSION_PROTON,FUSION_BORON11]) .and. &
         all(pb%product_ids == FUSION_HELIUM4), 'pB channel IDs and 3 alpha products')
    call fusion_nuclear_channel(FUSION_DT_ALPHAN, dt, status)
    call check(status == 0 .and. dt%product_count == 2 .and. &
         all(dt%reactant_ids == [FUSION_DEUTERON,FUSION_TRITON]) .and. &
         dt%product_ids(1) == FUSION_HELIUM4 .and. &
         dt%product_ids(2) == FUSION_MASS_NEUTRON .and. dt%product_ids(3) == -1, &
         'DT channel IDs and two product IDs')
    q_pb = p%rest_energy_J + b%rest_energy_J - 3.0_c_double*a%rest_energy_J
    q_dt = d%rest_energy_J + t%rest_energy_J - a%rest_energy_J - n%rest_energy_J
    call check(near(pb%q_J, q_pb) .and. near(dt%q_J, q_dt), 'Q uses returned masses')
  end subroutine data_checks

  subroutine two_body_checks()
    type(fusion_nuclear_channel_v1) :: dt
    type(fusion_particle_four_vector_v1) :: output(2)
    real(c_double) :: direction(3), ecm, total_kinetic
    integer(c_int) :: status
    direction = [0.0_c_double, 0.0_c_double, 1.0_c_double]
    ecm = 1.0e-14_c_double
    call fusion_nuclear_channel(FUSION_DT_ALPHAN, dt, status)
    call fusion_nuclear_two_body_cm(FUSION_DT_ALPHAN, ecm, direction, output, status)
    total_kinetic = output(1)%kinetic_energy_J + output(2)%kinetic_energy_J
    call check(status == 0 .and. near(total_kinetic, ecm + dt%q_J), &
         'DT two-body kinetic energy equals Ecm plus Q')
    output%mass_kg = -1.0_c_double
    call fusion_nuclear_two_body_cm(FUSION_PB11_3ALPHA, ecm, direction, output, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. zero_particles(output), &
         'pB three-body channel is rejected by two-body API')
  end subroutine two_body_checks

  subroutine extent_checks()
    type(fusion_particle_four_vector_v1) :: output(2), short_output(1)
    real(c_double) :: direction(3), short_direction(2)
    integer(c_int) :: status
    direction = 0.0_c_double; direction(3) = 1.0_c_double
    short_direction = 1.0_c_double; output%mass_kg = -2.0_c_double
    call fusion_nuclear_two_body_cm(FUSION_DT_ALPHAN, 1.0_c_double, &
         short_direction, output, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. zero_particles(output), &
         'bad direction extent clears output')
    short_output%mass_kg = -3.0_c_double
    call fusion_nuclear_two_body_cm(FUSION_DT_ALPHAN, 1.0_c_double, &
         direction, short_output, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. zero_particles(short_output), &
         'bad output extent clears output')
  end subroutine extent_checks
end program test_fusion_nuclear_data_fortran
