module fusion_laboratory_fortran
  !! ISO_C_BINDING wrappers for the laboratory two-body CM and Lorentz
  !! particle boost APIs.  All arrays passed to C are checked for their
  !! required extents and are CONTIGUOUS so C_LOC is safe for strided
  !! Fortran actual arguments through compiler temporaries.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  ! Exact C layout: five c_double fields, with the final three doubles
  ! represented by the interoperable length-three array.
  type, bind(C), public :: fusion_particle_four_vector_v1
     real(c_double) :: mass_kg
     real(c_double) :: kinetic_energy_J
     real(c_double) :: momentum_kg_m_s(3)
  end type fusion_particle_four_vector_v1

  ! Exact C layout: four c_double fields.
  type, bind(C), public :: fusion_boost_ledger_v1
     real(c_double) :: expected_kinetic_energy_J
     real(c_double) :: output_kinetic_energy_J
     real(c_double) :: energy_residual_J
     real(c_double) :: momentum_residual_kg_m_s
  end type fusion_boost_ledger_v1

  public :: fusion_two_body_cm
  public :: fusion_boost_particles

  interface
     function c_fusion_two_body_cm(mass_a_kg, mass_b_kg, available_energy_J, &
          direction, output) bind(C, name="fusion_c_two_body_cm") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: mass_a_kg
       real(c_double), value :: mass_b_kg
       real(c_double), value :: available_energy_J
       type(c_ptr), value :: direction
       type(c_ptr), value :: output
       integer(c_int) :: status
     end function c_fusion_two_body_cm

     function c_fusion_boost_particles(count, velocity_m_s, input, output, &
          ledger) bind(C, name="fusion_c_boost_particles") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: count
       type(c_ptr), value :: velocity_m_s
       type(c_ptr), value :: input
       type(c_ptr), value :: output
       type(c_ptr), value :: ledger
       integer(c_int) :: status
     end function c_fusion_boost_particles
  end interface

contains

  subroutine clear_particle(particle)
    type(fusion_particle_four_vector_v1), intent(out) :: particle

    particle%mass_kg = 0.0_c_double
    particle%kinetic_energy_J = 0.0_c_double
    particle%momentum_kg_m_s = 0.0_c_double
  end subroutine clear_particle

  subroutine clear_particles(particles)
    type(fusion_particle_four_vector_v1), intent(out) :: particles(:)
    integer :: i

    do i = 1, size(particles)
       call clear_particle(particles(i))
    end do
  end subroutine clear_particles

  subroutine clear_ledger(ledger)
    type(fusion_boost_ledger_v1), intent(out) :: ledger

    ledger%expected_kinetic_energy_J = 0.0_c_double
    ledger%output_kinetic_energy_J = 0.0_c_double
    ledger%energy_residual_J = 0.0_c_double
    ledger%momentum_residual_kg_m_s = 0.0_c_double
  end subroutine clear_ledger

  subroutine fusion_two_body_cm(mass_a_kg, mass_b_kg, available_energy_J, &
       direction, output, status)
    real(c_double), intent(in) :: mass_a_kg, mass_b_kg, available_energy_J
    real(c_double), intent(in), target, contiguous :: direction(:)
    type(fusion_particle_four_vector_v1), intent(out), target, contiguous :: &
         output(:)
    integer(c_int), intent(out) :: status

    call clear_particles(output)
    status = PB11_STATUS_INVALID_ARGUMENT

    if (size(direction, kind=c_int) /= 3_c_int) return
    if (size(output, kind=c_int) /= 2_c_int) return

    status = c_fusion_two_body_cm(mass_a_kg, mass_b_kg, available_energy_J, &
         c_loc(direction(1)), c_loc(output(1)))
    if (status /= PB11_STATUS_OK) call clear_particles(output)
  end subroutine fusion_two_body_cm

  subroutine fusion_boost_particles(count, velocity_m_s, input, output, &
       ledger, status)
    integer(c_int), intent(in) :: count
    real(c_double), intent(in), target, contiguous :: velocity_m_s(:)
    type(fusion_particle_four_vector_v1), intent(in), target, contiguous :: &
         input(:)
    type(fusion_particle_four_vector_v1), intent(out), target, contiguous :: &
         output(:)
    type(fusion_boost_ledger_v1), intent(out), target :: ledger
    integer(c_int), intent(out) :: status

    call clear_particles(output)
    call clear_ledger(ledger)
    status = PB11_STATUS_INVALID_ARGUMENT

    if (count < 1_c_int) return
    if (size(velocity_m_s, kind=c_int) /= 3_c_int) return
    if (size(input, kind=c_int) /= count) return
    if (size(output, kind=c_int) /= count) return

    status = c_fusion_boost_particles(count, c_loc(velocity_m_s(1)), &
         c_loc(input(1)), c_loc(output(1)), c_loc(ledger))
    if (status /= PB11_STATUS_OK) then
       call clear_particles(output)
       call clear_ledger(ledger)
    end if
  end subroutine fusion_boost_particles

end module fusion_laboratory_fortran
