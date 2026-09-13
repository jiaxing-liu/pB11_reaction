module fusion_reaction_event_fortran
  !! ISO_C_BINDING wrappers for the deterministic reaction-parent and
  !! laboratory-event APIs.
  !!
  !! The C layer owns the nuclear channel and kinematic policy.  This module
  !! supplies the exact C-interoperable parent record, reuses the laboratory
  !! particle and boost-ledger records, and validates every assumed-shape
  !! extent before forming a C_LOC.  A malformed call, or a C-layer error,
  !! leaves every available output cleared.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  use fusion_laboratory_fortran, only : fusion_particle_four_vector_v1, &
       fusion_boost_ledger_v1
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = 0_c_int
  integer(c_int), parameter, public :: FUSION_DD_TP = 1_c_int
  integer(c_int), parameter, public :: FUSION_DD_HE3N = 2_c_int
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = 3_c_int
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = 4_c_int
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = 5_c_int

  integer(c_int), parameter, public :: &
       FUSION_REACTANT_CLASSICAL_BUDGET = 0_c_int
  integer(c_int), parameter, public :: FUSION_REACTANT_ON_SHELL = 1_c_int

  ! Exact C layout: sixteen consecutive c_double values followed by two C
  ! ints.  The arrays retain their C order and are therefore interoperable
  ! with fusion_reaction_parent_v1 in fusion_reaction_event.h.
  type, bind(C), public :: fusion_reaction_parent_v1
     real(c_double) :: reactant_kinetic_J(2)
     real(c_double) :: classical_kinetic_J(2)
     real(c_double) :: on_shell_kinetic_J(2)
     real(c_double) :: relative_classical_energy_J
     real(c_double) :: available_cm_energy_J
     real(c_double) :: boost_velocity_m_s(3)
     real(c_double) :: expected_product_lab_kinetic_J
     real(c_double) :: classical_minus_on_shell_J
     real(c_double) :: momentum_sum_kg_m_s(3)
     integer(c_int) :: product_count
     integer(c_int) :: convention
  end type fusion_reaction_parent_v1

  public :: fusion_reaction_parent
  public :: fusion_reaction_lab_event

  interface
     function c_fusion_reaction_parent(channel, convention, pa, pb, out) &
          bind(C, name="fusion_c_reaction_parent") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       integer(c_int), value :: convention
       type(c_ptr), value :: pa
       type(c_ptr), value :: pb
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_reaction_parent

     function c_fusion_reaction_lab_event(channel, convention, pa, pb, &
          direction, q_J, cosine, azimuth, output, parent, ledger) &
          bind(C, name="fusion_c_reaction_lab_event") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       integer(c_int), value :: convention
       type(c_ptr), value :: pa
       type(c_ptr), value :: pb
       type(c_ptr), value :: direction
       real(c_double), value :: q_J
       real(c_double), value :: cosine
       real(c_double), value :: azimuth
       type(c_ptr), value :: output
       type(c_ptr), value :: parent
       type(c_ptr), value :: ledger
       integer(c_int) :: status
     end function c_fusion_reaction_lab_event
  end interface

contains

  subroutine clear_parent(out)
    type(fusion_reaction_parent_v1), intent(out) :: out

    out%reactant_kinetic_J = 0.0_c_double
    out%classical_kinetic_J = 0.0_c_double
    out%on_shell_kinetic_J = 0.0_c_double
    out%relative_classical_energy_J = 0.0_c_double
    out%available_cm_energy_J = 0.0_c_double
    out%boost_velocity_m_s = 0.0_c_double
    out%expected_product_lab_kinetic_J = 0.0_c_double
    out%classical_minus_on_shell_J = 0.0_c_double
    out%momentum_sum_kg_m_s = 0.0_c_double
    out%product_count = 0_c_int
    out%convention = 0_c_int
  end subroutine clear_parent

  subroutine clear_particles(output)
    type(fusion_particle_four_vector_v1), intent(out) :: output(:)
    integer :: i

    do i = 1, size(output)
       output(i)%mass_kg = 0.0_c_double
       output(i)%kinetic_energy_J = 0.0_c_double
       output(i)%momentum_kg_m_s = 0.0_c_double
    end do
  end subroutine clear_particles

  subroutine clear_ledger(ledger)
    type(fusion_boost_ledger_v1), intent(out) :: ledger

    ledger%expected_kinetic_energy_J = 0.0_c_double
    ledger%output_kinetic_energy_J = 0.0_c_double
    ledger%energy_residual_J = 0.0_c_double
    ledger%momentum_residual_kg_m_s = 0.0_c_double
  end subroutine clear_ledger

  subroutine fusion_reaction_parent(channel, convention, pa, pb, out, status)
    integer(c_int), intent(in) :: channel, convention
    real(c_double), intent(in), target, contiguous :: pa(:), pb(:)
    type(fusion_reaction_parent_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_parent(out)
    status = PB11_STATUS_INVALID_ARGUMENT

    ! The C ABI has exactly one three-vector for each reactant.  In
    ! particular, malformed zero-length actuals are rejected before C_LOC.
    if (size(pa, kind=c_int) /= 3_c_int) return
    if (size(pb, kind=c_int) /= 3_c_int) return

    status = c_fusion_reaction_parent(channel, convention, c_loc(pa(1)), &
         c_loc(pb(1)), c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_parent(out)
  end subroutine fusion_reaction_parent

  subroutine fusion_reaction_lab_event(channel, convention, pa, pb, &
       direction, q_J, cosine, azimuth, output, parent, ledger, status)
    integer(c_int), intent(in) :: channel, convention
    real(c_double), intent(in), target, contiguous :: pa(:), pb(:)
    real(c_double), intent(in), target, contiguous :: direction(:)
    real(c_double), intent(in) :: q_J, cosine, azimuth
    type(fusion_particle_four_vector_v1), intent(out), target, contiguous :: &
         output(:)
    type(fusion_reaction_parent_v1), intent(out), target :: parent
    type(fusion_boost_ledger_v1), intent(out), target :: ledger
    integer(c_int), intent(out) :: status

    call clear_particles(output)
    call clear_parent(parent)
    call clear_ledger(ledger)
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Every extent is checked before any C_LOC, including output(:).  This
    ! keeps zero-length and otherwise malformed Fortran actuals safe.
    if (size(pa, kind=c_int) /= 3_c_int) return
    if (size(pb, kind=c_int) /= 3_c_int) return
    if (size(direction, kind=c_int) /= 3_c_int) return
    if (size(output, kind=c_int) /= 3_c_int) return

    status = c_fusion_reaction_lab_event(channel, convention, c_loc(pa(1)), &
         c_loc(pb(1)), c_loc(direction(1)), q_J, cosine, azimuth, &
         c_loc(output(1)), c_loc(parent), c_loc(ledger))
    if (status /= PB11_STATUS_OK) then
       call clear_particles(output)
       call clear_parent(parent)
       call clear_ledger(ledger)
    end if
  end subroutine fusion_reaction_lab_event

end module fusion_reaction_event_fortran
