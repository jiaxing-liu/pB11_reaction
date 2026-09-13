module fusion_nuclear_data_fortran
  !! ISO_C_BINDING wrappers for the fixed nuclear mass/channel data and the
  !! deterministic two-body CM product constructor.
  !!
  !! Nuclear choices and validation remain in the C implementation.  This
  !! layer mirrors the frozen layouts, validates Fortran extents before
  !! C_LOC, and clears all nonnull outputs on every rejected call.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  use fusion_laboratory_fortran, only : fusion_particle_four_vector_v1
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: FUSION_PROTON = 0_c_int
  integer(c_int), parameter, public :: FUSION_DEUTERON = 1_c_int
  integer(c_int), parameter, public :: FUSION_TRITON = 2_c_int
  integer(c_int), parameter, public :: FUSION_HELIUM3 = 3_c_int
  integer(c_int), parameter, public :: FUSION_HELIUM4 = 4_c_int
  integer(c_int), parameter, public :: FUSION_BORON11 = 5_c_int
  integer(c_int), parameter, public :: FUSION_SPECIES_COUNT = 6_c_int
  integer(c_int), parameter, public :: FUSION_MASS_NEUTRON = 6_c_int
  integer(c_int), parameter, public :: FUSION_MASS_ELECTRON = 7_c_int

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = 0_c_int
  integer(c_int), parameter, public :: FUSION_DD_TP = 1_c_int
  integer(c_int), parameter, public :: FUSION_DD_HE3N = 2_c_int
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = 3_c_int
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = 4_c_int
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = 5_c_int

  type, bind(C), public :: fusion_nuclear_mass_v1
     real(c_double) :: mass_kg
     real(c_double) :: rest_energy_J
     real(c_double) :: reference_mass_u
     real(c_double) :: known_uncertainty_scale_kg
     integer(c_int) :: nuclear_charge
     integer(c_int) :: mass_number
  end type fusion_nuclear_mass_v1

  type, bind(C), public :: fusion_nuclear_channel_v1
     real(c_double) :: q_J
     real(c_double) :: known_uncertainty_scale_J
     integer(c_int) :: reactant_ids(2)
     integer(c_int) :: product_ids(3)
     integer(c_int) :: product_count
  end type fusion_nuclear_channel_v1

  public :: fusion_nuclear_mass
  public :: fusion_nuclear_channel
  public :: fusion_nuclear_two_body_cm

  interface
     function c_fusion_nuclear_mass(particle_id, out) bind(C, &
          name="fusion_c_nuclear_mass") result(status)
       import :: c_int, c_ptr
       integer(c_int), value :: particle_id
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_nuclear_mass

     function c_fusion_nuclear_channel(channel, out) bind(C, &
          name="fusion_c_nuclear_channel") result(status)
       import :: c_int, c_ptr
       integer(c_int), value :: channel
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_nuclear_channel

     function c_fusion_nuclear_two_body_cm(channel, reactant_energy, &
          direction, output) bind(C, &
          name="fusion_c_nuclear_two_body_cm") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       real(c_double), value :: reactant_energy
       type(c_ptr), value :: direction
       type(c_ptr), value :: output
       integer(c_int) :: status
     end function c_fusion_nuclear_two_body_cm
  end interface

contains

  subroutine clear_mass(out)
    type(fusion_nuclear_mass_v1), intent(out) :: out

    out%mass_kg = 0.0_c_double
    out%rest_energy_J = 0.0_c_double
    out%reference_mass_u = 0.0_c_double
    out%known_uncertainty_scale_kg = 0.0_c_double
    out%nuclear_charge = 0_c_int
    out%mass_number = 0_c_int
  end subroutine clear_mass

  subroutine clear_channel(out)
    type(fusion_nuclear_channel_v1), intent(out) :: out

    out%q_J = 0.0_c_double
    out%known_uncertainty_scale_J = 0.0_c_double
    out%reactant_ids = 0_c_int
    out%product_ids = 0_c_int
    out%product_count = 0_c_int
  end subroutine clear_channel

  subroutine clear_particles(output)
    type(fusion_particle_four_vector_v1), intent(out) :: output(:)
    integer :: i

    do i = 1, size(output)
       output(i)%mass_kg = 0.0_c_double
       output(i)%kinetic_energy_J = 0.0_c_double
       output(i)%momentum_kg_m_s = 0.0_c_double
    end do
  end subroutine clear_particles

  subroutine fusion_nuclear_mass(particle_id, out, status)
    integer(c_int), intent(in) :: particle_id
    type(fusion_nuclear_mass_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_mass(out)
    status = c_fusion_nuclear_mass(particle_id, c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_mass(out)
  end subroutine fusion_nuclear_mass

  subroutine fusion_nuclear_channel(channel, out, status)
    integer(c_int), intent(in) :: channel
    type(fusion_nuclear_channel_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_channel(out)
    status = c_fusion_nuclear_channel(channel, c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_channel(out)
  end subroutine fusion_nuclear_channel

  subroutine fusion_nuclear_two_body_cm(channel, reactant_CM_kinetic_J, &
       direction, output, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: reactant_CM_kinetic_J
    real(c_double), intent(in), target, contiguous :: direction(:)
    type(fusion_particle_four_vector_v1), intent(out), target, contiguous :: &
         output(:)
    integer(c_int), intent(out) :: status

    call clear_particles(output)
    status = PB11_STATUS_INVALID_ARGUMENT

    ! The C ABI has direction[3] and output[2].  Do not form C_LOC until both
    ! extents are exact, including for a malformed zero-length actual.
    if (size(direction, kind=c_int) /= 3_c_int) return
    if (size(output, kind=c_int) /= 2_c_int) return

    status = c_fusion_nuclear_two_body_cm(channel, reactant_CM_kinetic_J, &
         c_loc(direction(1)), c_loc(output(1)))
    if (status /= PB11_STATUS_OK) call clear_particles(output)
  end subroutine fusion_nuclear_two_body_cm

end module fusion_nuclear_data_fortran
