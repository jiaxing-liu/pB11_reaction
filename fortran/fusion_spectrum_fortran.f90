module fusion_spectrum_fortran
  !! ISO_C_BINDING wrappers for the nuclear Coulomb coefficient and
  !! sequential three-alpha amplitude C APIs.
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

  integer(c_int), parameter, public :: FUSION_ALPHA_BE8_L1 = 0_c_int
  integer(c_int), parameter, public :: FUSION_ALPHA_BE8_L2 = 1_c_int
  integer(c_int), parameter, public :: FUSION_ALPHA_BE8_L3 = 2_c_int
  integer(c_int), parameter, public :: FUSION_ALPHA_ALPHA_L0 = 3_c_int
  integer(c_int), parameter, public :: FUSION_ALPHA_ALPHA_L2 = 4_c_int
  integer(c_int), parameter, public :: FUSION_NUCLEAR_CHANNEL_COUNT = 5_c_int

  ! Exact C layout: five consecutive c_double fields.
  type, bind(C), public :: fusion_nuclear_coulomb_v1
     real(c_double) :: log_penetrability
     real(c_double) :: shift
     real(c_double) :: phase_real
     real(c_double) :: phase_imag
     real(c_double) :: rho
  end type fusion_nuclear_coulomb_v1

  ! Exact C layout: four length-five arrays followed by one c_double.
  type, bind(C), public :: fusion_alpha_amplitudes_v1
     real(c_double) :: unsym_real(5)
     real(c_double) :: unsym_imag(5)
     real(c_double) :: sym_real(5)
     real(c_double) :: sym_imag(5)
     real(c_double) :: phase_space_J
  end type fusion_alpha_amplitudes_v1

  ! Exact C layout: eleven consecutive c_double fields followed by two
  ! c_int fields.
  type, bind(C), public :: fusion_alpha_spectrum_v1
     real(c_double) :: mapped_number
     real(c_double) :: mapped_energy_J
     real(c_double) :: below_number
     real(c_double) :: below_energy_J
     real(c_double) :: above_number
     real(c_double) :: above_energy_J
     real(c_double) :: number_residual
     real(c_double) :: energy_residual_J
     real(c_double) :: normalization_J2
     real(c_double) :: l1_normalization_J2
     real(c_double) :: l3_normalization_J2
     integer(c_int) :: quadrature_events
     integer(c_int) :: pruned_events
  end type fusion_alpha_spectrum_v1

  public :: fusion_nuclear_coulomb
  public :: fusion_alpha_amplitudes
  public :: fusion_alpha_amplitudes_cutoff
  public :: fusion_alpha_spectrum_grid

  interface
     function c_fusion_nuclear_coulomb(channel, relative_energy_J, out) &
          bind(C, name="fusion_c_nuclear_coulomb") result(status)
       import :: c_double, c_int, fusion_nuclear_coulomb_v1
       integer(c_int), value :: channel
       real(c_double), value :: relative_energy_J
       type(fusion_nuclear_coulomb_v1), intent(out) :: out
       integer(c_int) :: status
     end function c_fusion_nuclear_coulomb

     function c_fusion_alpha_amplitudes(primary_l, available_energy_J, &
          intermediate_energy_J, cos_theta, out) &
          bind(C, name="fusion_c_alpha_amplitudes") result(status)
       import :: c_double, c_int, fusion_alpha_amplitudes_v1
       integer(c_int), value :: primary_l
       real(c_double), value :: available_energy_J
       real(c_double), value :: intermediate_energy_J
       real(c_double), value :: cos_theta
       type(fusion_alpha_amplitudes_v1), intent(out) :: out
       integer(c_int) :: status
     end function c_fusion_alpha_amplitudes

     function c_fusion_alpha_amplitudes_cutoff(primary_l, available_energy_J, &
          intermediate_energy_J, cos_theta, cutoff_J, out, pruned) &
          bind(C, name="fusion_c_alpha_amplitudes_cutoff") result(status)
       import :: c_double, c_int, fusion_alpha_amplitudes_v1
       integer(c_int), value :: primary_l
       real(c_double), value :: available_energy_J
       real(c_double), value :: intermediate_energy_J
       real(c_double), value :: cos_theta
       real(c_double), value :: cutoff_J
       type(fusion_alpha_amplitudes_v1), intent(out) :: out
       integer(c_int), intent(out) :: pruned
       integer(c_int) :: status
     end function c_fusion_alpha_amplitudes_cutoff

     function c_fusion_alpha_spectrum_grid(mode, available_energy_J, &
          cutoff_J, l1_fraction, relative_phase, nq, ncos, cells, edges, &
          birth, out) bind(C, name="fusion_c_alpha_spectrum_grid") &
          result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: mode
       real(c_double), value :: available_energy_J
       real(c_double), value :: cutoff_J
       real(c_double), value :: l1_fraction
       real(c_double), value :: relative_phase
       integer(c_int), value :: nq
       integer(c_int), value :: ncos
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_alpha_spectrum_grid
  end interface

contains

  subroutine clear_nuclear_coulomb(out)
    type(fusion_nuclear_coulomb_v1), intent(out) :: out

    out%log_penetrability = 0.0_c_double
    out%shift = 0.0_c_double
    out%phase_real = 0.0_c_double
    out%phase_imag = 0.0_c_double
    out%rho = 0.0_c_double
  end subroutine clear_nuclear_coulomb

  subroutine clear_alpha_amplitudes(out)
    type(fusion_alpha_amplitudes_v1), intent(out) :: out

    out%unsym_real = 0.0_c_double
    out%unsym_imag = 0.0_c_double
    out%sym_real = 0.0_c_double
    out%sym_imag = 0.0_c_double
    out%phase_space_J = 0.0_c_double
  end subroutine clear_alpha_amplitudes

  subroutine clear_alpha_spectrum(out)
    type(fusion_alpha_spectrum_v1), intent(out) :: out

    out%mapped_number = 0.0_c_double
    out%mapped_energy_J = 0.0_c_double
    out%below_number = 0.0_c_double
    out%below_energy_J = 0.0_c_double
    out%above_number = 0.0_c_double
    out%above_energy_J = 0.0_c_double
    out%number_residual = 0.0_c_double
    out%energy_residual_J = 0.0_c_double
    out%normalization_J2 = 0.0_c_double
    out%l1_normalization_J2 = 0.0_c_double
    out%l3_normalization_J2 = 0.0_c_double
    out%quadrature_events = 0_c_int
    out%pruned_events = 0_c_int
  end subroutine clear_alpha_spectrum

  subroutine fusion_nuclear_coulomb(channel, relative_energy_J, out, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: relative_energy_J
    type(fusion_nuclear_coulomb_v1), intent(out) :: out
    integer(c_int), intent(out) :: status

    call clear_nuclear_coulomb(out)
    status = c_fusion_nuclear_coulomb(channel, relative_energy_J, out)
    if (status /= PB11_STATUS_OK) call clear_nuclear_coulomb(out)
  end subroutine fusion_nuclear_coulomb

  subroutine fusion_alpha_amplitudes(primary_l, available_energy_J, &
       intermediate_energy_J, cos_theta, out, status)
    integer(c_int), intent(in) :: primary_l
    real(c_double), intent(in) :: available_energy_J
    real(c_double), intent(in) :: intermediate_energy_J
    real(c_double), intent(in) :: cos_theta
    type(fusion_alpha_amplitudes_v1), intent(out) :: out
    integer(c_int), intent(out) :: status

    call clear_alpha_amplitudes(out)
    status = c_fusion_alpha_amplitudes(primary_l, available_energy_J, &
         intermediate_energy_J, cos_theta, out)
    if (status /= PB11_STATUS_OK) call clear_alpha_amplitudes(out)
  end subroutine fusion_alpha_amplitudes

  subroutine fusion_alpha_amplitudes_cutoff(primary_l, available_energy_J, &
       intermediate_energy_J, cos_theta, cutoff_J, out, pruned, status)
    integer(c_int), intent(in) :: primary_l
    real(c_double), intent(in) :: available_energy_J
    real(c_double), intent(in) :: intermediate_energy_J, cos_theta, cutoff_J
    type(fusion_alpha_amplitudes_v1), intent(out) :: out
    integer(c_int), intent(out) :: pruned, status

    call clear_alpha_amplitudes(out)
    pruned = 0_c_int
    status = PB11_STATUS_INVALID_ARGUMENT
    status = c_fusion_alpha_amplitudes_cutoff(primary_l, available_energy_J, &
         intermediate_energy_J, cos_theta, cutoff_J, out, pruned)
    if (status /= PB11_STATUS_OK) then
       call clear_alpha_amplitudes(out)
       pruned = 0_c_int
    end if
  end subroutine fusion_alpha_amplitudes_cutoff

  subroutine fusion_alpha_spectrum_grid(mode, available_energy_J, cutoff_J, &
       l1_fraction, relative_phase, nq, ncos, cells, edges_J, &
       birth_per_event, out, status)
    integer(c_int), intent(in) :: mode, nq, ncos, cells
    real(c_double), intent(in) :: available_energy_J, cutoff_J
    real(c_double), intent(in) :: l1_fraction, relative_phase
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(out), target, contiguous :: birth_per_event(:)
    type(fusion_alpha_spectrum_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges
    type(c_ptr) :: edges_ptr, birth_ptr, out_ptr

    call clear_alpha_spectrum(out)
    birth_per_event = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Validate extents before C_LOC.  The C API accepts only one or more
    ! cells, so neither array has a valid zero-size pointer here.
    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(birth_per_event, kind=c_int) /= cells) return

    edges_ptr = c_loc(edges_J(1))
    birth_ptr = c_loc(birth_per_event(1))
    out_ptr = c_loc(out)
    status = c_fusion_alpha_spectrum_grid(mode, available_energy_J, cutoff_J, &
         l1_fraction, relative_phase, nq, ncos, cells, edges_ptr, birth_ptr, &
         out_ptr)
    if (status /= PB11_STATUS_OK) then
       birth_per_event = 0.0_c_double
       call clear_alpha_spectrum(out)
    end if
  end subroutine fusion_alpha_spectrum_grid

end module fusion_spectrum_fortran
