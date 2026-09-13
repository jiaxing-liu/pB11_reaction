module fusion_spectrum_fortran
  !! ISO_C_BINDING wrappers for the nuclear Coulomb coefficient and
  !! sequential three-alpha amplitude C APIs.
  use, intrinsic :: iso_c_binding, only : c_double, c_int
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

  public :: fusion_nuclear_coulomb
  public :: fusion_alpha_amplitudes

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

end module fusion_spectrum_fortran
