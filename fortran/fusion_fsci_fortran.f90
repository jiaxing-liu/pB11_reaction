module fusion_fsci_fortran
  !! ISO_C_BINDING wrappers for the FSCI Coulomb/amplitude/spectrum APIs.
  !!
  !! The public result types are owned by fusion_spectrum_fortran.  This
  !! module only adds the policy-aware entry points and their policy enum.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  use fusion_spectrum_fortran, only : fusion_nuclear_coulomb_v1, &
       fusion_alpha_amplitudes_v1, fusion_alpha_spectrum_v1
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_int), parameter, public :: FUSION_ALPHA_FSCI_NONE = 0_c_int
  integer(c_int), parameter, public :: FUSION_ALPHA_FSCI_REFS2018_R16 = 1_c_int

  public :: fusion_nuclear_coulomb_radius16
  public :: fusion_alpha_amplitudes_fsci_cutoff
  public :: fusion_alpha_spectrum_model_grid

  interface
     function c_fusion_nuclear_coulomb_radius16(channel, relative_energy_J, &
          out) bind(C, name="fusion_c_nuclear_coulomb_radius16") result(status)
       import :: c_double, c_int, fusion_nuclear_coulomb_v1
       integer(c_int), value :: channel
       real(c_double), value :: relative_energy_J
       type(fusion_nuclear_coulomb_v1), intent(out) :: out
       integer(c_int) :: status
     end function c_fusion_nuclear_coulomb_radius16

     function c_fusion_alpha_amplitudes_fsci_cutoff(primary_l, policy, &
          available_energy_J, intermediate_energy_J, cos_theta, cutoff_J, &
          out, pruned) bind(C, &
          name="fusion_c_alpha_amplitudes_fsci_cutoff") result(status)
       import :: c_double, c_int, fusion_alpha_amplitudes_v1
       integer(c_int), value :: primary_l
       integer(c_int), value :: policy
       real(c_double), value :: available_energy_J
       real(c_double), value :: intermediate_energy_J
       real(c_double), value :: cos_theta
       real(c_double), value :: cutoff_J
       type(fusion_alpha_amplitudes_v1), intent(out) :: out
       integer(c_int), intent(out) :: pruned
       integer(c_int) :: status
     end function c_fusion_alpha_amplitudes_fsci_cutoff

     function c_fusion_alpha_spectrum_model_grid(mode, policy, &
          available_energy_J, cutoff_J, l1_fraction, relative_phase, nq, nc, &
          cells, edges, birth, out) bind(C, &
          name="fusion_c_alpha_spectrum_model_grid") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: mode
       integer(c_int), value :: policy
       real(c_double), value :: available_energy_J
       real(c_double), value :: cutoff_J
       real(c_double), value :: l1_fraction
       real(c_double), value :: relative_phase
       integer(c_int), value :: nq
       integer(c_int), value :: nc
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_alpha_spectrum_model_grid
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

  subroutine fusion_nuclear_coulomb_radius16(channel, relative_energy_J, &
       out, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: relative_energy_J
    type(fusion_nuclear_coulomb_v1), intent(out) :: out
    integer(c_int), intent(out) :: status

    call clear_nuclear_coulomb(out)
    status = c_fusion_nuclear_coulomb_radius16(channel, relative_energy_J, out)
    if (status /= PB11_STATUS_OK) call clear_nuclear_coulomb(out)
  end subroutine fusion_nuclear_coulomb_radius16

  subroutine fusion_alpha_amplitudes_fsci_cutoff(primary_l, policy, &
       available_energy_J, intermediate_energy_J, cos_theta, cutoff_J, out, &
       pruned, status)
    integer(c_int), intent(in) :: primary_l, policy
    real(c_double), intent(in) :: available_energy_J, intermediate_energy_J
    real(c_double), intent(in) :: cos_theta, cutoff_J
    type(fusion_alpha_amplitudes_v1), intent(out) :: out
    integer(c_int), intent(out) :: pruned, status

    call clear_alpha_amplitudes(out)
    pruned = 0_c_int
    status = c_fusion_alpha_amplitudes_fsci_cutoff(primary_l, policy, &
         available_energy_J, intermediate_energy_J, cos_theta, cutoff_J, out, &
         pruned)
    if (status /= PB11_STATUS_OK) then
       call clear_alpha_amplitudes(out)
       pruned = 0_c_int
    end if
  end subroutine fusion_alpha_amplitudes_fsci_cutoff

  subroutine fusion_alpha_spectrum_model_grid(mode, policy, available_energy_J, &
       cutoff_J, l1_fraction, relative_phase, nq, nc, cells, edges_J, &
       birth_per_event, out, status)
    integer(c_int), intent(in) :: mode, policy, nq, nc, cells
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

    ! Validate extents before C_LOC.  The C API has no zero-cell pointer form.
    if (cells < 1_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(birth_per_event, kind=c_int) /= cells) return

    edges_ptr = c_loc(edges_J(1))
    birth_ptr = c_loc(birth_per_event(1))
    out_ptr = c_loc(out)
    status = c_fusion_alpha_spectrum_model_grid(mode, policy, &
         available_energy_J, cutoff_J, l1_fraction, relative_phase, nq, nc, &
         cells, edges_ptr, birth_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       birth_per_event = 0.0_c_double
       call clear_alpha_spectrum(out)
    end if
  end subroutine fusion_alpha_spectrum_model_grid

end module fusion_fsci_fortran
