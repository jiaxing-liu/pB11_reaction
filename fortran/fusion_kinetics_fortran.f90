module fusion_kinetics_fortran
  !! ISO_C_BINDING wrappers for the energy-space trial step and Coulomb
  !! coefficients.  The wrappers validate Fortran array extents before
  !! taking C addresses and clear all outputs on every rejected call.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_ptr, &
       c_loc
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  ! These fields are in the exact order of fusion_kinetic_ledger_v1.
  type, bind(C), public :: fusion_kinetic_ledger_v1
     real(c_double) :: initial_number_m3
     real(c_double) :: final_number_m3
     real(c_double) :: initial_energy_J_m3
     real(c_double) :: final_energy_J_m3
     real(c_double) :: born_number_m3
     real(c_double) :: born_energy_J_m3
     real(c_double) :: escaped_number_m3
     real(c_double) :: escaped_energy_J_m3
     real(c_double) :: thermalized_number_m3
     real(c_double) :: thermalized_energy_J_m3
     real(c_double) :: particle_balance_error_m3
     real(c_double) :: energy_balance_error_J_m3
  end type fusion_kinetic_ledger_v1

  ! These fields are in the exact order of fusion_maxwellian_bath_v1.
  type, bind(C), public :: fusion_maxwellian_bath_v1
     real(c_double) :: density_m3
     real(c_double) :: mass_kg
     real(c_double) :: mean_charge_squared
     real(c_double) :: kT_J
     real(c_double) :: coulomb_log
  end type fusion_maxwellian_bath_v1

  ! These fields are in the exact order of fusion_coulomb_energy_v1.
  type, bind(C), public :: fusion_coulomb_energy_v1
     real(c_double) :: diffusion_J2_s
     real(c_double) :: mean_energy_rate_J_s
  end type fusion_coulomb_energy_v1

  public :: fusion_energy_fp_trial
  public :: fusion_coulomb_energy

  interface
     function c_fusion_energy_fp_trial(cells, baths, dt_s, edges, old_number, &
          bath_kT, diffusion, birth, escape, thermalization, trial, heat, &
          ledger) bind(C, name="fusion_c_energy_fp_trial") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: cells
       integer(c_int), value :: baths
       real(c_double), value :: dt_s
       type(c_ptr), value :: edges
       type(c_ptr), value :: old_number
       type(c_ptr), value :: bath_kT
       type(c_ptr), value :: diffusion
       type(c_ptr), value :: birth
       type(c_ptr), value :: escape
       real(c_double), value :: thermalization
       type(c_ptr), value :: trial
       type(c_ptr), value :: heat
       type(c_ptr), value :: ledger
       integer(c_int) :: status
     end function c_fusion_energy_fp_trial

     function c_fusion_coulomb_energy(energy_J, mass_kg, charge_number, &
          bath, out) bind(C, name="fusion_c_coulomb_energy") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: energy_J
       real(c_double), value :: mass_kg
       real(c_double), value :: charge_number
       type(c_ptr), value :: bath
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_coulomb_energy
  end interface

contains

  subroutine clear_ledger(ledger)
    type(fusion_kinetic_ledger_v1), intent(out) :: ledger

    ledger%initial_number_m3 = 0.0_c_double
    ledger%final_number_m3 = 0.0_c_double
    ledger%initial_energy_J_m3 = 0.0_c_double
    ledger%final_energy_J_m3 = 0.0_c_double
    ledger%born_number_m3 = 0.0_c_double
    ledger%born_energy_J_m3 = 0.0_c_double
    ledger%escaped_number_m3 = 0.0_c_double
    ledger%escaped_energy_J_m3 = 0.0_c_double
    ledger%thermalized_number_m3 = 0.0_c_double
    ledger%thermalized_energy_J_m3 = 0.0_c_double
    ledger%particle_balance_error_m3 = 0.0_c_double
    ledger%energy_balance_error_J_m3 = 0.0_c_double
  end subroutine clear_ledger

  subroutine clear_coulomb(out)
    type(fusion_coulomb_energy_v1), intent(out) :: out

    out%diffusion_J2_s = 0.0_c_double
    out%mean_energy_rate_J_s = 0.0_c_double
  end subroutine clear_coulomb

  subroutine fusion_energy_fp_trial(cells, baths, dt_s, edges_J, old_number_m3, &
       bath_kT_J, diffusion_J2_s, birth_m3_s, escape_s_inv, &
       thermalization_s_inv, trial_number_m3, heat_to_bath_J_m3, ledger, &
       status)
    integer(c_int), intent(in) :: cells, baths
    real(c_double), intent(in) :: dt_s
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: old_number_m3(:)
    real(c_double), intent(in), target, contiguous :: bath_kT_J(:)
    real(c_double), intent(in), target, contiguous :: diffusion_J2_s(:)
    real(c_double), intent(in), target, contiguous :: birth_m3_s(:)
    real(c_double), intent(in), target, contiguous :: escape_s_inv(:)
    real(c_double), intent(in) :: thermalization_s_inv
    real(c_double), intent(out), target, contiguous :: trial_number_m3(:)
    real(c_double), intent(out), target, contiguous :: heat_to_bath_J_m3(:)
    type(fusion_kinetic_ledger_v1), intent(out), target :: ledger
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges, expected_diff
    type(c_ptr) :: edges_ptr, old_ptr, bath_ptr, diffusion_ptr
    type(c_ptr) :: birth_ptr, escape_ptr, trial_ptr, heat_ptr, ledger_ptr

    call clear_ledger(ledger)
    trial_number_m3 = 0.0_c_double
    heat_to_bath_J_m3 = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! The C API uses int dimensions.  Check additions and products before
    ! converting them to array extents, so malformed dimensions cannot wrap.
    if (cells < 1_c_int .or. baths < 0_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(old_number_m3, kind=c_int) /= cells) return
    if (size(birth_m3_s, kind=c_int) /= cells) return
    if (size(escape_s_inv, kind=c_int) /= cells) return
    if (size(trial_number_m3, kind=c_int) /= cells) return
    if (size(bath_kT_J, kind=c_int) /= baths) return
    if (size(heat_to_bath_J_m3, kind=c_int) /= baths) return

    expected_diff = 0_c_int
    if (cells > 1_c_int) then
       if (baths > huge(expected_diff) / (cells - 1_c_int)) return
       expected_diff = baths * (cells - 1_c_int)
    end if
    if (size(diffusion_J2_s, kind=c_int) /= expected_diff) return

    ! A zero-extent Fortran array is valid for no baths and for the
    ! one-cell diffusion vector.  In those cases pass C_NULL_PTR instead of
    ! taking C_LOC of element one.  diffusion is flattened bath-major:
    ! diffusion((b-1)*(cells-1)+f) is C diffusion[(b-1)*(cells-1)+(f-1)],
    ! where Fortran b=1..baths and f=1..cells-1.  C bath and face IDs start
    ! at zero; the flat storage is therefore bath-major in both interfaces.
    edges_ptr = c_loc(edges_J(1))
    old_ptr = c_loc(old_number_m3(1))
    birth_ptr = c_loc(birth_m3_s(1))
    escape_ptr = c_loc(escape_s_inv(1))
    trial_ptr = c_loc(trial_number_m3(1))
    ledger_ptr = c_loc(ledger)
    bath_ptr = c_null_ptr
    diffusion_ptr = c_null_ptr
    heat_ptr = c_null_ptr
    if (baths > 0_c_int) then
       bath_ptr = c_loc(bath_kT_J(1))
       heat_ptr = c_loc(heat_to_bath_J_m3(1))
    end if
    if (expected_diff > 0_c_int) diffusion_ptr = c_loc(diffusion_J2_s(1))

    status = c_fusion_energy_fp_trial(cells, baths, dt_s, edges_ptr, old_ptr, &
         bath_ptr, diffusion_ptr, birth_ptr, escape_ptr, &
         thermalization_s_inv, trial_ptr, heat_ptr, ledger_ptr)
    if (status /= PB11_STATUS_OK) then
       trial_number_m3 = 0.0_c_double
       heat_to_bath_J_m3 = 0.0_c_double
       call clear_ledger(ledger)
    end if
  end subroutine fusion_energy_fp_trial

  subroutine fusion_coulomb_energy(energy_J, mass_kg, charge_number, bath, &
       out, status)
    real(c_double), intent(in) :: energy_J, mass_kg, charge_number
    type(fusion_maxwellian_bath_v1), intent(in), target :: bath
    type(fusion_coulomb_energy_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_coulomb(out)
    status = c_fusion_coulomb_energy(energy_J, mass_kg, charge_number, &
         c_loc(bath), c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_coulomb(out)
  end subroutine fusion_coulomb_energy

end module fusion_kinetics_fortran
