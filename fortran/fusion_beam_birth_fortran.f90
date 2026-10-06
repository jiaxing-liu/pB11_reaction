module fusion_beam_birth_fortran
  !! ISO_C_BINDING wrapper for the finite-temperature beam-product-birth grid API.
  !!
  !! The result embeds the existing thermal-birth C result layout.  BIRTH is
  !! declared (cells,7); its column-major storage is the C species-major
  !! layout [species][cell] required by fusion_c_beam_birth_grid.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  use fusion_thermal_birth_fortran, only : fusion_thermal_birth_v1, &
       pb11_status_ok_base => PB11_STATUS_OK, &
       pb11_status_null_output_base => PB11_STATUS_NULL_OUTPUT, &
       pb11_status_invalid_argument_base => PB11_STATUS_INVALID_ARGUMENT, &
       pb11_status_out_of_range_base => PB11_STATUS_OUT_OF_RANGE, &
       pb11_status_numerical_failure_base => PB11_STATUS_NUMERICAL_FAILURE, &
       pb11_status_exception_base => PB11_STATUS_EXCEPTION, &
       pb11_status_unknown_method_base => PB11_STATUS_UNKNOWN_METHOD, &
       fusion_pb11_3alpha_base => FUSION_PB11_3ALPHA, &
       fusion_dd_tp_base => FUSION_DD_TP, &
       fusion_dd_he3n_base => FUSION_DD_HE3N, &
       fusion_dt_alphan_base => FUSION_DT_ALPHAN, &
       fusion_dhe3_alphap_base => FUSION_DHE3_ALPHAP, &
       fusion_channel_count_base => FUSION_CHANNEL_COUNT, &
       fusion_thermal_birth_species_base => FUSION_THERMAL_BIRTH_SPECIES, &
       fusion_endpoint_s_base => FUSION_ENDPOINT_S, &
       fusion_high_flat_base => FUSION_HIGH_FLAT, &
       fusion_pb_low_tb_base => FUSION_PB_LOW_TB, &
       fusion_pb_low_ns_base => FUSION_PB_LOW_NS, &
       fusion_pb_low_c0_minus12_base => FUSION_PB_LOW_C0_MINUS12, &
       fusion_pb_low_c0_plus12_base => FUSION_PB_LOW_C0_PLUS12, &
       fusion_pb_remainder_entrance_proxy_base => &
            FUSION_PB_REMAINDER_ENTRANCE_PROXY, &
       fusion_pb_remainder_all_low_base => FUSION_PB_REMAINDER_ALL_LOW, &
       fusion_pb_remainder_all_broad_base => FUSION_PB_REMAINDER_ALL_BROAD
  implicit none
  private

  ! Re-export the common status and channel values used by the beam API.
  integer(c_int), parameter, public :: PB11_STATUS_OK = pb11_status_ok_base
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = &
       pb11_status_null_output_base
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = &
       pb11_status_invalid_argument_base
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = &
       pb11_status_out_of_range_base
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = &
       pb11_status_numerical_failure_base
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = &
       pb11_status_exception_base
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = &
       pb11_status_unknown_method_base

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = &
       fusion_pb11_3alpha_base
  integer(c_int), parameter, public :: FUSION_DD_TP = fusion_dd_tp_base
  integer(c_int), parameter, public :: FUSION_DD_HE3N = fusion_dd_he3n_base
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = fusion_dt_alphan_base
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = &
       fusion_dhe3_alphap_base
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = &
       fusion_channel_count_base
  integer(c_int), parameter, public :: FUSION_BEAM_BIRTH_SPECIES = &
       fusion_thermal_birth_species_base
  ! Keep the shared name available as well; the beam grid has the same
  ! seven-species product layout as the thermal-birth grid.
  integer(c_int), parameter, public :: FUSION_THERMAL_BIRTH_SPECIES = &
       fusion_thermal_birth_species_base

  integer(c_int), parameter, public :: FUSION_ENDPOINT_S = fusion_endpoint_s_base
  integer(c_int), parameter, public :: FUSION_HIGH_FLAT = fusion_high_flat_base
  integer(c_int), parameter, public :: FUSION_PB_LOW_TB = fusion_pb_low_tb_base
  integer(c_int), parameter, public :: FUSION_PB_LOW_NS = fusion_pb_low_ns_base
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_MINUS12 = &
       fusion_pb_low_c0_minus12_base
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_PLUS12 = &
       fusion_pb_low_c0_plus12_base

  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ENTRANCE_PROXY = &
       fusion_pb_remainder_entrance_proxy_base
  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ALL_LOW = &
       fusion_pb_remainder_all_low_base
  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ALL_BROAD = &
       fusion_pb_remainder_all_broad_base

  ! Exact C layout: eight consecutive c_double fields followed by nine
  ! c_int fields.  The trailing ABI padding is supplied by BIND(C).
  type, bind(C), public :: fusion_beam_birth_options_v1
     real(c_double) :: relative_max_J
     real(c_double) :: angular_max_exponent
     real(c_double) :: ground_state_q_J
     real(c_double) :: cutoff_J
     real(c_double) :: l1_fraction
     real(c_double) :: relative_phase
     real(c_double) :: narrow_peak_fraction
     real(c_double) :: continuum_peak_scale
     integer(c_int) :: continuation
     integer(c_int) :: pb_low
     integer(c_int) :: remainder_policy
     integer(c_int) :: broad_mode
     integer(c_int) :: fsci_policy
     integer(c_int) :: relative_order
     integer(c_int) :: angular_order
     integer(c_int) :: nq
     integer(c_int) :: ncos
  end type fusion_beam_birth_options_v1

  ! Exact C layout: the thermal-birth result followed by three c_double
  ! diagnostics.  fusion_thermal_birth_v1 is itself a BIND(C) type.
  type, bind(C), public :: fusion_beam_birth_v1
     type(fusion_thermal_birth_v1) :: spectrum
     real(c_double) :: relative_retained_probability
     real(c_double) :: retained_pair_probability
     real(c_double) :: angular_omitted_pair_probability
  end type fusion_beam_birth_v1

  public :: fusion_thermal_birth_v1
  public :: fusion_beam_birth_grid
  public :: fusion_beam_birth_resolve_support

  interface
     function c_fusion_beam_birth_grid(channel, projectile_slot, &
          projectile_energy_J, target_kT_J, options, cells, edges, birth, &
          out) bind(C, name="fusion_c_beam_birth_grid") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       integer(c_int), value :: projectile_slot
       real(c_double), value :: projectile_energy_J
       real(c_double), value :: target_kT_J
       type(c_ptr), value :: options
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_birth_grid
     function c_fusion_beam_birth_resolve_support(projectile_energy_J, &
          base, out) bind(C, name="fusion_c_beam_birth_resolve_support") &
          result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: projectile_energy_J
       type(c_ptr), value :: base, out
       integer(c_int) :: status
     end function c_fusion_beam_birth_resolve_support
  end interface

contains

  subroutine fusion_beam_birth_resolve_support(projectile_energy_J, base, &
       out, status)
    !! Resolve the beam relative-energy support through the C library.
    !! BASE and OUT are distinct Fortran actual arguments; the C API also
    !! supports aliasing, but Fortran INTENT(IN)/INTENT(OUT) arguments do not.
    real(c_double), intent(in) :: projectile_energy_J
    type(fusion_beam_birth_options_v1), intent(in), target :: base
    type(fusion_beam_birth_options_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    status = c_fusion_beam_birth_resolve_support(projectile_energy_J, &
         c_loc(base), c_loc(out))
  end subroutine fusion_beam_birth_resolve_support

  subroutine clear_thermal_birth_fields(out)
    type(fusion_thermal_birth_v1), intent(out) :: out

    out%reactivity_m3_s = 0.0_c_double
    out%reference_reactivity_m3_s = 0.0_c_double
    out%reactant_energy_moment_J_m3_s = 0.0_c_double
    out%reference_reactant_energy_moment_J_m3_s = 0.0_c_double
    out%product_energy_moment_J_m3_s = 0.0_c_double
    out%number_residual_m3_s = 0.0_c_double
    out%energy_residual_J_m3_s = 0.0_c_double
    out%relative_rate_discrepancy = 0.0_c_double
    out%relative_reactant_energy_discrepancy = 0.0_c_double
    out%cm_retained_probability = 0.0_c_double
    out%cm_tail_probability = 0.0_c_double
    out%cm_tail_energy_moment_J = 0.0_c_double
    out%max_cm_energy_shift_fraction = 0.0_c_double
    out%max_shell_remap_fraction = 0.0_c_double
    out%below_number_m3_s = 0.0_c_double
    out%below_energy_J_m3_s = 0.0_c_double
    out%above_number_m3_s = 0.0_c_double
    out%above_energy_J_m3_s = 0.0_c_double
  end subroutine clear_thermal_birth_fields

  subroutine clear_beam_birth(out)
    type(fusion_beam_birth_v1), intent(out) :: out

    call clear_thermal_birth_fields(out%spectrum)
    out%relative_retained_probability = 0.0_c_double
    out%retained_pair_probability = 0.0_c_double
    out%angular_omitted_pair_probability = 0.0_c_double
  end subroutine clear_beam_birth

  subroutine fusion_beam_birth_grid(channel, projectile_slot, projectile_energy_J, &
       target_kT_J, options, edges_J, birth, out, status)
    integer(c_int), intent(in) :: channel, projectile_slot
    real(c_double), intent(in) :: projectile_energy_J, target_kT_J
    type(fusion_beam_birth_options_v1), intent(in), target :: options
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(out), target, contiguous :: birth(:,:)
    type(fusion_beam_birth_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, expected_edges
    type(c_ptr) :: options_ptr, edges_ptr, birth_ptr, out_ptr

    call clear_beam_birth(out)
    birth = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Infer the C cell count from the first birth extent and validate every
    ! assumed-shape extent before taking C_LOC; malformed zero-size actual
    ! arrays remain safe.
    if (size(birth, 1, kind=c_int) < 1_c_int) return
    if (size(birth, 1, kind=c_int) >= huge(cells)) return
    cells = size(birth, 1, kind=c_int)
    if (size(birth, 2, kind=c_int) /= FUSION_BEAM_BIRTH_SPECIES) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return

    options_ptr = c_loc(options)
    edges_ptr = c_loc(edges_J(1))
    birth_ptr = c_loc(birth(1,1))
    out_ptr = c_loc(out)
    status = c_fusion_beam_birth_grid(channel, projectile_slot, &
         projectile_energy_J, target_kT_J, options_ptr, cells, edges_ptr, &
         birth_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       birth = 0.0_c_double
       call clear_beam_birth(out)
    end if
  end subroutine fusion_beam_birth_grid

end module fusion_beam_birth_fortran
