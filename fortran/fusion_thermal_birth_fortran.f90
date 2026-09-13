module fusion_thermal_birth_fortran
  !! ISO_C_BINDING wrapper for the finite-temperature product-birth grid API.
  !!
  !! The C API owns the reaction and kinematic model.  This module exposes
  !! the two C structs with their exact field kinds/order and validates the
  !! assumed-shape Fortran arrays before forming C pointers.  BIRTH is
  !! declared (cells,7); its column-major storage is the C species-major
  !! layout [species][cell] required by fusion_c_thermal_birth_grid.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc
  use pb11_fortran, only : pb11_status_ok_base => PB11_STATUS_OK, &
       pb11_status_null_output_base => PB11_STATUS_NULL_OUTPUT, &
       pb11_status_invalid_argument_base => PB11_STATUS_INVALID_ARGUMENT, &
       pb11_status_out_of_range_base => PB11_STATUS_OUT_OF_RANGE, &
       pb11_status_numerical_failure_base => PB11_STATUS_NUMERICAL_FAILURE, &
       pb11_status_exception_base => PB11_STATUS_EXCEPTION, &
       pb11_status_unknown_method_base => PB11_STATUS_UNKNOWN_METHOD
  implicit none
  private

  ! Re-export the common status values without importing unrelated wrappers.
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

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = 0_c_int
  integer(c_int), parameter, public :: FUSION_DD_TP = 1_c_int
  integer(c_int), parameter, public :: FUSION_DD_HE3N = 2_c_int
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = 3_c_int
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = 4_c_int
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = 5_c_int

  integer(c_int), parameter, public :: FUSION_THERMAL_BIRTH_SPECIES = 7_c_int

  integer(c_int), parameter, public :: FUSION_ENDPOINT_S = 1_c_int
  integer(c_int), parameter, public :: FUSION_HIGH_FLAT = 2_c_int
  integer(c_int), parameter, public :: FUSION_PB_LOW_TB = 0_c_int
  integer(c_int), parameter, public :: FUSION_PB_LOW_NS = 1_c_int
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_MINUS12 = 2_c_int
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_PLUS12 = 3_c_int

  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ENTRANCE_PROXY = 0_c_int
  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ALL_LOW = 1_c_int
  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ALL_BROAD = 2_c_int

  ! Exact C layout: eight consecutive c_double fields followed by nine
  ! c_int fields.  The trailing ABI padding is supplied by BIND(C).
  type, bind(C), public :: fusion_thermal_birth_options_v1
     real(c_double) :: relative_max_J
     real(c_double) :: cm_max_kT
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
     integer(c_int) :: cm_order
     integer(c_int) :: nq
     integer(c_int) :: ncos
  end type fusion_thermal_birth_options_v1

  ! Exact C layout: sixteen scalar/array double slots followed by four
  ! seven-element arrays (44 c_double values total).
  type, bind(C), public :: fusion_thermal_birth_v1
     real(c_double) :: reactivity_m3_s
     real(c_double) :: reference_reactivity_m3_s
     real(c_double) :: reactant_energy_moment_J_m3_s(2)
     real(c_double) :: reference_reactant_energy_moment_J_m3_s(2)
     real(c_double) :: product_energy_moment_J_m3_s
     real(c_double) :: number_residual_m3_s
     real(c_double) :: energy_residual_J_m3_s
     real(c_double) :: relative_rate_discrepancy
     real(c_double) :: relative_reactant_energy_discrepancy
     real(c_double) :: cm_retained_probability
     real(c_double) :: cm_tail_probability
     real(c_double) :: cm_tail_energy_moment_J
     real(c_double) :: max_cm_energy_shift_fraction
     real(c_double) :: max_shell_remap_fraction
     real(c_double) :: below_number_m3_s(7)
     real(c_double) :: below_energy_J_m3_s(7)
     real(c_double) :: above_number_m3_s(7)
     real(c_double) :: above_energy_J_m3_s(7)
  end type fusion_thermal_birth_v1

  public :: fusion_thermal_birth_grid

  interface
     function c_fusion_thermal_birth_grid(channel, kT_J, options, cells, &
          edges, birth, out) bind(C, name="fusion_c_thermal_birth_grid") &
          result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       real(c_double), value :: kT_J
       type(c_ptr), value :: options
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_thermal_birth_grid
  end interface

contains

  subroutine clear_thermal_birth(out)
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
  end subroutine clear_thermal_birth

  subroutine fusion_thermal_birth_grid(channel, kT_J, options, edges_J, &
       birth, out, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: kT_J
    type(fusion_thermal_birth_options_v1), intent(in), target :: options
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(out), target, contiguous :: birth(:,:)
    type(fusion_thermal_birth_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: cells, expected_edges
    type(c_ptr) :: options_ptr, edges_ptr, birth_ptr, out_ptr

    call clear_thermal_birth(out)
    birth = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! The C ABI uses a positive cell count and a fixed seven-species row.
    ! Check every extent before taking C_LOC; this also makes malformed
    ! zero-size actual arrays safe.
    if (size(birth, 1, kind=c_int) < 1_c_int) return
    if (size(birth, 1, kind=c_int) >= huge(cells)) return
    cells = size(birth, 1, kind=c_int)
    if (size(birth, 2, kind=c_int) /= FUSION_THERMAL_BIRTH_SPECIES) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return

    options_ptr = c_loc(options)
    edges_ptr = c_loc(edges_J(1))
    birth_ptr = c_loc(birth(1,1))
    out_ptr = c_loc(out)
    status = c_fusion_thermal_birth_grid(channel, kT_J, options_ptr, cells, &
         edges_ptr, birth_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       birth = 0.0_c_double
       call clear_thermal_birth(out)
    end if
  end subroutine fusion_thermal_birth_grid

end module fusion_thermal_birth_fortran
