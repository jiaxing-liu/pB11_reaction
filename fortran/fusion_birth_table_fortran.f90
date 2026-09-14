module fusion_birth_table_fortran
  !! ISO_C_BINDING wrappers for the immutable thermal-birth coefficient table.
  !!
  !! A table handle is an opaque C pointer with caller-owned lifetime.  This
  !! module deliberately has no finalizer: successful handles are destroyed
  !! exactly once by fusion_birth_table_destroy.  Fortran birth arrays use
  !! shape (cells,7), whose column-major storage is the C species-major
  !! layout [species][cell].
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_null_ptr, c_ptr, c_size_t, c_loc
  use fusion_thermal_birth_fortran, only : &
       fusion_thermal_birth_options_v1, &
       thermal_status_ok => PB11_STATUS_OK, &
       thermal_status_null_output => PB11_STATUS_NULL_OUTPUT, &
       thermal_status_invalid_argument => PB11_STATUS_INVALID_ARGUMENT, &
       thermal_status_out_of_range => PB11_STATUS_OUT_OF_RANGE, &
       thermal_status_numerical_failure => PB11_STATUS_NUMERICAL_FAILURE, &
       thermal_status_exception => PB11_STATUS_EXCEPTION, &
       thermal_status_unknown_method => PB11_STATUS_UNKNOWN_METHOD, &
       thermal_birth_species => FUSION_THERMAL_BIRTH_SPECIES, &
       thermal_pb11_3alpha => FUSION_PB11_3ALPHA, &
       thermal_dd_tp => FUSION_DD_TP, thermal_dd_he3n => FUSION_DD_HE3N, &
       thermal_dt_alphan => FUSION_DT_ALPHAN, &
       thermal_dhe3_alphap => FUSION_DHE3_ALPHAP, &
       thermal_channel_count => FUSION_CHANNEL_COUNT
  implicit none
  private

  ! Re-export only the common status and channel constants needed by callers.
  integer(c_int), parameter, public :: PB11_STATUS_OK = thermal_status_ok
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = &
       thermal_status_null_output
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = &
       thermal_status_invalid_argument
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = &
       thermal_status_out_of_range
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = &
       thermal_status_numerical_failure
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = &
       thermal_status_exception
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = &
       thermal_status_unknown_method

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = &
       thermal_pb11_3alpha
  integer(c_int), parameter, public :: FUSION_DD_TP = thermal_dd_tp
  integer(c_int), parameter, public :: FUSION_DD_HE3N = thermal_dd_he3n
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = thermal_dt_alphan
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = &
       thermal_dhe3_alphap
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = &
       thermal_channel_count
  integer(c_int), parameter, public :: FUSION_BIRTH_TABLE_SPECIES = &
       thermal_birth_species

  ! Exact C layout: six c_double fields followed by three c_int fields and
  ! the ABI's trailing alignment padding.
  type, bind(C), public :: fusion_birth_table_control_v1
     real(c_double) :: max_rate_error
     real(c_double) :: max_debit_error
     real(c_double) :: max_number_L1
     real(c_double) :: max_energy_L1
     real(c_double) :: max_direct_rate_discrepancy
     real(c_double) :: max_direct_debit_discrepancy
     integer(c_int) :: max_knots
     integer(c_int) :: max_evaluations
     integer(c_int) :: max_depth
  end type fusion_birth_table_control_v1

  ! Exact C layout: eight doubles, four c_int fields, then the nested source
  ! and control structs in their declared C order.
  type, bind(C), public :: fusion_birth_table_info_v1
     real(c_double) :: lower_kT_J
     real(c_double) :: upper_kT_J
     real(c_double) :: max_validated_rate_error
     real(c_double) :: max_validated_debit_error
     real(c_double) :: max_validated_number_L1
     real(c_double) :: max_validated_energy_L1
     real(c_double) :: max_sampled_direct_rate_discrepancy
     real(c_double) :: max_sampled_direct_debit_discrepancy
     integer(c_int) :: channel
     integer(c_int) :: cells
     integer(c_int) :: knots
     integer(c_int) :: direct_evaluations
     type(fusion_thermal_birth_options_v1) :: source
     type(fusion_birth_table_control_v1) :: control
  end type fusion_birth_table_info_v1

  ! Exact C layout: one scalar, one two-element array, and four seven-element
  ! arrays (31 c_double values total).
  type, bind(C), public :: fusion_birth_coefficients_v1
     real(c_double) :: reactivity_m3_s
     real(c_double) :: reactant_energy_moment_J_m3_s(2)
     real(c_double) :: below_number_m3_s(7)
     real(c_double) :: below_energy_J_m3_s(7)
     real(c_double) :: above_number_m3_s(7)
     real(c_double) :: above_energy_J_m3_s(7)
  end type fusion_birth_coefficients_v1

  public :: fusion_birth_table_create
  public :: fusion_birth_table_destroy
  public :: fusion_birth_table_info
  public :: fusion_birth_table_evaluate

  interface
     function c_fusion_birth_table_create(channel, lower_kT_J, upper_kT_J, &
          source, control, cells, edges, out) bind(C, &
          name="fusion_c_birth_table_create") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       real(c_double), value :: lower_kT_J
       real(c_double), value :: upper_kT_J
       type(c_ptr), value :: source
       type(c_ptr), value :: control
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_birth_table_create

     subroutine c_fusion_birth_table_destroy(table) bind(C, &
          name="fusion_c_birth_table_destroy")
       import :: c_ptr
       type(c_ptr), value :: table
     end subroutine c_fusion_birth_table_destroy

     function c_fusion_birth_table_info(table, out) bind(C, &
          name="fusion_c_birth_table_info") result(status)
       import :: c_int, c_ptr
       type(c_ptr), value :: table
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_birth_table_info

     function c_fusion_birth_table_evaluate(table, kT_J, cells, birth, out) &
          bind(C, name="fusion_c_birth_table_evaluate") result(status)
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: table
       real(c_double), value :: kT_J
       integer(c_int), value :: cells
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_birth_table_evaluate
  end interface

contains

  subroutine clear_birth_table_info(info)
    type(fusion_birth_table_info_v1), intent(out) :: info

    info%lower_kT_J = 0.0_c_double
    info%upper_kT_J = 0.0_c_double
    info%max_validated_rate_error = 0.0_c_double
    info%max_validated_debit_error = 0.0_c_double
    info%max_validated_number_L1 = 0.0_c_double
    info%max_validated_energy_L1 = 0.0_c_double
    info%max_sampled_direct_rate_discrepancy = 0.0_c_double
    info%max_sampled_direct_debit_discrepancy = 0.0_c_double
    info%channel = 0_c_int
    info%cells = 0_c_int
    info%knots = 0_c_int
    info%direct_evaluations = 0_c_int
    info%source%relative_max_J = 0.0_c_double
    info%source%cm_max_kT = 0.0_c_double
    info%source%ground_state_q_J = 0.0_c_double
    info%source%cutoff_J = 0.0_c_double
    info%source%l1_fraction = 0.0_c_double
    info%source%relative_phase = 0.0_c_double
    info%source%narrow_peak_fraction = 0.0_c_double
    info%source%continuum_peak_scale = 0.0_c_double
    info%source%continuation = 0_c_int
    info%source%pb_low = 0_c_int
    info%source%remainder_policy = 0_c_int
    info%source%broad_mode = 0_c_int
    info%source%fsci_policy = 0_c_int
    info%source%relative_order = 0_c_int
    info%source%cm_order = 0_c_int
    info%source%nq = 0_c_int
    info%source%ncos = 0_c_int
    info%control%max_rate_error = 0.0_c_double
    info%control%max_debit_error = 0.0_c_double
    info%control%max_number_L1 = 0.0_c_double
    info%control%max_energy_L1 = 0.0_c_double
    info%control%max_direct_rate_discrepancy = 0.0_c_double
    info%control%max_direct_debit_discrepancy = 0.0_c_double
    info%control%max_knots = 0_c_int
    info%control%max_evaluations = 0_c_int
    info%control%max_depth = 0_c_int
  end subroutine clear_birth_table_info

  subroutine clear_birth_coefficients(out)
    type(fusion_birth_coefficients_v1), intent(out) :: out

    out%reactivity_m3_s = 0.0_c_double
    out%reactant_energy_moment_J_m3_s = 0.0_c_double
    out%below_number_m3_s = 0.0_c_double
    out%below_energy_J_m3_s = 0.0_c_double
    out%above_number_m3_s = 0.0_c_double
    out%above_energy_J_m3_s = 0.0_c_double
  end subroutine clear_birth_coefficients

  subroutine fusion_birth_table_create(channel, lower_kT_J, upper_kT_J, &
       source, control, edges_J, table, status)
    integer(c_int), intent(in) :: channel
    real(c_double), intent(in) :: lower_kT_J, upper_kT_J
    type(fusion_thermal_birth_options_v1), intent(in), target :: source
    type(fusion_birth_table_control_v1), intent(in), target :: control
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    type(c_ptr), intent(out), target :: table
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: edge_count
    integer(c_int) :: cells
    type(c_ptr), target :: created

    ! Assignment intentionally replaces a caller handle without destroying
    ! it.  Ownership remains explicit: callers destroy each live handle once.
    table = c_null_ptr
    created = c_null_ptr
    status = PB11_STATUS_INVALID_ARGUMENT

    edge_count = size(edges_J, kind=c_size_t)
    if (edge_count < 2_c_size_t .or. edge_count > 100001_c_size_t) return
    cells = int(edge_count - 1_c_size_t, c_int)

    status = c_fusion_birth_table_create(channel, lower_kT_J, upper_kT_J, &
         c_loc(source), c_loc(control), cells, c_loc(edges_J(1)), &
         c_loc(created))
    if (status /= PB11_STATUS_OK) then
       table = c_null_ptr
       return
    end if
    table = created
  end subroutine fusion_birth_table_create

  subroutine fusion_birth_table_destroy(table)
    type(c_ptr), intent(inout) :: table

    if (c_associated(table)) call c_fusion_birth_table_destroy(table)
    table = c_null_ptr
  end subroutine fusion_birth_table_destroy

  subroutine fusion_birth_table_info(table, info, status)
    type(c_ptr), intent(in) :: table
    type(fusion_birth_table_info_v1), intent(out), target :: info
    integer(c_int), intent(out) :: status

    call clear_birth_table_info(info)
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(table)) return

    status = c_fusion_birth_table_info(table, c_loc(info))
    if (status /= PB11_STATUS_OK) call clear_birth_table_info(info)
  end subroutine fusion_birth_table_info

  subroutine fusion_birth_table_evaluate(table, kT_J, birth, out, status)
    type(c_ptr), intent(in) :: table
    real(c_double), intent(in) :: kT_J
    real(c_double), intent(out), target, contiguous :: birth(:,:)
    type(fusion_birth_coefficients_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    type(fusion_birth_table_info_v1), target :: info
    integer(c_int) :: cells, query_status

    birth = 0.0_c_double
    call clear_birth_coefficients(out)
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(table)) return

    ! Query the immutable C-side extent before checking the Fortran array or
    ! forming C_LOC.  The caller cannot accidentally choose a C cell count.
    call fusion_birth_table_info(table, info, query_status)
    if (query_status /= PB11_STATUS_OK) then
       status = query_status
       return
    end if
    if (info%cells < 1_c_int .or. info%cells > 100000_c_int) return
    cells = info%cells
    if (size(birth, 1, kind=c_size_t) /= int(cells, c_size_t)) return
    if (size(birth, 2, kind=c_size_t) /= int(FUSION_BIRTH_TABLE_SPECIES, &
         c_size_t)) return

    status = c_fusion_birth_table_evaluate(table, kT_J, cells, &
         c_loc(birth(1,1)), c_loc(out))
    if (status /= PB11_STATUS_OK) then
       birth = 0.0_c_double
       call clear_birth_coefficients(out)
    end if
  end subroutine fusion_birth_table_evaluate

end module fusion_birth_table_fortran
