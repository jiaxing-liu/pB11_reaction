module fusion_beam_birth_table_fortran
  !! ISO_C_BINDING wrappers for the immutable fixed-energy beam-birth table.
  !!
  !! A table handle is an opaque C pointer with caller-owned lifetime.  This
  !! module deliberately has no finalizer: successful handles are destroyed
  !! exactly once by fusion_beam_birth_table_destroy.  Fortran birth arrays
  !! use shape (cells,7), whose column-major storage is the C species-major
  !! layout [species][cell].
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_int8_t, c_int64_t, c_null_ptr, c_ptr, c_size_t, c_loc
  use fusion_beam_birth_fortran, only : &
       fusion_beam_birth_options_v1, &
       beam_status_ok => PB11_STATUS_OK, &
       beam_status_null_output => PB11_STATUS_NULL_OUTPUT, &
       beam_status_invalid_argument => PB11_STATUS_INVALID_ARGUMENT, &
       beam_status_out_of_range => PB11_STATUS_OUT_OF_RANGE, &
       beam_status_numerical_failure => PB11_STATUS_NUMERICAL_FAILURE, &
       beam_status_exception => PB11_STATUS_EXCEPTION, &
       beam_status_unknown_method => PB11_STATUS_UNKNOWN_METHOD, &
       beam_pb11_3alpha => FUSION_PB11_3ALPHA, &
       beam_dd_tp => FUSION_DD_TP, beam_dd_he3n => FUSION_DD_HE3N, &
       beam_dt_alphan => FUSION_DT_ALPHAN, &
       beam_dhe3_alphap => FUSION_DHE3_ALPHAP, &
       beam_channel_count => FUSION_CHANNEL_COUNT, &
       beam_species => FUSION_BEAM_BIRTH_SPECIES, &
       beam_endpoint_s => FUSION_ENDPOINT_S, &
       beam_high_flat => FUSION_HIGH_FLAT, &
       beam_pb_low_tb => FUSION_PB_LOW_TB, &
       beam_pb_low_ns => FUSION_PB_LOW_NS, &
       beam_pb_low_c0_minus12 => FUSION_PB_LOW_C0_MINUS12, &
       beam_pb_low_c0_plus12 => FUSION_PB_LOW_C0_PLUS12, &
       beam_remainder_entrance_proxy => FUSION_PB_REMAINDER_ENTRANCE_PROXY, &
       beam_remainder_all_low => FUSION_PB_REMAINDER_ALL_LOW, &
       beam_remainder_all_broad => FUSION_PB_REMAINDER_ALL_BROAD
  use fusion_birth_table_fortran, only : fusion_birth_table_control_v1, &
       fusion_birth_coefficients_v1
  implicit none
  private

  ! Re-export the common status and channel values used by the table API.
  integer(c_int), parameter, public :: PB11_STATUS_OK = beam_status_ok
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = &
       beam_status_null_output
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = &
       beam_status_invalid_argument
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = &
       beam_status_out_of_range
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = &
       beam_status_numerical_failure
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = &
       beam_status_exception
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = &
       beam_status_unknown_method

  integer(c_int), parameter, public :: FUSION_PB11_3ALPHA = beam_pb11_3alpha
  integer(c_int), parameter, public :: FUSION_DD_TP = beam_dd_tp
  integer(c_int), parameter, public :: FUSION_DD_HE3N = beam_dd_he3n
  integer(c_int), parameter, public :: FUSION_DT_ALPHAN = beam_dt_alphan
  integer(c_int), parameter, public :: FUSION_DHE3_ALPHAP = beam_dhe3_alphap
  integer(c_int), parameter, public :: FUSION_CHANNEL_COUNT = beam_channel_count
  integer(c_int), parameter, public :: FUSION_BEAM_BIRTH_SPECIES = beam_species
  integer(c_int), parameter, public :: FUSION_BEAM_BIRTH_TABLE_SPECIES = &
       beam_species

  integer(c_int), parameter, public :: FUSION_ENDPOINT_S = beam_endpoint_s
  integer(c_int), parameter, public :: FUSION_HIGH_FLAT = beam_high_flat
  integer(c_int), parameter, public :: FUSION_PB_LOW_TB = beam_pb_low_tb
  integer(c_int), parameter, public :: FUSION_PB_LOW_NS = beam_pb_low_ns
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_MINUS12 = &
       beam_pb_low_c0_minus12
  integer(c_int), parameter, public :: FUSION_PB_LOW_C0_PLUS12 = &
       beam_pb_low_c0_plus12
  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ENTRANCE_PROXY = &
       beam_remainder_entrance_proxy
  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ALL_LOW = &
       beam_remainder_all_low
  integer(c_int), parameter, public :: FUSION_PB_REMAINDER_ALL_BROAD = &
       beam_remainder_all_broad

  ! Exact C layout: nine doubles, five c_int fields, two uint64_t counters,
  ! followed by the nested source and control structs in their C order.
  type, bind(C), public :: fusion_beam_birth_table_info_v1
     real(c_double) :: projectile_energy_J
     real(c_double) :: lower_kT_J
     real(c_double) :: upper_kT_J
     real(c_double) :: max_validated_rate_error
     real(c_double) :: max_validated_debit_error
     real(c_double) :: max_validated_number_L1
     real(c_double) :: max_validated_energy_L1
     real(c_double) :: max_sampled_direct_rate_discrepancy
     real(c_double) :: max_sampled_direct_debit_discrepancy
     integer(c_int) :: channel
     integer(c_int) :: projectile_slot
     integer(c_int) :: cells
     integer(c_int) :: knots
     integer(c_int) :: direct_evaluations
     integer(c_int64_t) :: spectral_entries_evaluated
     integer(c_int64_t) :: stored_spectral_entries
     type(fusion_beam_birth_options_v1) :: source
     type(fusion_birth_table_control_v1) :: control
  end type fusion_beam_birth_table_info_v1

  ! The coefficient layout is shared with fusion_birth_table_fortran.  Keep
  ! that type use-associated rather than duplicating the C struct here.
  public :: fusion_beam_birth_options_v1
  public :: fusion_birth_table_control_v1
  public :: fusion_birth_coefficients_v1

  public :: fusion_beam_birth_table_matches_request
  public :: fusion_beam_birth_table_create
  public :: fusion_beam_birth_table_destroy
  public :: fusion_beam_birth_table_info
  public :: fusion_beam_birth_table_evaluate
  public :: fusion_beam_birth_table_pack_size
  public :: fusion_beam_birth_table_pack
  public :: fusion_beam_birth_table_unpack

  interface
    function c_beam_matches(table,channel,slot,energy,lower,upper,source,control,cells,edges,matches) &
        bind(C,name="fusion_c_beam_birth_table_matches_request") result(status)
      import :: c_ptr,c_int,c_double
      type(c_ptr),value :: table,source,control,edges,matches
      integer(c_int),value :: channel,slot,cells
      real(c_double),value :: energy,lower,upper
      integer(c_int) :: status
    end function
     function c_fusion_beam_birth_table_create(channel, projectile_slot, &
          projectile_energy_J, lower_kT_J, upper_kT_J, source, control, &
          cells, edges, out) bind(C, &
          name="fusion_c_beam_birth_table_create") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: channel
       integer(c_int), value :: projectile_slot
       real(c_double), value :: projectile_energy_J
       real(c_double), value :: lower_kT_J
       real(c_double), value :: upper_kT_J
       type(c_ptr), value :: source
       type(c_ptr), value :: control
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_birth_table_create

     subroutine c_fusion_beam_birth_table_destroy(table) bind(C, &
          name="fusion_c_beam_birth_table_destroy")
       import :: c_ptr
       type(c_ptr), value :: table
     end subroutine c_fusion_beam_birth_table_destroy

     function c_fusion_beam_birth_table_info(table, out) bind(C, &
          name="fusion_c_beam_birth_table_info") result(status)
       import :: c_int, c_ptr
       type(c_ptr), value :: table
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_birth_table_info

     function c_fusion_beam_birth_table_evaluate(table, target_kT_J, cells, &
          birth, out) bind(C, &
          name="fusion_c_beam_birth_table_evaluate") result(status)
       import :: c_double, c_int, c_ptr
       type(c_ptr), value :: table
       real(c_double), value :: target_kT_J
       integer(c_int), value :: cells
       type(c_ptr), value :: birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_birth_table_evaluate

     function c_fusion_beam_birth_table_pack_size(table, required) bind(C, &
          name="fusion_c_beam_birth_table_pack_size") result(status)
       import :: c_int, c_ptr, c_size_t
       type(c_ptr), value :: table
       type(c_ptr), value :: required
       integer(c_int) :: status
     end function c_fusion_beam_birth_table_pack_size

     function c_fusion_beam_birth_table_pack(table, buffer, capacity, written) &
          bind(C, name="fusion_c_beam_birth_table_pack") result(status)
       import :: c_int, c_ptr, c_size_t
       type(c_ptr), value :: table
       type(c_ptr), value :: buffer
       integer(c_size_t), value :: capacity
       type(c_ptr), value :: written
       integer(c_int) :: status
     end function c_fusion_beam_birth_table_pack

     function c_fusion_beam_birth_table_unpack(buffer, length, out) bind(C, &
          name="fusion_c_beam_birth_table_unpack") result(status)
       import :: c_int, c_ptr, c_size_t
       type(c_ptr), value :: buffer
       integer(c_size_t), value :: length
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_beam_birth_table_unpack
  end interface

contains

  subroutine fusion_beam_birth_table_matches_request(table,channel,slot,energy,lower,upper, &
       source,control,edges,matches,status)
    type(c_ptr),intent(in) :: table
    integer(c_int),intent(in) :: channel,slot
    real(c_double),intent(in) :: energy,lower,upper
    type(fusion_beam_birth_options_v1),intent(in),target :: source
    type(fusion_birth_table_control_v1),intent(in),target :: control
    real(c_double),intent(in),contiguous,target :: edges(:)
    logical,intent(out) :: matches
    integer(c_int),intent(out) :: status
    integer(c_int),target :: answer
    integer(c_size_t) :: cells
    matches=.false.
    status=PB11_STATUS_INVALID_ARGUMENT
    cells=size(edges,kind=c_size_t)-1_c_size_t
    if(cells<1_c_size_t.or.cells>100000_c_size_t)return
    answer=0
    status=c_beam_matches(table,channel,slot,energy,lower,upper,c_loc(source),c_loc(control), &
      int(cells,c_int),c_loc(edges(1)),c_loc(answer))
    if(status==PB11_STATUS_OK)matches=answer==1
  end subroutine


  subroutine clear_beam_birth_table_info(info)
    type(fusion_beam_birth_table_info_v1), intent(out) :: info

    info%projectile_energy_J = 0.0_c_double
    info%lower_kT_J = 0.0_c_double
    info%upper_kT_J = 0.0_c_double
    info%max_validated_rate_error = 0.0_c_double
    info%max_validated_debit_error = 0.0_c_double
    info%max_validated_number_L1 = 0.0_c_double
    info%max_validated_energy_L1 = 0.0_c_double
    info%max_sampled_direct_rate_discrepancy = 0.0_c_double
    info%max_sampled_direct_debit_discrepancy = 0.0_c_double
    info%channel = 0_c_int
    info%projectile_slot = 0_c_int
    info%cells = 0_c_int
    info%knots = 0_c_int
    info%direct_evaluations = 0_c_int
    info%spectral_entries_evaluated = 0_c_int64_t
    info%stored_spectral_entries = 0_c_int64_t
    info%source%relative_max_J = 0.0_c_double
    info%source%angular_max_exponent = 0.0_c_double
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
    info%source%angular_order = 0_c_int
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
  end subroutine clear_beam_birth_table_info

  subroutine clear_birth_coefficients(out)
    type(fusion_birth_coefficients_v1), intent(out) :: out

    out%reactivity_m3_s = 0.0_c_double
    out%reactant_energy_moment_J_m3_s = 0.0_c_double
    out%below_number_m3_s = 0.0_c_double
    out%below_energy_J_m3_s = 0.0_c_double
    out%above_number_m3_s = 0.0_c_double
    out%above_energy_J_m3_s = 0.0_c_double
  end subroutine clear_birth_coefficients

  subroutine fusion_beam_birth_table_create(channel, projectile_slot, &
       projectile_energy_J, lower_kT_J, upper_kT_J, source, control, edges_J, &
       table, status)
    integer(c_int), intent(in) :: channel, projectile_slot
    real(c_double), intent(in) :: projectile_energy_J, lower_kT_J, upper_kT_J
    type(fusion_beam_birth_options_v1), intent(in), target :: source
    type(fusion_birth_table_control_v1), intent(in), target :: control
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    type(c_ptr), intent(out), target :: table
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: edge_count
    integer(c_int) :: cells
    type(c_ptr), target :: created

    ! Replacing a caller handle is intentional; ownership remains explicit.
    table = c_null_ptr
    created = c_null_ptr
    status = PB11_STATUS_INVALID_ARGUMENT

    edge_count = size(edges_J, kind=c_size_t)
    if (edge_count < 2_c_size_t .or. edge_count > 100001_c_size_t) return
    cells = int(edge_count - 1_c_size_t, c_int)

    status = c_fusion_beam_birth_table_create(channel, projectile_slot, &
         projectile_energy_J, lower_kT_J, upper_kT_J, c_loc(source), &
         c_loc(control), cells, c_loc(edges_J(1)), c_loc(created))
    if (status /= PB11_STATUS_OK) then
       table = c_null_ptr
       return
    end if
    table = created
  end subroutine fusion_beam_birth_table_create

  subroutine fusion_beam_birth_table_destroy(table)
    type(c_ptr), intent(inout) :: table

    if (c_associated(table)) call c_fusion_beam_birth_table_destroy(table)
    table = c_null_ptr
  end subroutine fusion_beam_birth_table_destroy

  subroutine fusion_beam_birth_table_info(table, info, status)
    type(c_ptr), intent(in) :: table
    type(fusion_beam_birth_table_info_v1), intent(out), target :: info
    integer(c_int), intent(out) :: status

    call clear_beam_birth_table_info(info)
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(table)) return

    status = c_fusion_beam_birth_table_info(table, c_loc(info))
    if (status /= PB11_STATUS_OK) call clear_beam_birth_table_info(info)
  end subroutine fusion_beam_birth_table_info

  subroutine fusion_beam_birth_table_evaluate(table, target_kT_J, birth, out, &
       status)
    type(c_ptr), intent(in) :: table
    real(c_double), intent(in) :: target_kT_J
    real(c_double), intent(out), target, contiguous :: birth(:,:)
    type(fusion_birth_coefficients_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    type(fusion_beam_birth_table_info_v1), target :: info
    integer(c_int) :: cells, query_status

    birth = 0.0_c_double
    call clear_birth_coefficients(out)
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(table)) return

    ! Obtain the immutable C-side extent before checking the Fortran array or
    ! forming C_LOC.  The caller cannot choose a mismatched C cell count.
    call fusion_beam_birth_table_info(table, info, query_status)
    if (query_status /= PB11_STATUS_OK) then
       status = query_status
       return
    end if
    if (info%cells < 1_c_int .or. info%cells > 100000_c_int) return
    cells = info%cells
    if (size(birth, 1, kind=c_size_t) /= int(cells, c_size_t)) return
    if (size(birth, 2, kind=c_size_t) /= int(FUSION_BEAM_BIRTH_TABLE_SPECIES, &
         c_size_t)) return

    status = c_fusion_beam_birth_table_evaluate(table, target_kT_J, cells, &
         c_loc(birth(1,1)), c_loc(out))
    if (status /= PB11_STATUS_OK) then
       birth = 0.0_c_double
       call clear_birth_coefficients(out)
    end if
  end subroutine fusion_beam_birth_table_evaluate

  subroutine fusion_beam_birth_table_pack_size(table, required, status)
    type(c_ptr), intent(in) :: table
    integer(c_size_t), intent(out), target :: required
    integer(c_int), intent(out) :: status

    required = 0_c_size_t
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(table)) return

    status = c_fusion_beam_birth_table_pack_size(table, c_loc(required))
    if (status /= PB11_STATUS_OK) required = 0_c_size_t
  end subroutine fusion_beam_birth_table_pack_size

  subroutine fusion_beam_birth_table_pack(table, buffer, capacity, &
       bytes_written, status)
    type(c_ptr), intent(in) :: table
    integer(c_int8_t), intent(inout), target, contiguous :: buffer(:)
    integer(c_size_t), intent(in) :: capacity
    integer(c_size_t), intent(out), target :: bytes_written
    integer(c_int), intent(out) :: status

    bytes_written = 0_c_size_t
    status = PB11_STATUS_INVALID_ARGUMENT
    if (.not. c_associated(table)) return
    ! C requires a nonnull buffer even when capacity is zero.  Avoid C_LOC for
    ! zero-length Fortran actuals and reject a capacity beyond the actual
    ! caller-owned extent before entering the ABI.
    if (capacity <= 0_c_size_t) return
    if (size(buffer, kind=c_size_t) < capacity) return

    status = c_fusion_beam_birth_table_pack(table, c_loc(buffer(1)), &
         capacity, c_loc(bytes_written))
    if (status /= PB11_STATUS_OK) bytes_written = 0_c_size_t
  end subroutine fusion_beam_birth_table_pack

  subroutine fusion_beam_birth_table_unpack(buffer, length, table, status)
    integer(c_int8_t), intent(in), target, contiguous :: buffer(:)
    integer(c_size_t), intent(in) :: length
    type(c_ptr), intent(out), target :: table
    integer(c_int), intent(out) :: status

    type(c_ptr), target :: unpacked

    table = c_null_ptr
    unpacked = c_null_ptr
    status = PB11_STATUS_INVALID_ARGUMENT
    ! Do not form C_LOC for a zero-length actual or for a requested length
    ! beyond the actual buffer.  The C ABI receives exactly the requested
    ! prefix when a larger Fortran buffer is supplied.
    if (length <= 0_c_size_t) return
    if (size(buffer, kind=c_size_t) < length) return

    status = c_fusion_beam_birth_table_unpack(c_loc(buffer(1)), length, &
         c_loc(unpacked))
    if (status /= PB11_STATUS_OK) then
       table = c_null_ptr
       return
    end if
    table = unpacked
  end subroutine fusion_beam_birth_table_unpack

end module fusion_beam_birth_table_fortran
