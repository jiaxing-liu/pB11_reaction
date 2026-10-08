! Optional diagnostics ABI only; build independently from fusion_fortran.
! C uint64_t uses c_int64_t: callers must stay in the positive signed range.
! size_t lengths use c_size_t. Borrowed physics objects remain opaque c_ptrs.
module fusion_capture_fortran
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: fusion_capture_trial_input_v1, fusion_capture_replay_report_v1
  public :: fusion_capture_table_identity_v1
  public :: fusion_capture_write_trial_v1, fusion_capture_replay_file_v1
  public :: fusion_capture_context_create_v1, fusion_capture_context_destroy_v1
  public :: fusion_capture_thermal_unpack_v1, fusion_capture_beam_unpack_v1
  public :: fusion_capture_thermal_destroy_v1, fusion_capture_beam_destroy_v1
  public :: fusion_capture_table_identity_v1_get

  type, bind(C) :: fusion_capture_trial_input_v1
    integer(c_int) :: entry_point
    integer(c_int) :: original_status
    integer(c_int) :: host_zone
    real(c_double) :: host_time_s
    real(c_double) :: dt_s
    type(c_ptr) :: options
    type(c_ptr) :: fast_options
    type(c_ptr) :: thermal_tables
    integer(c_int) :: beam_table_count
    type(c_ptr) :: beam_tables
    integer(c_int) :: effective_charge
    integer(c_int) :: cells
    type(c_ptr) :: edges_J
    type(c_ptr) :: thermal_number_m3
    real(c_double) :: electron_energy_J_m3
    real(c_double) :: ion_energy_J_m3
    real(c_double) :: electron_density_m3
    type(c_ptr) :: thermal_charge_squared
    integer(c_int) :: inert_count
    type(c_ptr) :: inert
    type(c_ptr) :: coulomb_logs
    type(c_ptr) :: old_s_m3
    type(c_ptr) :: old_t_m3
    type(c_ptr) :: external_birth_m3_s
    type(c_ptr) :: escape_s_inv
    type(c_ptr) :: floor_limits
    integer(c_int) :: table_domain_policy
  end type

  type, bind(C) :: fusion_capture_replay_report_v1
    integer(c_int) :: original_status
    integer(c_int) :: replay_status
    integer(c_int) :: status_matches
    integer(c_int) :: entry_point
    integer(c_int) :: host_zone
    integer(c_int) :: cells
    real(c_double) :: host_time_s
    real(c_double) :: dt_s
    character(kind=c_char) :: kernel_identity(65)
    character(kind=c_char) :: output_sha256(65)
  end type

  type, bind(C) :: fusion_capture_table_identity_v1
    integer(c_int) :: kind
    integer(c_int64_t) :: packed_bytes
    character(kind=c_char) :: content_sha256(65)
    character(kind=c_char) :: kernel_identity(65)
    character(kind=c_char) :: qualified_path(4097)
  end type

  interface
    integer(c_int) function fusion_capture_write_trial_v1(context, path, input, max_bytes) bind(C)
      import :: c_int, c_ptr, c_char, c_int64_t, fusion_capture_trial_input_v1
      type(c_ptr), value :: context
      character(kind=c_char), intent(in) :: path(*)
      type(fusion_capture_trial_input_v1), intent(in) :: input
      integer(c_int64_t), value :: max_bytes
    end function

    integer(c_int) function fusion_capture_replay_file_v1(path, max_bytes, report) bind(C)
      import :: c_int, c_char, c_int64_t, fusion_capture_replay_report_v1
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int64_t), value :: max_bytes
      type(fusion_capture_replay_report_v1), intent(out) :: report
    end function

    integer(c_int) function fusion_capture_context_create_v1(capacity, out) bind(C)
      import :: c_int, c_int64_t, c_ptr
      integer(c_int64_t), value :: capacity
      type(c_ptr), intent(out) :: out
    end function

    ! The C context destroy returns void; it is deliberately a subroutine.
    subroutine fusion_capture_context_destroy_v1(context) bind(C)
      import :: c_ptr
      type(c_ptr), value :: context
    end subroutine

    integer(c_int) function fusion_capture_thermal_unpack_v1(context, bytes, length, qualified_path, &
        out, unpack_status) bind(C)
      import :: c_int, c_ptr, c_size_t, c_char
      type(c_ptr), value :: context
      type(c_ptr), value :: bytes
      integer(c_size_t), value :: length
      character(kind=c_char), intent(in) :: qualified_path(*)
      type(c_ptr), intent(out) :: out
      integer(c_int), intent(out) :: unpack_status
    end function

    integer(c_int) function fusion_capture_thermal_destroy_v1(context, handle) bind(C)
      import :: c_int, c_ptr
      type(c_ptr), value :: context
      type(c_ptr), intent(inout) :: handle
    end function

    integer(c_int) function fusion_capture_beam_unpack_v1(context, bytes, length, qualified_path, &
        out, unpack_status) bind(C)
      import :: c_int, c_ptr, c_size_t, c_char
      type(c_ptr), value :: context
      type(c_ptr), value :: bytes
      integer(c_size_t), value :: length
      character(kind=c_char), intent(in) :: qualified_path(*)
      type(c_ptr), intent(out) :: out
      integer(c_int), intent(out) :: unpack_status
    end function

    integer(c_int) function fusion_capture_beam_destroy_v1(context, handle) bind(C)
      import :: c_int, c_ptr
      type(c_ptr), value :: context
      type(c_ptr), intent(inout) :: handle
    end function

    integer(c_int) function fusion_capture_table_identity_v1_get(context, kind, handle, out) bind(C)
      import :: c_int, c_ptr, fusion_capture_table_identity_v1
      type(c_ptr), value :: context
      integer(c_int), value :: kind
      type(c_ptr), value :: handle
      type(fusion_capture_table_identity_v1), intent(out) :: out
    end function
  end interface
end module
