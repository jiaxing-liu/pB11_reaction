program test_fusion_birth_table_bytes
  use, intrinsic :: iso_c_binding, only : c_associated, c_double, c_int, &
       c_int8_t, c_null_ptr, c_ptr, c_size_t
  use fusion_thermal_birth_fortran, only : &
       fusion_thermal_birth_options_v1, FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, &
       FUSION_PB_LOW_TB, FUSION_PB_REMAINDER_ENTRANCE_PROXY
  use fusion_birth_table_fortran, only : &
       PB11_STATUS_OK, fusion_birth_table_control_v1, &
       fusion_birth_table_create, fusion_birth_table_destroy, &
       fusion_birth_table_kernel_identity, fusion_birth_table_matches_request, &
       fusion_birth_table_pack_size, fusion_birth_table_pack, &
       fusion_birth_table_unpack
  implicit none

  real(c_double), parameter :: kev_j = 1.602176634e-16_c_double
  real(c_double), parameter :: mev_j = 1.602176634e-13_c_double
  real(c_double), parameter :: lower_kT_j = 9.0_c_double * kev_j
  real(c_double), parameter :: upper_kT_j = 11.0_c_double * kev_j
  integer(c_int), parameter :: cells = 16_c_int
  integer(c_int8_t), parameter :: sentinel = -85_c_int8_t

  type(fusion_thermal_birth_options_v1), target :: source, altered_source
  type(fusion_birth_table_control_v1), target :: control, altered_control
  real(c_double), target :: edges(cells + 1_c_int), altered_edges(cells + 1_c_int)
  type(c_ptr) :: table, restored, bad_table
  integer(c_int8_t), allocatable, target :: bytes(:), repacked(:), corrupt(:)
  integer(c_int8_t), allocatable, target :: short_bytes(:), empty(:)
  integer(c_size_t) :: required, packed_length, written
  integer(c_size_t) :: short_capacity
  integer(c_int) :: status
  logical :: matches
  character(len=64) :: identity
  integer :: failures, n

  failures = 0
  call make_source(source)
  call make_control(control)
  call make_edges(edges)
  table = c_null_ptr
  restored = c_null_ptr
  bad_table = c_null_ptr

  identity = fusion_birth_table_kernel_identity()
  call check(len_trim(identity) == 64, &
       'kernel identity exposes the complete 64-character string')

  call fusion_birth_table_create(FUSION_DT_ALPHAN, lower_kT_j, upper_kT_j, &
       source, control, edges, table, status)
  call check(status == PB11_STATUS_OK .and. c_associated(table), &
       'small real-DT table construction succeeds')
  if (.not. c_associated(table)) then
     write(*, '(A, I0)') 'create status: ', status
     error stop 1
  end if

  required = 0_c_size_t
  call fusion_birth_table_pack_size(table, required, status)
  call check(status == PB11_STATUS_OK .and. required > 0_c_size_t, &
       'pack size returns a positive byte count')
  if (status /= PB11_STATUS_OK .or. required == 0_c_size_t) then
     call fusion_birth_table_destroy(table)
     error stop 1
  end if
  packed_length = required
  n = int(required)
  allocate(bytes(n), repacked(n), corrupt(n), short_bytes(n))

  bytes = sentinel
  call fusion_birth_table_pack(table, bytes, required, written, status)
  call check(status == PB11_STATUS_OK .and. written == required, &
       'exact buffer pack succeeds')

  repacked = sentinel
  call fusion_birth_table_unpack(bytes, packed_length, restored, status)
  call check(status == PB11_STATUS_OK .and. c_associated(restored), &
       'exact buffer unpack succeeds')
  if (c_associated(restored)) then
     call fusion_birth_table_pack(restored, repacked, required, written, status)
     call check(status == PB11_STATUS_OK .and. written == required .and. &
          all(repacked == bytes), 'unpack and repack preserve exact bytes')
  end if

  call fusion_birth_table_matches_request(restored, FUSION_DT_ALPHAN, &
       lower_kT_j, upper_kT_j, source, control, edges, matches, status)
  call check(status == PB11_STATUS_OK .and. matches, &
       'full constructor request matches after unpack')

  altered_source = source
  altered_source%cm_max_kT = source%cm_max_kT - 1.0_c_double
  call fusion_birth_table_matches_request(restored, FUSION_DT_ALPHAN, &
       lower_kT_j, upper_kT_j, altered_source, control, edges, matches, status)
  call check(status == PB11_STATUS_OK .and. .not. matches, &
       'source option mismatch is reported without an error')

  altered_control = control
  altered_control%max_rate_error = control%max_rate_error * 0.9_c_double
  call fusion_birth_table_matches_request(restored, FUSION_DT_ALPHAN, &
       lower_kT_j, upper_kT_j, source, altered_control, edges, matches, status)
  call check(status == PB11_STATUS_OK .and. .not. matches, &
       'construction control mismatch is reported without an error')

  altered_edges = edges
  altered_edges(2) = nearest(altered_edges(2), 1.0_c_double)
  call fusion_birth_table_matches_request(restored, FUSION_DT_ALPHAN, &
       lower_kT_j, upper_kT_j, source, control, altered_edges, matches, status)
  call check(status == PB11_STATUS_OK .and. .not. matches, &
       'edge mismatch is reported without an error')

  short_bytes = sentinel
  short_capacity = required - 1_c_size_t
  call fusion_birth_table_pack(table, short_bytes, short_capacity, written, status)
  call check(status /= PB11_STATUS_OK .and. written == 0_c_size_t .and. &
       all(short_bytes == sentinel), &
       'short pack capacity leaves caller bytes and count unchanged')

  required = 77_c_size_t
  call fusion_birth_table_pack_size(c_null_ptr, required, status)
  call check(status /= PB11_STATUS_OK .and. required == 0_c_size_t, &
       'null table pack size clears the required count')

  allocate(empty(0))
  bad_table = c_null_ptr
  call fusion_birth_table_unpack(empty, 0_c_size_t, bad_table, status)
  call check(status /= PB11_STATUS_OK .and. .not. c_associated(bad_table), &
       'zero-length buffer is rejected without forming C_LOC')

  corrupt = bytes
  corrupt(1) = ieor(corrupt(1), 1_c_int8_t)
  bad_table = c_null_ptr
  call fusion_birth_table_unpack(corrupt, packed_length, bad_table, status)
  call check(status /= PB11_STATUS_OK .and. .not. c_associated(bad_table), &
       'corrupt buffer is rejected atomically')

  call fusion_birth_table_destroy(restored)
  call fusion_birth_table_destroy(table)
  call check(.not. c_associated(restored) .and. .not. c_associated(table), &
       'destroy clears both caller handles')

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran thermal byte test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'Fortran thermal birth-table byte test passed'

contains

  subroutine make_source(o)
    type(fusion_thermal_birth_options_v1), intent(out) :: o

    o%relative_max_J = 5.0_c_double * mev_j
    o%cm_max_kT = 40.0_c_double
    o%ground_state_q_J = 91.84_c_double * kev_j
    o%cutoff_J = 0.001_c_double * mev_j
    o%l1_fraction = 0.76_c_double
    o%relative_phase = 0.0_c_double
    o%narrow_peak_fraction = 0.051_c_double
    o%continuum_peak_scale = 1.0_c_double
    o%continuation = FUSION_ENDPOINT_S
    o%pb_low = FUSION_PB_LOW_TB
    o%remainder_policy = FUSION_PB_REMAINDER_ENTRANCE_PROXY
    o%broad_mode = 13_c_int
    o%fsci_policy = 0_c_int
    o%relative_order = 16_c_int
    o%cm_order = 12_c_int
    o%nq = 4_c_int
    o%ncos = 4_c_int
  end subroutine make_source

  subroutine make_control(o)
    type(fusion_birth_table_control_v1), intent(out) :: o

    o%max_rate_error = 1.0e-2_c_double
    o%max_debit_error = 1.0e-2_c_double
    o%max_number_L1 = 1.0e-2_c_double
    o%max_energy_L1 = 1.0e-2_c_double
    o%max_direct_rate_discrepancy = 1.0e-5_c_double
    o%max_direct_debit_discrepancy = 1.0e-5_c_double
    o%max_knots = 32_c_int
    o%max_evaluations = 256_c_int
    o%max_depth = 10_c_int
  end subroutine make_control

  subroutine make_edges(e)
    real(c_double), intent(out) :: e(:)
    integer :: i

    do i = 1, size(e)
       e(i) = 25.0_c_double * mev_j * real(i - 1, c_double) / &
            real(size(e) - 1, c_double)
    end do
  end subroutine make_edges

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

end program test_fusion_birth_table_bytes
