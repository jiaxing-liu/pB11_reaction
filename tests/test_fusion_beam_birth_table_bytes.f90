program test_fusion_beam_birth_table_bytes
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int8_t, c_size_t, &
       c_null_ptr, c_ptr
  use fusion_beam_birth_table_fortran, only : &
       fusion_beam_birth_options_v1, fusion_birth_table_control_v1, &
       FUSION_DT_ALPHAN, PB11_STATUS_OK, &
       fusion_beam_birth_table_create, fusion_beam_birth_table_destroy, &
       fusion_beam_birth_table_pack_size, fusion_beam_birth_table_pack, &
       fusion_beam_birth_table_unpack
  implicit none

  real(c_double), parameter :: kev = 1.602176634e-16_c_double
  real(c_double), parameter :: mev = 1.602176634e-13_c_double
  real(c_double), parameter :: projectile_energy = 2.1399859254691007e-13_c_double
  integer(c_int), parameter :: cells = 2_c_int
  integer(c_int8_t), parameter :: sentinel = -85_c_int8_t

  type(fusion_beam_birth_options_v1), target :: source
  type(fusion_birth_table_control_v1), target :: control
  real(c_double), target :: edges(cells + 1)
  type(c_ptr) :: table, restored, bad_table
  integer(c_int8_t), allocatable, target :: bytes(:), repacked(:), bad(:), empty(:)
  integer(c_int8_t), allocatable, target :: short_bytes(:)
  integer(c_size_t) :: required, packed_length, written
  integer(c_size_t) :: short_capacity, oversize_capacity
  integer(c_int) :: status
  integer :: n

  source%relative_max_J = 0.0_c_double
  source%angular_max_exponent = 0.0_c_double
  source%ground_state_q_J = 0.0_c_double
  source%cutoff_J = 0.0_c_double
  source%l1_fraction = 0.0_c_double
  source%relative_phase = 0.0_c_double
  source%narrow_peak_fraction = 0.0_c_double
  source%continuum_peak_scale = 0.0_c_double
  source%continuation = 0_c_int
  source%pb_low = 0_c_int
  source%remainder_policy = 0_c_int
  source%broad_mode = 0_c_int
  source%fsci_policy = 0_c_int
  source%relative_order = 0_c_int
  source%angular_order = 0_c_int
  source%nq = 0_c_int
  source%ncos = 0_c_int
  source%relative_max_J = 2.5_c_double * mev
  source%angular_max_exponent = 40.0_c_double
  source%ground_state_q_J = 91.84_c_double * kev
  source%cutoff_J = 0.001_c_double * mev
  source%l1_fraction = 0.76_c_double
  source%narrow_peak_fraction = 0.051_c_double
  source%continuum_peak_scale = 1.0_c_double
  source%continuation = 1_c_int
  source%broad_mode = 13_c_int
  source%relative_order = 16_c_int
  source%angular_order = 8_c_int
  source%nq = 4_c_int
  source%ncos = 4_c_int

  control%max_rate_error = 0.0_c_double
  control%max_debit_error = 0.0_c_double
  control%max_number_L1 = 0.0_c_double
  control%max_energy_L1 = 0.0_c_double
  control%max_direct_rate_discrepancy = 0.0_c_double
  control%max_direct_debit_discrepancy = 0.0_c_double
  control%max_knots = 0_c_int
  control%max_evaluations = 0_c_int
  control%max_depth = 0_c_int
  control%max_rate_error = 0.002_c_double
  control%max_debit_error = 0.002_c_double
  control%max_number_L1 = 0.01_c_double
  control%max_energy_L1 = 0.01_c_double
  control%max_direct_rate_discrepancy = 1.0e-4_c_double
  control%max_direct_debit_discrepancy = 1.0e-4_c_double
  control%max_knots = 512_c_int
  control%max_evaluations = 4096_c_int
  control%max_depth = 16_c_int

  edges = [0.0_c_double, 25.0_c_double * kev, 50.0_c_double * kev]
  table = c_null_ptr
  call fusion_beam_birth_table_create(FUSION_DT_ALPHAN, 1_c_int, &
       projectile_energy, 0.05_c_double * kev, 0.2_c_double * kev, source, &
       control, edges, table, status)
  call require(status == PB11_STATUS_OK, "DT table create")

  required = 99_c_size_t
  call fusion_beam_birth_table_pack_size(table, required, status)
  call require(status == PB11_STATUS_OK .and. required > 0_c_size_t, &
       "pack size")
  packed_length = required
  n = int(required)
  allocate(bytes(n), repacked(n), bad(n), short_bytes(n))

  bytes = sentinel
  call fusion_beam_birth_table_pack(table, bytes, required, written, status)
  call require(status == PB11_STATUS_OK .and. written == required, &
       "exact pack")
  repacked = sentinel
  restored = c_null_ptr
  call fusion_beam_birth_table_unpack(bytes, required, restored, status)
  call require(status == PB11_STATUS_OK, "exact unpack")
  call fusion_beam_birth_table_pack(restored, repacked, required, written, status)
  call require(status == PB11_STATUS_OK .and. written == required .and. &
       all(repacked == bytes), "exact byte roundtrip")

  short_bytes = sentinel
  short_capacity = required - 1_c_size_t
  call fusion_beam_birth_table_pack(table, short_bytes, short_capacity, written, &
       status)
  call require(status /= PB11_STATUS_OK .and. written == 0_c_size_t .and. &
       all(short_bytes == sentinel), "short capacity remains untouched")

  bytes = sentinel
  oversize_capacity = required + 1_c_size_t
  call fusion_beam_birth_table_pack(table, bytes, oversize_capacity, written, &
       status)
  call require(status /= PB11_STATUS_OK .and. written == 0_c_size_t .and. &
       all(bytes == sentinel), "oversize requested capacity rejected")

  required = 55_c_size_t
  call fusion_beam_birth_table_pack_size(c_null_ptr, required, status)
  call require(status /= PB11_STATUS_OK .and. required == 0_c_size_t, &
       "null table pack size clears required")
  written = 77_c_size_t
  bytes = sentinel
  call fusion_beam_birth_table_pack(c_null_ptr, bytes, int(size(bytes), c_size_t), &
       written, status)
  call require(status /= PB11_STATUS_OK .and. written == 0_c_size_t .and. &
       all(bytes == sentinel), "null table pack clears written")

  allocate(empty(0))
  bad_table = c_null_ptr
  call fusion_beam_birth_table_unpack(empty, 0_c_size_t, bad_table, status)
  call require(status /= PB11_STATUS_OK .and. .not. associated_c(bad_table), &
       "zero length input rejected")

  bad = repacked
  bad(1) = ieor(bad(1), 1_c_int8_t)
  bad_table = c_null_ptr
  call fusion_beam_birth_table_unpack(bad, packed_length, bad_table, status)
  call require(status /= PB11_STATUS_OK .and. .not. associated_c(bad_table), &
       "bad cache rejected")

  call fusion_beam_birth_table_destroy(restored)
  call fusion_beam_birth_table_destroy(table)
  print '(a)', "PASS: Fortran beam-table byte wrappers"

contains

  logical function associated_c(value)
    use, intrinsic :: iso_c_binding, only : c_associated
    type(c_ptr), intent(in) :: value
    associated_c = c_associated(value)
  end function associated_c

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
       print '(a,1x,a)', "FAIL:", label
       error stop 1
    end if
  end subroutine require

end program test_fusion_beam_birth_table_bytes
