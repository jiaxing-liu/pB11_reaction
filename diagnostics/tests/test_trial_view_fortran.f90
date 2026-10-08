! ABI field/layout checks only; this is not a serialization or physics test.
program test_trial_view_fortran
  use, intrinsic :: iso_c_binding
  use fusion_capture_fortran
  implicit none
  type(fusion_capture_trial_input_v1) :: view
  type(fusion_capture_replay_report_v1) :: report
  integer(c_int), target :: sentinels(14)
  type(c_ptr) :: pointers(14)
  integer(c_size_t) :: view_size
  integer(c_size_t) :: report_size
  integer(c_int) :: status
  integer :: i
  interface
    integer(c_int) function trial_view_probe(view, pointers, view_size, report, report_size) bind(C)
      import :: c_int, c_ptr, c_size_t, fusion_capture_trial_input_v1, fusion_capture_replay_report_v1
      type(fusion_capture_trial_input_v1), intent(in) :: view
      type(c_ptr), intent(in) :: pointers(*)
      integer(c_size_t), intent(out) :: view_size
      type(fusion_capture_replay_report_v1), intent(out) :: report
      integer(c_size_t), intent(out) :: report_size
    end function
  end interface
  do i = 1, 14
    sentinels(i) = i
    pointers(i) = c_loc(sentinels(i))
  end do
  view%entry_point = 101_c_int
  view%original_status = 102_c_int
  view%host_zone = 103_c_int
  view%host_time_s = 4.25_c_double
  view%dt_s = 5.25_c_double
  view%options = pointers(1)
  view%fast_options = pointers(2)
  view%thermal_tables = pointers(3)
  view%beam_table_count = 109_c_int
  view%beam_tables = pointers(4)
  view%effective_charge = 111_c_int
  view%cells = 112_c_int
  view%edges_J = pointers(5)
  view%thermal_number_m3 = pointers(6)
  view%electron_energy_J_m3 = 15.25_c_double
  view%ion_energy_J_m3 = 16.25_c_double
  view%electron_density_m3 = 17.25_c_double
  view%thermal_charge_squared = pointers(7)
  view%inert_count = 119_c_int
  view%inert = pointers(8)
  view%coulomb_logs = pointers(9)
  view%old_s_m3 = pointers(10)
  view%old_t_m3 = pointers(11)
  view%external_birth_m3_s = pointers(12)
  view%escape_s_inv = pointers(13)
  view%floor_limits = pointers(14)
  view%table_domain_policy = 127_c_int
  status = trial_view_probe(view, pointers, view_size, report, report_size)
  if (status /= 0_c_int) stop 1
  if (view_size /= c_sizeof(view)) stop 2
  if (report_size /= c_sizeof(report)) stop 3
  if (report%original_status /= 201_c_int) stop 11
  if (report%replay_status /= 202_c_int) stop 12
  if (report%status_matches /= 203_c_int) stop 13
  if (report%entry_point /= 204_c_int) stop 14
  if (report%host_zone /= 205_c_int) stop 15
  if (report%cells /= 206_c_int) stop 16
  if (report%host_time_s /= 7.75_c_double) stop 17
  if (report%dt_s /= 8.75_c_double) stop 18
  do i = 1, 64
    if (report%kernel_identity(i) /= 'K') stop 30
    if (report%output_sha256(i) /= 'a') stop 31
  end do
  if (report%kernel_identity(65) /= c_null_char) stop 32
  if (report%output_sha256(65) /= c_null_char) stop 33
  print *, 'trial view/report ABI PASS (no physics or serialization exercised)'
end program
