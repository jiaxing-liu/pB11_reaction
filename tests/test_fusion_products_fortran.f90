program test_fusion_products_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, &
       ieee_value
  use fusion_products_fortran
  implicit none

  integer :: failures

  failures = 0
  call verify_mapping_and_spill()
  call verify_noncontiguous_mapping()
  call verify_empty_packets()
  call verify_three_body_conservation()
  call verify_errors()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran fusion-products test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran fusion-products tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function close_relative(actual, expected, tolerance)
    real(c_double), intent(in) :: actual, expected, tolerance
    real(c_double) :: scale

    if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
       close_relative = .false.
       return
    end if
    if (actual == expected) then
       close_relative = .true.
       return
    end if
    scale = max(abs(actual), abs(expected))
    close_relative = scale > 0.0_c_double .and. &
         abs(actual - expected) <= tolerance * scale
  end function close_relative

  subroutine fill_mapping(out, value)
    type(fusion_birth_mapping_v1), intent(out) :: out
    real(c_double), intent(in) :: value

    out%input_number_m3_s = value
    out%input_energy_W_m3 = value
    out%mapped_number_m3_s = value
    out%mapped_energy_W_m3 = value
    out%below_number_m3_s = value
    out%below_energy_W_m3 = value
    out%above_number_m3_s = value
    out%above_energy_W_m3 = value
    out%number_residual_m3_s = value
    out%energy_residual_W_m3 = value
  end subroutine fill_mapping

  logical function mapping_is_zero(out)
    type(fusion_birth_mapping_v1), intent(in) :: out

    mapping_is_zero = out%input_number_m3_s == 0.0_c_double .and. &
         out%input_energy_W_m3 == 0.0_c_double .and. &
         out%mapped_number_m3_s == 0.0_c_double .and. &
         out%mapped_energy_W_m3 == 0.0_c_double .and. &
         out%below_number_m3_s == 0.0_c_double .and. &
         out%below_energy_W_m3 == 0.0_c_double .and. &
         out%above_number_m3_s == 0.0_c_double .and. &
         out%above_energy_W_m3 == 0.0_c_double .and. &
         out%number_residual_m3_s == 0.0_c_double .and. &
         out%energy_residual_W_m3 == 0.0_c_double
  end function mapping_is_zero

  logical function mappings_match(left, right, tolerance)
    type(fusion_birth_mapping_v1), intent(in) :: left, right
    real(c_double), intent(in) :: tolerance

    mappings_match = close_relative(left%input_number_m3_s, &
         right%input_number_m3_s, tolerance) .and. &
         close_relative(left%input_energy_W_m3, right%input_energy_W_m3, &
         tolerance) .and. &
         close_relative(left%mapped_number_m3_s, right%mapped_number_m3_s, &
         tolerance) .and. &
         close_relative(left%mapped_energy_W_m3, right%mapped_energy_W_m3, &
         tolerance) .and. &
         close_relative(left%below_number_m3_s, right%below_number_m3_s, &
         tolerance) .and. &
         close_relative(left%below_energy_W_m3, right%below_energy_W_m3, &
         tolerance) .and. &
         close_relative(left%above_number_m3_s, right%above_number_m3_s, &
         tolerance) .and. &
         close_relative(left%above_energy_W_m3, right%above_energy_W_m3, &
         tolerance) .and. &
         close_relative(left%number_residual_m3_s, right%number_residual_m3_s, &
         tolerance) .and. &
         close_relative(left%energy_residual_W_m3, right%energy_residual_W_m3, &
         tolerance)
  end function mappings_match

  subroutine fill_three_body(out, value)
    type(fusion_three_body_cm_v1), intent(out) :: out
    real(c_double), intent(in) :: value

    out%kinetic_energy_J = value
    out%momentum_x_kg_m_s = value
    out%momentum_z_kg_m_s = value
    out%energy_residual_J = value
    out%momentum_residual_kg_m_s = value
  end subroutine fill_three_body

  logical function three_body_is_zero(out)
    type(fusion_three_body_cm_v1), intent(in) :: out

    three_body_is_zero = all(out%kinetic_energy_J == 0.0_c_double) .and. &
         all(out%momentum_x_kg_m_s == 0.0_c_double) .and. &
         all(out%momentum_z_kg_m_s == 0.0_c_double) .and. &
         out%energy_residual_J == 0.0_c_double .and. &
         out%momentum_residual_kg_m_s == 0.0_c_double
  end function three_body_is_zero

  subroutine verify_mapping_and_spill()
    real(c_double) :: edges(4), packet_energy(3), packet_rate(3)
    real(c_double) :: cell_birth(3)
    type(fusion_birth_mapping_v1) :: out
    integer(c_int) :: status

    edges = [1.0e-16_c_double, 3.0e-16_c_double, 5.0e-16_c_double, &
         9.0e-16_c_double]
    packet_energy = [3.0e-16_c_double, 1.0e-16_c_double, 1.0e-15_c_double]
    packet_rate = [4.0e20_c_double, 2.0e20_c_double, 3.0e20_c_double]
    cell_birth = -1.0_c_double
    call fill_mapping(out, -1.0_c_double)

    call fusion_map_birth_packets(3_c_int, edges, 3_c_int, packet_energy, &
         packet_rate, cell_birth, out, status)
    call check(status == PB11_STATUS_OK, 'birth mapping returns OK')
    call check(close_relative(cell_birth(1), 2.0e20_c_double, &
         1.0e-14_c_double) .and. &
         close_relative(cell_birth(2), 2.0e20_c_double, 1.0e-14_c_double) &
         .and. cell_birth(3) == 0.0_c_double, &
         'between-center packet is linearly interpolated')
    call check(close_relative(out%input_number_m3_s, 9.0e20_c_double, &
         1.0e-14_c_double) .and. &
         close_relative(out%input_energy_W_m3, 4.4e5_c_double, &
         1.0e-14_c_double), 'mapping ledger records packet input')
    call check(close_relative(out%mapped_number_m3_s, 4.0e20_c_double, &
         1.0e-14_c_double) .and. &
         close_relative(out%mapped_energy_W_m3, 1.2e5_c_double, &
         1.0e-14_c_double), 'mapping ledger records center interpolation')
    call check(close_relative(out%below_number_m3_s, 2.0e20_c_double, &
         1.0e-14_c_double) .and. &
         close_relative(out%below_energy_W_m3, 2.0e4_c_double, &
         1.0e-14_c_double), 'below-center spill is explicit')
    call check(close_relative(out%above_number_m3_s, 3.0e20_c_double, &
         1.0e-14_c_double) .and. &
         close_relative(out%above_energy_W_m3, 3.0e5_c_double, &
         1.0e-14_c_double), 'above-center spill is explicit')
    call check(abs(out%number_residual_m3_s) < 1.0e8_c_double .and. &
         abs(out%energy_residual_W_m3) < 1.0e-5_c_double, &
         'mapping number and energy residuals are small')
    call check(close_relative(out%input_number_m3_s, &
         out%mapped_number_m3_s + out%below_number_m3_s + &
         out%above_number_m3_s + out%number_residual_m3_s, 1.0e-12_c_double) &
         .and. close_relative(out%input_energy_W_m3, &
         out%mapped_energy_W_m3 + out%below_energy_W_m3 + &
         out%above_energy_W_m3 + out%energy_residual_W_m3, 1.0e-12_c_double), &
         'mapping ledger closes without dropping spill')
  end subroutine verify_mapping_and_spill

  subroutine verify_noncontiguous_mapping()
    real(c_double) :: edges_store(7), energy_store(5), rate_store(5)
    real(c_double) :: birth_store(5)
    real(c_double) :: edges_ref(4), energy_ref(3), rate_ref(3), birth_ref(3)
    type(fusion_birth_mapping_v1) :: out, out_ref
    integer(c_int) :: status, status_ref

    edges_store = [1.0e-16_c_double, -1.0_c_double, 3.0e-16_c_double, &
         -1.0_c_double, 5.0e-16_c_double, -1.0_c_double, 9.0e-16_c_double]
    energy_store = [3.0e-16_c_double, -1.0_c_double, 1.0e-16_c_double, &
         -1.0_c_double, 1.0e-15_c_double]
    rate_store = [4.0e20_c_double, -1.0_c_double, 2.0e20_c_double, &
         -1.0_c_double, 3.0e20_c_double]
    birth_store = -1.0_c_double
    edges_ref = [1.0e-16_c_double, 3.0e-16_c_double, 5.0e-16_c_double, &
         9.0e-16_c_double]
    energy_ref = [3.0e-16_c_double, 1.0e-16_c_double, 1.0e-15_c_double]
    rate_ref = [4.0e20_c_double, 2.0e20_c_double, 3.0e20_c_double]
    birth_ref = -1.0_c_double
    call fill_mapping(out, -1.0_c_double)
    call fill_mapping(out_ref, -1.0_c_double)

    call fusion_map_birth_packets(3_c_int, edges_store(1:7:2), 3_c_int, &
         energy_store(1:5:2), rate_store(1:5:2), birth_store(1:5:2), out, &
         status)
    call fusion_map_birth_packets(3_c_int, edges_ref, 3_c_int, energy_ref, &
         rate_ref, birth_ref, out_ref, status_ref)
    call check(status == PB11_STATUS_OK .and. status_ref == PB11_STATUS_OK, &
         'noncontiguous and contiguous mappings return OK')
    call check(close_relative(birth_store(1), birth_ref(1), 1.0e-14_c_double) &
         .and. close_relative(birth_store(3), birth_ref(2), 1.0e-14_c_double) &
         .and. close_relative(birth_store(5), birth_ref(3), 1.0e-14_c_double), &
         'noncontiguous birth output matches contiguous reference')
    call check(mappings_match(out, out_ref, 1.0e-14_c_double), &
         'noncontiguous packet inputs match contiguous ledger')
  end subroutine verify_noncontiguous_mapping

  subroutine verify_empty_packets()
    real(c_double) :: edges(4), empty_energy(0), empty_rate(0), birth(3)
    type(fusion_birth_mapping_v1) :: out
    integer(c_int) :: status

    edges = [1.0e-16_c_double, 3.0e-16_c_double, 5.0e-16_c_double, &
         9.0e-16_c_double]
    birth = -1.0_c_double
    call fill_mapping(out, -1.0_c_double)
    call fusion_map_birth_packets(3_c_int, edges, 0_c_int, empty_energy, &
         empty_rate, birth, out, status)
    call check(status == PB11_STATUS_OK, 'empty packet mapping returns OK')
    call check(all(birth == 0.0_c_double) .and. mapping_is_zero(out), &
         'empty packet arrays safely produce zero mapping')
  end subroutine verify_empty_packets

  subroutine verify_three_body_conservation()
    real(c_double), parameter :: mass = 1.0e-27_c_double
    real(c_double), parameter :: available = 1.0e-13_c_double
    real(c_double), parameter :: intermediate = 4.0e-14_c_double
    real(c_double), parameter :: cosine = 0.3_c_double
    type(fusion_three_body_cm_v1) :: out
    real(c_double) :: momentum_scale
    integer(c_int) :: status

    call fill_three_body(out, -1.0_c_double)
    call fusion_three_equal_sequential_cm(mass, available, intermediate, &
         cosine, out, status)
    call check(status == PB11_STATUS_OK, 'three-body kinematics returns OK')
    call check(all(ieee_is_finite(out%kinetic_energy_J)) .and. &
         all(out%kinetic_energy_J >= 0.0_c_double), &
         'three-body kinetic energies are finite and nonnegative')
    call check(close_relative(sum(out%kinetic_energy_J), available, &
         1.0e-10_c_double), 'three-body kinetic energy is conserved')
    momentum_scale = sum(abs(out%momentum_x_kg_m_s)) + &
         sum(abs(out%momentum_z_kg_m_s))
    call check(abs(sum(out%momentum_x_kg_m_s)) <= &
         1.0e-10_c_double * momentum_scale .and. &
         abs(sum(out%momentum_z_kg_m_s)) <= &
         1.0e-10_c_double * momentum_scale, &
         'three-body CM momentum is conserved')
    call check(abs(out%energy_residual_J) <= 1.0e-10_c_double * available .and. &
         abs(out%momentum_residual_kg_m_s) <= &
         1.0e-10_c_double * momentum_scale, &
         'three-body reported conservation residuals are small')
    call check(out%momentum_z_kg_m_s(1) > 0.0_c_double, &
         'primary three-body product follows the documented +z direction')
  end subroutine verify_three_body_conservation

  subroutine verify_errors()
    real(c_double) :: edges(3), packet_energy(1), packet_rate(1), birth(2)
    real(c_double) :: empty_energy(0), empty_rate(0)
    real(c_double) :: nan_rate
    type(fusion_birth_mapping_v1) :: mapping
    type(fusion_three_body_cm_v1) :: three_body
    integer(c_int) :: status

    edges = [1.0e-16_c_double, 3.0e-16_c_double, 5.0e-16_c_double]
    packet_energy = [2.0e-16_c_double]
    packet_rate = [1.0e20_c_double]
    birth = -1.0_c_double
    call fill_mapping(mapping, -1.0_c_double)
    call fusion_map_birth_packets(2_c_int, edges, 1_c_int, packet_energy, &
         packet_rate, birth, mapping, status)
    call check(status == PB11_STATUS_OK, 'valid two-cell mapping setup returns OK')

    nan_rate = ieee_value(0.0_c_double, ieee_quiet_nan)
    packet_rate(1) = nan_rate
    birth = -1.0_c_double
    call fusion_map_birth_packets(2_c_int, edges, 1_c_int, packet_energy, &
         packet_rate, birth, mapping, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'NaN packet rate returns invalid-argument status')
    call check(all(birth == 0.0_c_double) .and. mapping_is_zero(mapping), &
         'NaN packet rate clears mapping outputs')

    birth = -1.0_c_double
    call fusion_map_birth_packets(2_c_int, edges, 1_c_int, empty_energy, &
         empty_rate, birth, mapping, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         all(birth == 0.0_c_double) .and. mapping_is_zero(mapping), &
         'packet extent mismatch is rejected and cleared')

    call fill_three_body(three_body, -1.0_c_double)
    call fusion_three_equal_sequential_cm(-1.0e-27_c_double, 1.0e-13_c_double, &
         4.0e-14_c_double, 0.0_c_double, three_body, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         three_body_is_zero(three_body), &
         'negative product mass clears three-body output')

    call fill_three_body(three_body, -1.0_c_double)
    call fusion_three_equal_sequential_cm(1.0e-27_c_double, 1.0e-13_c_double, &
         4.0e-14_c_double, 1.2_c_double, three_body, status)
    call check(status == PB11_STATUS_OUT_OF_RANGE .and. &
         three_body_is_zero(three_body), 'invalid angle clears three-body output')
  end subroutine verify_errors

end program test_fusion_products_fortran
