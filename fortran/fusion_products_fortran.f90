module fusion_products_fortran
  !! ISO_C_BINDING wrappers for product-birth mapping and exact sequential
  !! three-body CM kinematics.  No spectrum, branching, or angular model is
  !! selected here; all such inputs remain the caller's responsibility.
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

  type, bind(C), public :: fusion_birth_mapping_v1
     real(c_double) :: input_number_m3_s
     real(c_double) :: input_energy_W_m3
     real(c_double) :: mapped_number_m3_s
     real(c_double) :: mapped_energy_W_m3
     real(c_double) :: below_number_m3_s
     real(c_double) :: below_energy_W_m3
     real(c_double) :: above_number_m3_s
     real(c_double) :: above_energy_W_m3
     real(c_double) :: number_residual_m3_s
     real(c_double) :: energy_residual_W_m3
  end type fusion_birth_mapping_v1

  type, bind(C), public :: fusion_three_body_cm_v1
     real(c_double) :: kinetic_energy_J(3)
     real(c_double) :: momentum_x_kg_m_s(3)
     real(c_double) :: momentum_z_kg_m_s(3)
     real(c_double) :: energy_residual_J
     real(c_double) :: momentum_residual_kg_m_s
  end type fusion_three_body_cm_v1

  public :: fusion_map_birth_packets
  public :: fusion_three_equal_sequential_cm

  interface
     function c_fusion_map_birth_packets(cells, edges, packets, &
          packet_energy, packet_rate, cell_birth, out) &
          bind(C, name="fusion_c_map_birth_packets") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       integer(c_int), value :: packets
       type(c_ptr), value :: packet_energy
       type(c_ptr), value :: packet_rate
       type(c_ptr), value :: cell_birth
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_map_birth_packets

     function c_fusion_three_equal_sequential_cm(product_mass, available, &
          intermediate, cosine, out) &
          bind(C, name="fusion_c_three_equal_sequential_cm") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: product_mass
       real(c_double), value :: available
       real(c_double), value :: intermediate
       real(c_double), value :: cosine
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_three_equal_sequential_cm
  end interface

contains

  subroutine clear_mapping(out)
    type(fusion_birth_mapping_v1), intent(out) :: out

    out%input_number_m3_s = 0.0_c_double
    out%input_energy_W_m3 = 0.0_c_double
    out%mapped_number_m3_s = 0.0_c_double
    out%mapped_energy_W_m3 = 0.0_c_double
    out%below_number_m3_s = 0.0_c_double
    out%below_energy_W_m3 = 0.0_c_double
    out%above_number_m3_s = 0.0_c_double
    out%above_energy_W_m3 = 0.0_c_double
    out%number_residual_m3_s = 0.0_c_double
    out%energy_residual_W_m3 = 0.0_c_double
  end subroutine clear_mapping

  subroutine clear_three_body(out)
    type(fusion_three_body_cm_v1), intent(out) :: out

    out%kinetic_energy_J = 0.0_c_double
    out%momentum_x_kg_m_s = 0.0_c_double
    out%momentum_z_kg_m_s = 0.0_c_double
    out%energy_residual_J = 0.0_c_double
    out%momentum_residual_kg_m_s = 0.0_c_double
  end subroutine clear_three_body

  subroutine fusion_map_birth_packets(cells, edges_J, packets, &
       packet_energy_J, packet_rate_m3_s, cell_birth_m3_s, out, status)
    integer(c_int), intent(in) :: cells, packets
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: packet_energy_J(:)
    real(c_double), intent(in), target, contiguous :: packet_rate_m3_s(:)
    real(c_double), intent(out), target, contiguous :: cell_birth_m3_s(:)
    type(fusion_birth_mapping_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    integer(c_int) :: expected_edges
    type(c_ptr) :: edges_ptr, packet_energy_ptr, packet_rate_ptr
    type(c_ptr) :: cell_birth_ptr, out_ptr

    call clear_mapping(out)
    cell_birth_m3_s = 0.0_c_double
    status = PB11_STATUS_INVALID_ARGUMENT

    ! Validate dimensions before C_LOC.  Zero-extent packet arrays are valid
    ! when packets=0 and are passed as C_NULL_PTR to the C API.
    if (cells < 1_c_int .or. packets < 0_c_int) return
    if (cells >= huge(cells)) return
    expected_edges = cells + 1_c_int
    if (size(edges_J, kind=c_int) /= expected_edges) return
    if (size(packet_energy_J, kind=c_int) /= packets) return
    if (size(packet_rate_m3_s, kind=c_int) /= packets) return
    if (size(cell_birth_m3_s, kind=c_int) /= cells) return

    edges_ptr = c_loc(edges_J(1))
    cell_birth_ptr = c_loc(cell_birth_m3_s(1))
    out_ptr = c_loc(out)
    packet_energy_ptr = c_null_ptr
    packet_rate_ptr = c_null_ptr
    if (packets > 0_c_int) then
       packet_energy_ptr = c_loc(packet_energy_J(1))
       packet_rate_ptr = c_loc(packet_rate_m3_s(1))
    end if

    status = c_fusion_map_birth_packets(cells, edges_ptr, packets, &
         packet_energy_ptr, packet_rate_ptr, cell_birth_ptr, out_ptr)
    if (status /= PB11_STATUS_OK) then
       cell_birth_m3_s = 0.0_c_double
       call clear_mapping(out)
    end if
  end subroutine fusion_map_birth_packets

  subroutine fusion_three_equal_sequential_cm(product_mass_kg, &
       available_energy_J, intermediate_relative_energy_J, cos_theta, out, &
       status)
    real(c_double), intent(in) :: product_mass_kg, available_energy_J
    real(c_double), intent(in) :: intermediate_relative_energy_J, cos_theta
    type(fusion_three_body_cm_v1), intent(out), target :: out
    integer(c_int), intent(out) :: status

    call clear_three_body(out)
    status = c_fusion_three_equal_sequential_cm(product_mass_kg, &
         available_energy_J, intermediate_relative_energy_J, cos_theta, &
         c_loc(out))
    if (status /= PB11_STATUS_OK) call clear_three_body(out)
  end subroutine fusion_three_equal_sequential_cm

end module fusion_products_fortran
