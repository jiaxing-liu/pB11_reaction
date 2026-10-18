module fusion_thermal_birth_budget_fortran
  !! Isolated explicit finite-source budget binding; not host/cache policy.
  use, intrinsic :: iso_c_binding, only : c_double,c_int,c_int64_t,c_size_t,c_ptr,c_loc
  use fusion_thermal_birth_fortran, only : fusion_thermal_birth_options_v1, &
       fusion_thermal_birth_v1,fusion_thermal_birth_grid, &
       PB11_STATUS_OK,PB11_STATUS_INVALID_ARGUMENT,PB11_STATUS_OUT_OF_RANGE
  implicit none
  private
  public :: fusion_thermal_birth_options_v1,fusion_thermal_birth_v1
  public :: fusion_thermal_birth_grid
  public :: PB11_STATUS_OK,PB11_STATUS_INVALID_ARGUMENT,PB11_STATUS_OUT_OF_RANGE
  public :: pb_thermal_birth_budget_grid

  type,bind(C),public :: fusion_thermal_birth_budget_policy_v1
    integer(c_int) :: version,coarse_nq,coarse_ncos,reserved
    real(c_double) :: number_L1_tolerance,energy_L1_tolerance
  end type

  type,bind(C),public :: fusion_thermal_birth_budget_report_v1
    integer(c_int) :: version=0,full_nq=0,full_ncos=0,coarse_nq=0,coarse_ncos=0,mixed=0
    ! C uint64_t uses the interoperable signed64 storage kind. Reported counts
    ! must fit the positive signed range when interpreted as Fortran integers.
    integer(c_int64_t) :: relative_nodes=0,cm_nodes=0,joint_nodes=0,high_pairs=0,low_pairs=0
    integer(c_int64_t) :: high_relative_nodes=0,low_relative_nodes=0
    real(c_double) :: number_L1_bound=0,energy_L1_bound=0
    real(c_double) :: number_roundoff_allowance=0,energy_roundoff_allowance=0
    real(c_double) :: reaction_weight_m3_s=0,product_energy_weight_J_m3_s=0
  end type

  interface
    function c_pb_budget(kT_J,options,policy,cells,edges,birth,out,report) &
         bind(C,name="fusion_c_pb_thermal_birth_budget_grid") result(status)
      import :: c_double,c_int,c_ptr
      real(c_double),value :: kT_J
      type(c_ptr),value :: options,policy
      integer(c_int),value :: cells
      type(c_ptr),value :: edges,birth,out,report
      integer(c_int) :: status
    end function
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

  subroutine pb_thermal_birth_budget_grid(kT_J,options,policy,edges_J,birth,out,report,status)
    real(c_double),intent(in) :: kT_J
    type(fusion_thermal_birth_options_v1),intent(in),target :: options
    type(fusion_thermal_birth_budget_policy_v1),intent(in),target :: policy
    real(c_double),intent(in),target,contiguous :: edges_J(:)
    real(c_double),intent(out),target,contiguous :: birth(:,:)
    type(fusion_thermal_birth_v1),intent(out),target :: out
    type(fusion_thermal_birth_budget_report_v1),intent(out),target :: report
    integer(c_int),intent(out) :: status
    integer(c_size_t) :: extent
    integer(c_int) :: cells

    birth=0.0_c_double
    call clear_thermal_birth(out)
    report=fusion_thermal_birth_budget_report_v1()
    status=PB11_STATUS_INVALID_ARGUMENT
    ! Inspect all extents with a wide size kind before narrowing or C_LOC.
    extent=size(birth,1,kind=c_size_t)
    if(extent<1_c_size_t.or.extent>100000_c_size_t)return
    if(size(birth,2,kind=c_size_t)/=7_c_size_t)return
    if(size(edges_J,kind=c_size_t)/=extent+1_c_size_t)return
    cells=int(extent,c_int)
    status=c_pb_budget(kT_J,c_loc(options),c_loc(policy),cells,c_loc(edges_J(1)), &
         c_loc(birth(1,1)),c_loc(out),c_loc(report))
    if(status/=PB11_STATUS_OK)then
      birth=0.0_c_double
      call clear_thermal_birth(out)
      report=fusion_thermal_birth_budget_report_v1()
    endif
  end subroutine
end module
