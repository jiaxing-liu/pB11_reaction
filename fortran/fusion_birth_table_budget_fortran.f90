module fusion_birth_table_budget_fortran
  use, intrinsic :: iso_c_binding, only : c_double,c_int,c_ptr,c_loc,c_null_ptr,c_size_t
  use fusion_birth_table_fortran, only : fusion_birth_table_control_v1,fusion_birth_table_info_v1, &
       fusion_birth_coefficients_v1,fusion_birth_table_create,fusion_birth_table_evaluate, &
       fusion_birth_table_destroy,fusion_birth_table_info,fusion_birth_table_pack_size, &
       fusion_birth_table_pack,fusion_birth_table_unpack,fusion_birth_table_content_digest, &
       fusion_birth_table_matches_request,fusion_birth_table_kernel_identity, &
       PB11_STATUS_OK,PB11_STATUS_INVALID_ARGUMENT,PB11_STATUS_OUT_OF_RANGE
  use fusion_thermal_birth_budget_fortran, only : fusion_thermal_birth_options_v1, &
       fusion_thermal_birth_budget_policy_v1,fusion_thermal_birth_budget_report_v1
  implicit none
  private
  public :: fusion_birth_table_control_v1,fusion_birth_table_info_v1,fusion_birth_coefficients_v1
  public :: fusion_thermal_birth_options_v1,fusion_thermal_birth_budget_policy_v1
  public :: fusion_thermal_birth_budget_report_v1
  public :: fusion_birth_table_create,fusion_birth_table_evaluate,fusion_birth_table_destroy
  public :: fusion_birth_table_info,fusion_birth_table_pack_size,fusion_birth_table_pack
  public :: fusion_birth_table_unpack,fusion_birth_table_content_digest
  public :: fusion_birth_table_matches_request,fusion_birth_table_kernel_identity
  public :: PB11_STATUS_OK,PB11_STATUS_INVALID_ARGUMENT,PB11_STATUS_OUT_OF_RANGE
  public :: pb_birth_table_create_budget,pb_birth_table_budget_info
  public :: pb_birth_table_budget_matches_request,pb_birth_table_budget_knot_report
  type,bind(C),public :: fusion_birth_table_budget_info_v1
    type(fusion_birth_table_info_v1) :: base
    type(fusion_thermal_birth_budget_policy_v1) :: policy
    real(c_double) :: max_sampled_source_number_bound,max_sampled_source_energy_bound
    real(c_double) :: max_sampled_number_roundoff_allowance,max_sampled_energy_roundoff_allowance
    real(c_double) :: max_validated_interpolation_number_L1,max_validated_interpolation_energy_L1
  end type
  interface
    function c_create(lo,hi,source,policy,control,cells,edges,out) &
         bind(C,name="fusion_c_pb_birth_table_create_budget") result(status)
      import :: c_double,c_int,c_ptr
      real(c_double),value :: lo,hi
      type(c_ptr),value :: source,policy,control
      integer(c_int),value :: cells
      type(c_ptr),value :: edges,out
      integer(c_int) :: status
    end function
    function c_info(table,out) bind(C,name="fusion_c_pb_birth_table_budget_info") result(status)
      import :: c_ptr,c_int
      type(c_ptr),value :: table,out
      integer(c_int) :: status
    end function
    function c_matches(table,lo,hi,source,policy,control,cells,edges,answer) &
         bind(C,name="fusion_c_pb_birth_table_budget_matches_request") result(status)
      import :: c_double,c_int,c_ptr
      type(c_ptr),value :: table
      real(c_double),value :: lo,hi
      type(c_ptr),value :: source,policy,control
      integer(c_int),value :: cells
      type(c_ptr),value :: edges,answer
      integer(c_int) :: status
    end function
    function c_report(table,index,temperature,out) &
         bind(C,name="fusion_c_pb_birth_table_budget_knot_report") result(status)
      import :: c_int,c_ptr
      type(c_ptr),value :: table
      integer(c_int),value :: index
      type(c_ptr),value :: temperature,out
      integer(c_int) :: status
    end function
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
  subroutine clear_budget_info(info)
    type(fusion_birth_table_budget_info_v1),intent(out)::info
    call clear_birth_table_info(info%base)
    info%policy=fusion_thermal_birth_budget_policy_v1(0_c_int,0_c_int,0_c_int,0_c_int,0._c_double,0._c_double)
    info%max_sampled_source_number_bound=0;info%max_sampled_source_energy_bound=0
    info%max_sampled_number_roundoff_allowance=0;info%max_sampled_energy_roundoff_allowance=0
    info%max_validated_interpolation_number_L1=0;info%max_validated_interpolation_energy_L1=0
  end subroutine
  subroutine pb_birth_table_create_budget(lo,hi,source,policy,control,edges_J,table,status)
    real(c_double),intent(in)::lo,hi
    type(fusion_thermal_birth_options_v1),intent(in),target::source
    type(fusion_thermal_birth_budget_policy_v1),intent(in),target::policy
    type(fusion_birth_table_control_v1),intent(in),target::control
    real(c_double),intent(in),target,contiguous::edges_J(:)
    type(c_ptr),intent(out),target::table
    integer(c_int),intent(out)::status
    type(c_ptr),target::candidate
    integer(c_size_t)::count
    table=c_null_ptr;candidate=c_null_ptr;status=PB11_STATUS_INVALID_ARGUMENT
    count=size(edges_J,kind=c_size_t)
    if(count<2_c_size_t.or.count>100001_c_size_t)return
    status=c_create(lo,hi,c_loc(source),c_loc(policy),c_loc(control),int(count-1,c_int), &
         c_loc(edges_J(1)),c_loc(candidate))
    if(status==PB11_STATUS_OK)table=candidate
  end subroutine
  subroutine pb_birth_table_budget_info(table,info,status)
    type(c_ptr),intent(in)::table
    type(fusion_birth_table_budget_info_v1),intent(out),target::info
    integer(c_int),intent(out)::status
    call clear_budget_info(info)
    status=c_info(table,c_loc(info))
    if(status/=PB11_STATUS_OK)call clear_budget_info(info)
  end subroutine
  subroutine pb_birth_table_budget_matches_request(table,lo,hi,source,policy,control,edges_J,matches,status)
    type(c_ptr),intent(in)::table
    real(c_double),intent(in)::lo,hi
    type(fusion_thermal_birth_options_v1),intent(in),target::source
    type(fusion_thermal_birth_budget_policy_v1),intent(in),target::policy
    type(fusion_birth_table_control_v1),intent(in),target::control
    real(c_double),intent(in),target,contiguous::edges_J(:)
    logical,intent(out)::matches
    integer(c_int),intent(out)::status
    integer(c_int),target::answer
    integer(c_size_t)::count
    matches=.false.;answer=0;status=PB11_STATUS_INVALID_ARGUMENT
    count=size(edges_J,kind=c_size_t)
    if(count<2_c_size_t.or.count>100001_c_size_t)return
    status=c_matches(table,lo,hi,c_loc(source),c_loc(policy),c_loc(control),int(count-1,c_int), &
         c_loc(edges_J(1)),c_loc(answer))
    if(status==PB11_STATUS_OK)matches=answer==1_c_int
  end subroutine
  subroutine pb_birth_table_budget_knot_report(table,index,kT_J,report,status)
    type(c_ptr),intent(in)::table
    integer(c_int),intent(in)::index
    real(c_double),intent(out),target::kT_J
    type(fusion_thermal_birth_budget_report_v1),intent(out),target::report
    integer(c_int),intent(out)::status
    kT_J=0;report=fusion_thermal_birth_budget_report_v1()
    status=c_report(table,index,c_loc(kT_J),c_loc(report))
    if(status/=PB11_STATUS_OK)then
      kT_J=0;report=fusion_thermal_birth_budget_report_v1()
    endif
  end subroutine
end module
