program test_fusion_birth_table_budget_fortran
  use, intrinsic :: iso_c_binding, only:c_double,c_int,c_size_t,c_int8_t,c_sizeof,c_ptr,c_null_ptr,c_associated
  use fusion_birth_table_budget_fortran
  implicit none
  integer,parameter::n=8
  real(c_double),parameter::MeV=1.602176634e-13_c_double,keV=1.602176634e-16_c_double
  type(fusion_thermal_birth_options_v1)::source
  type(fusion_thermal_birth_budget_policy_v1)::policy,wrong
  type(fusion_birth_table_control_v1)::control
  type(fusion_birth_table_budget_info_v1)::info
  type(fusion_birth_table_info_v1)::base
  type(fusion_birth_coefficients_v1)::coeff,copycoeff
  type(fusion_thermal_birth_budget_report_v1)::report
  type(c_ptr)::table,copy,bad
  real(c_double)::edges(n+1),empty(0),birth(n,7),copied(n,7),lo,hi,T
  integer(c_int8_t),allocatable::bytes(:)
  integer(c_int8_t)::digest(32),copydigest(32)
  integer(c_size_t)::required,written,length
  integer(c_int)::status
  integer::i,failures
  logical::matches
  failures=0;table=c_null_ptr;copy=c_null_ptr
  call check(c_sizeof(info)==c_sizeof(base)+c_sizeof(policy)+6*c_sizeof(0._c_double),'nested budget info size')
  call check(c_sizeof(info)==328_c_size_t,'expected C info layout')
  source=fusion_thermal_birth_options_v1(2*MeV,40._c_double,.09184_c_double*MeV,.001_c_double*MeV, &
       .76_c_double,0._c_double,.051_c_double,1._c_double,1,0,0,13,0,4,4,4,4)
  policy=fusion_thermal_birth_budget_policy_v1(1,4,4,0,1.e-4_c_double,1.e-4_c_double)
  ! A deliberately loose tiny-domain numerical contract fixture, not physics qualification.
  control=fusion_birth_table_control_v1(.01_c_double,.01_c_double,.01_c_double,.01_c_double, &
       .05_c_double,.05_c_double,64,512,10)
  lo=.5_c_double*keV;hi=lo*1.000001_c_double
  do i=1,n+1
    edges(i)=20*MeV*(i-1)/n
  enddo
  call pb_birth_table_create_budget(lo,hi,source,policy,control,edges,table,status)
  call check(status==0.and.c_associated(table),'create explicit budget table')
  if(status/=0)then
    write(*,*)'create failed status ',status
    error stop 1
  endif
  call pb_birth_table_budget_info(table,info,status)
  call check(status==0.and.info%base%cells==n.and.info%base%channel==0,'nested base fields')
  call check(info%policy%version==1.and.info%policy%coarse_nq==4,'nested policy fields')
  call check(info%policy%energy_L1_tolerance==policy%energy_L1_tolerance,'policy double layout')
  call check(info%max_sampled_source_number_bound>=info%max_sampled_number_roundoff_allowance, &
       'effective source bounds include allowance')
  call check(info%base%max_validated_number_L1>=info%max_validated_interpolation_number_L1, &
       'combined interpolation/source accounting')
  call fusion_birth_table_info(table,base,status)
  call check(status==0.and.base%knots==info%base%knots,'generic info remains usable')
  call pb_birth_table_budget_matches_request(table,lo,hi,source,policy,control,edges,matches,status)
  call check(status==0.and.matches,'explicit request match')
  call fusion_birth_table_matches_request(table,0_c_int,lo,hi,source,control,edges,matches,status)
  call check(status==0.and..not.matches,'old matcher rejects budget table')
  wrong=policy;wrong%number_L1_tolerance=.0002_c_double
  call pb_birth_table_budget_matches_request(table,lo,hi,source,wrong,control,edges,matches,status)
  call check(status==0.and..not.matches,'different numerical policy rejects match')
  call pb_birth_table_budget_knot_report(table,0_c_int,T,report,status)
  call check(status==0.and.T==lo.and.report%version==1,'zero-based first knot report')
  call check(report%full_nq==4.and.report%coarse_ncos==4,'knot report layout')
  call fusion_birth_table_evaluate(table,lo,birth,coeff,status)
  call check(status==0,'generic evaluate budget table')
  call fusion_birth_table_pack_size(table,required,status)
  call check(status==0.and.required>0,'generic pack size')
  allocate(bytes(required))
  call fusion_birth_table_pack(table,bytes,required,written,status)
  call check(status==0.and.written==required,'generic format2 pack')
  call fusion_birth_table_unpack(bytes,written,copy,status)
  call check(status==0.and.c_associated(copy),'generic format2 unpack')
  if(c_associated(copy))then
    call pb_birth_table_budget_matches_request(copy,lo,hi,source,policy,control,edges,matches,status)
    call check(status==0.and.matches,'policy persists after unpack')
    call fusion_birth_table_evaluate(copy,lo,copied,copycoeff,status)
    call check(status==0.and.all(copied==birth),'unpack evaluate identical')
    call fusion_birth_table_content_digest(table,digest,length,status)
    call check(status==0.and.length==required,'generic digest packed-size metadata')
    call fusion_birth_table_content_digest(copy,copydigest,length,status)
    call check(status==0.and.all(copydigest==digest),'unpack digest identical')
    call fusion_birth_table_destroy(copy)
  endif
  deallocate(bytes)
  matches=.true.
  call pb_birth_table_budget_matches_request(table,lo,hi,source,policy,control,empty,matches,status)
  call check(status/=0.and..not.matches,'bad shape clears matches')
  bad=table
  call pb_birth_table_create_budget(lo,hi,source,policy,control,empty,bad,status)
  call check(status/=0.and..not.c_associated(bad),'bad shape clears output handle')
  report%version=9;report%relative_nodes=9;report%energy_L1_bound=9;T=9
  call pb_birth_table_budget_knot_report(table,-1_c_int,T,report,status)
  call check(status/=0.and.T==0.and.report%version==0.and.report%relative_nodes==0.and. &
       report%energy_L1_bound==0,'invalid index clears temperature/report')
  info%base%cells=9;info%policy%version=9;info%max_sampled_source_energy_bound=9
  call pb_birth_table_budget_info(c_null_ptr,info,status)
  call check(status/=0.and.info%base%cells==0.and.info%policy%version==0.and. &
       info%max_sampled_source_energy_bound==0,'null table clears nested info')
  call fusion_birth_table_destroy(table)
  ! Explicit budget creation is still a distinct contract when budgets are zero.
  policy%number_L1_tolerance=0;policy%energy_L1_tolerance=0
  call pb_birth_table_create_budget(lo,hi,source,policy,control,edges,table,status)
  call check(status==0.and.c_associated(table),'zero-budget explicit table')
  if(c_associated(table))then
    call fusion_birth_table_matches_request(table,0_c_int,lo,hi,source,control,edges,matches,status)
    call check(status==0.and..not.matches,'old matcher rejects explicit zero budget')
    call fusion_birth_table_destroy(table)
  endif
  if(failures/=0)error stop 1
  print *,'PASS Fortran budget table layout, lifecycle, format2, matching and failure clearing'
contains
  subroutine check(ok,message)
    logical,intent(in)::ok
    character(*),intent(in)::message
    if(.not.ok)then
      failures=failures+1;print *,'FAIL ',message
    endif
  end subroutine
end program
