program test_fusion_thermal_birth_budget_fortran
  use, intrinsic :: iso_c_binding, only : c_double,c_int,c_int64_t,c_size_t,c_sizeof,int8=>c_int8_t
  use fusion_thermal_birth_budget_fortran
  implicit none
  integer,parameter :: n=8
  real(c_double),parameter :: MeV=1.602176634e-13_c_double,keV=1.602176634e-16_c_double
  type(fusion_thermal_birth_options_v1) :: options
  type(fusion_thermal_birth_budget_policy_v1) :: policy
  type(fusion_thermal_birth_budget_report_v1) :: report
  type(fusion_thermal_birth_v1) :: out,full
  real(c_double) :: edges(n+1),birth(n,7),reference(n,7),bad_species(n,6),bad_edges(n)
  real(c_double) :: empty(0,7),empty_edges(0)
  real(c_double),allocatable :: too_many(:,:),large_edges(:)
  integer(c_int) :: status
  integer :: i,failures
  failures=0
  call check(c_sizeof(policy)==4*c_sizeof(0_c_int)+2*c_sizeof(0._c_double),'policy ABI size')
  call check(c_sizeof(report)==6*c_sizeof(0_c_int)+7*c_sizeof(0_c_int64_t)+ &
       6*c_sizeof(0._c_double),'report ABI size')
  call check(c_sizeof(options)==104_c_size_t.and.c_sizeof(out)==352_c_size_t,'reused original ABI types')
  options=fusion_thermal_birth_options_v1(2*MeV,40._c_double,.09184_c_double*MeV, &
       .001_c_double*MeV,.76_c_double,0._c_double,.051_c_double,1._c_double, &
       1_c_int,0_c_int,0_c_int,13_c_int,0_c_int,4_c_int,4_c_int,4_c_int,4_c_int)
  policy=fusion_thermal_birth_budget_policy_v1(1_c_int,4_c_int,4_c_int,0_c_int,0._c_double,0._c_double)
  do i=1,n+1
    edges(i)=20*MeV*(i-1)/n
  enddo
  call fusion_thermal_birth_grid(0_c_int,.5_c_double*keV,options,edges,reference,full,status)
  call check(status==PB11_STATUS_OK,'original pB full source call')
  call pb_thermal_birth_budget_grid(.5_c_double*keV,options,policy,edges,birth,out,report,status)
  call check(status==PB11_STATUS_OK,'zero budget call')
  call check(all(transfer(birth,[0_int8],size=56*8)==transfer(reference,[0_int8],size=56*8)), &
       'zero budget birth bitwise equals full')
  call check(all(transfer(out,[0_int8],size=352)==transfer(full,[0_int8],size=352)), &
       'zero budget result bitwise equals full')
  call check(report%version==1.and.report%full_nq==4.and.report%full_ncos==4,'report int fields')
  call check(report%coarse_nq==4.and.report%coarse_ncos==4.and.report%mixed==0,'coarse/report layout')
  call check(report%number_L1_bound==0.and.report%energy_L1_bound==0,'zero bounds')
  call check(report%relative_nodes>0_c_int64_t.and.report%cm_nodes>0_c_int64_t,'uint64 count layout')
  call check(report%joint_nodes==report%high_pairs.and.report%low_pairs==0_c_int64_t,'all parents full')
  call check(report%reaction_weight_m3_s==out%reactivity_m3_s,'last report double layout')

  bad_edges=edges(1:n)
  call dirty();birth=-7
  call pb_thermal_birth_budget_grid(.5_c_double*keV,options,policy,bad_edges,birth,out,report,status)
  call rejected('bad edge extent');call check(all(birth==0),'bad edges clears birth')
  call dirty();bad_species=-7
  call pb_thermal_birth_budget_grid(.5_c_double*keV,options,policy,edges,bad_species,out,report,status)
  call rejected('bad species extent');call check(all(bad_species==0),'bad species clears all birth')
  call dirty()
  call pb_thermal_birth_budget_grid(.5_c_double*keV,options,policy,edges,empty,out,report,status)
  call rejected('zero cells safe')
  call dirty();birth=-7
  call pb_thermal_birth_budget_grid(.5_c_double*keV,options,policy,empty_edges,birth,out,report,status)
  call rejected('zero edges safe');call check(all(birth==0),'zero edges clears birth')
  allocate(too_many(100001,7),large_edges(100002));large_edges=0;too_many=-7
  call dirty()
  call pb_thermal_birth_budget_grid(.5_c_double*keV,options,policy,large_edges,too_many,out,report,status)
  call rejected('cell upper bound');call check(all(too_many==0),'oversized birth cleared')
  deallocate(too_many,large_edges)
  call dirty();birth=-7;policy%reserved=1
  call pb_thermal_birth_budget_grid(.5_c_double*keV,options,policy,edges,birth,out,report,status)
  call rejected('C policy failure');call check(all(birth==0),'C failure clears birth')
  if(failures/=0)error stop 1
  write(*,'(A)')'PASS Fortran budget layout, original full/zero budget and atomic shape rejection'
contains
  subroutine check(ok,message)
    logical,intent(in)::ok
    character(*),intent(in)::message
    if(.not.ok)then
      failures=failures+1;write(*,'(A)')'FAIL '//message
    endif
  end subroutine
  subroutine dirty()
    out=transfer(spread(-7._c_double,1,44),out)
    report%version=-7;report%full_nq=-7;report%full_ncos=-7
    report%coarse_nq=-7;report%coarse_ncos=-7;report%mixed=-7
    report%relative_nodes=-7;report%cm_nodes=-7;report%joint_nodes=-7
    report%high_pairs=-7;report%low_pairs=-7
    report%high_relative_nodes=-7;report%low_relative_nodes=-7
    report%number_L1_bound=-7;report%energy_L1_bound=-7
    report%number_roundoff_allowance=-7;report%energy_roundoff_allowance=-7
    report%reaction_weight_m3_s=-7;report%product_energy_weight_J_m3_s=-7
  end subroutine
  subroutine rejected(message)
    character(*),intent(in)::message
    call check(status/=PB11_STATUS_OK,message//' rejected')
    call check(all(transfer(out,[0._c_double],size=44)==0),message//' out cleared')
    call check(report%version==0.and.report%full_nq==0.and.report%full_ncos==0.and. &
         report%coarse_nq==0.and.report%coarse_ncos==0.and.report%mixed==0,message//' report ints clear')
    call check(report%relative_nodes==0.and.report%cm_nodes==0.and.report%joint_nodes==0.and. &
         report%high_pairs==0.and.report%low_pairs==0.and.report%high_relative_nodes==0.and. &
         report%low_relative_nodes==0,message//' report counts clear')
    call check(report%number_L1_bound==0.and.report%energy_L1_bound==0.and. &
         report%number_roundoff_allowance==0.and.report%energy_roundoff_allowance==0.and. &
         report%reaction_weight_m3_s==0.and.report%product_energy_weight_J_m3_s==0, &
         message//' report doubles clear')
  end subroutine
end program
