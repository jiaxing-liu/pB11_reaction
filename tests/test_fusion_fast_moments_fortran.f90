program test_moments
 use iso_c_binding, only:c_double,c_int,c_sizeof
 use fusion_fast_moments_fortran
 implicit none
 real(c_double) :: edges(2),s(6),t(6)
 type(fusion_fast_moments_v1) :: out
 integer(c_int) :: status
 edges=[0._c_double,2.e-15_c_double];s=0;t=0;s(5)=1;t(5)=2
 call fusion_fast_moments(edges,s,t,out,status)
 if(status/=0.or.c_sizeof(out)/=120) stop 1
 if(out%number_m3(5)/=3.or.out%charge_number_m3/=6.or.out%charge_squared_number_m3/=12) stop 2
 if(abs(out%pressure_Pa-2.e-15_c_double)>1.e-28_c_double) stop 3
 call fusion_fast_moments(edges,s(1:5),t,out,status)
 if(status==0.or.any(out%number_m3/=0).or.any(out%energy_J_m3/=0)) stop 4
 if(out%charge_number_m3/=0.or.out%charge_squared_number_m3/=0.or.out%pressure_Pa/=0) stop 5
 print *, 'PASS fast moments Fortran ABI, S/T alpha charge and pressure, bad extents'
end program
