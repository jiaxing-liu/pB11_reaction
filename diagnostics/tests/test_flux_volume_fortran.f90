program test_flux_volume_fortran
 use,intrinsic::iso_c_binding,only:c_int,c_double
 use fusion_flux_volume_fortran
 implicit none
 real(c_double)::x(3)=[0._c_double,.5_c_double,1._c_double]
 real(c_double)::a(9),d(9),f(3)=1._c_double,c(3)=[0._c_double,.5_c_double,1._c_double]
 real(c_double)::expected,pi
 type(fusion_flux_volume_diagnostic_v1)::out
 integer(c_int)::status
 integer::j
 do j=1,3
  a(3*j-2:3*j)=[.7_c_double,.4_c_double*x(j),.3_c_double*x(j)]
  d(3*j-2:3*j)=[0._c_double,.4_c_double,.3_c_double]
 enddo
 pi=acos(-1._c_double)
 status=fusion_c_flux_core_enclosed_volume_diagnostic_v1(3_c_int,1_c_int,x,a,d,f,c, &
                                                      1_c_int,.75_c_double,out)
 expected=2*pi*pi*.7_c_double*.4_c_double*.3_c_double*.75_c_double**2
 if(status/=0.or.out%version/=1.or.out%angular_nodes/=4.or.out%boundary_orientation/=1)stop 1
 if(abs(out%signed_volume_m3/expected-1)>1.e-12_c_double)stop 2
 status=fusion_c_flux_core_enclosed_volume_diagnostic_v1(3_c_int,1_c_int,x,a,d,f,c, &
                                                      1_c_int,-1._c_double,out)
 if(status==0.or.out%version/=0.or.out%absolute_volume_m3/=0)stop 3
 print *,'Fortran C ABI volume and error-clearing checks pass'
end program
