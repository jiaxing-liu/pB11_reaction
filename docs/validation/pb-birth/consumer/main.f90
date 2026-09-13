program main
 use iso_c_binding
 use fusion_pb_birth_fortran
 implicit none
 real(c_double),parameter::u=1.602176634e-13_c_double
 real(c_double)::e(3)=[0._c_double,4*u,8*u],b(2)
 type(fusion_pb_birth_v1)::r
 integer(c_int)::s
 call fusion_pb_cm_source_grid(8.84_c_double*u,.09184_c_double*u,.05_c_double,.95_c_double, &
  13_c_int,1_c_int,.001_c_double*u,.76_c_double,0._c_double,32_c_int,32_c_int,2_c_int,e,b,r,s)
 if(s/=0)stop 1
 if(abs(r%mapped_number+r%below_number+r%above_number-3)>1e-12_c_double)stop 2
 if(abs((r%mapped_energy_J+r%below_energy_J+r%above_energy_J)/(8.84_c_double*u)-1)>1e-12_c_double)stop 3
end program
