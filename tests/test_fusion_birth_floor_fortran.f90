program test_birth_floor_fortran
 use iso_c_binding
 use fusion_birth_floor_fortran
 implicit none
 type(fusion_birth_floor_options_v1)::o
 type(fusion_birth_floor_ledger_v1)::l
 real(c_double)::n(6),u(6)
 integer(c_int)::status
 o=fusion_birth_floor_options_v1(.25_c_double,1._c_double,8._c_double,.5_c_double,.5_c_double)
 n=0;u=0;n(5)=4;u(5)=.5_c_double
 status=fusion_birth_floor_project(o,n,u,l)
 if(status/=0)error stop 'binding rejected valid input'
 if(any(l%mapped_number_m3/=n))error stop 'binding particle count changed'
 if(l%mapped_energy_J_m3(5)/=1._c_double)error stop 'binding mapped energy'
 if(l%ion_energy_correction_J_m3(5)/=-.5_c_double)error stop 'binding correction'
 if(l%remaining_ion_energy_J_m3/=7.5_c_double)error stop 'binding remaining bath'
 o%max_center_over_ion_kT=.125_c_double
 status=fusion_birth_floor_project(o,n,u,l)
 if(status/=3.or.any(l%mapped_number_m3/=0).or.l%remaining_ion_energy_J_m3/=0) &
   error stop 'binding invalid output not cleared'
 print *, 'PASS birth floor Fortran ABI, energy correction and failure atomicity'
end program
