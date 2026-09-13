program main
 use iso_c_binding
 use fusion_pb_population_fortran
 use fusion_nuclear_data_fortran,only:fusion_nuclear_mass_v1,fusion_nuclear_mass
 implicit none
 type(fusion_nuclear_mass_v1)::p,b
 type(fusion_pb_population_rates_v1)::r
 integer(c_int)::s
 call fusion_nuclear_mass(0_c_int,p,s)
 if(s/=0)stop 1
 call fusion_nuclear_mass(5_c_int,b,s)
 if(s/=0)stop 2
 call fusion_pb_thermal_population_rates(1_c_int,0_c_int,p%mass_kg,b%mass_kg, &
  1.602176634e-15_c_double,1.1215236438e-15_c_double,.051_c_double,1._c_double,r,s)
 if(s/=0)stop 3
 if(abs((r%alpha0_peak%resolved_reactivity_m3_s+r%narrow_remainder%resolved_reactivity_m3_s+ &
  r%other_remainder%resolved_reactivity_m3_s)/r%total%resolved_reactivity_m3_s-1)>1e-8_c_double)stop 4
end program
