program main
 use iso_c_binding
 use fusion_reaction_event_fortran
 use fusion_laboratory_fortran,only:fusion_particle_four_vector_v1,fusion_boost_ledger_v1
 implicit none
 real(c_double)::pa(3)=[1e-21_c_double,0._c_double,0._c_double],pb(3)=0._c_double
 real(c_double)::d(3)=[0._c_double,0._c_double,1._c_double]
 type(fusion_reaction_parent_v1)::p
 type(fusion_particle_four_vector_v1)::o(3)
 type(fusion_boost_ledger_v1)::l
 integer(c_int)::s
 call fusion_reaction_lab_event(0_c_int,0_c_int,pa,pb,d,91.84_c_double*1.602176634e-16_c_double, &
  .2_c_double,1._c_double,o,p,l,s)
 if(s/=0)stop 1
 if(abs(l%energy_residual_J)>1e-10_c_double*l%expected_kinetic_energy_J)stop 2
end program
