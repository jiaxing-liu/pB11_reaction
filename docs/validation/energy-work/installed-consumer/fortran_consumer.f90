program test_energy_work
 use iso_c_binding,only:c_int,c_double,c_sizeof,c_size_t
 use fusion_energy_work_fortran
 implicit none
 real(c_double)::edges(2),old(1),trial(1),face(2)
 type(fusion_energy_work_ledger_v1)::ledger
 integer(c_int)::status
 edges=[1._c_double,3._c_double];old=4._c_double
 if(c_sizeof(ledger)/=88_c_size_t)stop 1
 call fusion_energy_work_trial(1_c_int,.5_c_double,1._c_double,edges,old,trial,face,ledger,status)
 if(status/=0.or.abs(trial(1)-3.2_c_double)>1.e-13_c_double)stop 2
 if(abs(ledger%work_on_particles_J_m3+.8_c_double)>1.e-13_c_double)stop 3
 if(abs(ledger%lower_energy_J_m3-.8_c_double)>1.e-13_c_double)stop 4
 trial=9;face=9
 call fusion_energy_work_trial(1_c_int,.5_c_double,1._c_double,edges(:1),old,trial,face,ledger,status)
 if(status==0.or.any(trial/=0).or.any(face/=0).or.ledger%initial_number_m3/=0)stop 5
 print *, 'PASS Fortran energy-work ABI, signed cooling and extent rejection'
end program
