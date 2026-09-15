program test_source_grid
 use iso_c_binding
 use fusion_source_state_fortran
 implicit none
 type(c_ptr)::state
 real(c_double)::edges(3)=[0._c_double,1._c_double,2._c_double],s(12)=0,t(12)=0
 integer(c_int)::rc
 logical::matches
 call fusion_source_state_create(2_c_int,edges,s,t,0._c_double,42_c_int64_t,state,rc)
 if(rc/=0)stop 1
 call fusion_source_state_grid_matches(state,edges,matches,rc)
 if(rc/=0.or..not.matches)stop 2
 edges(2)=.5_c_double
 call fusion_source_state_grid_matches(state,edges,matches,rc)
 if(rc/=0.or.matches)stop 3
 call fusion_source_state_grid_matches(state,edges(:1),matches,rc)
 if(rc==0.or.matches)stop 4
 call fusion_source_state_destroy(state)
 call fusion_source_state_grid_matches(state,edges,matches,rc)
 if(rc==0.or.matches)stop 5
 print *, 'PASS source grid Fortran matching and guards'
end program
