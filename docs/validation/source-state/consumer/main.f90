program consumer
 use iso_c_binding
 use fusion_source_state_fortran
 implicit none
 type(c_ptr) :: state
 real(c_double) :: edges(2),s(6),t(6),time
 type(fusion_source_ledger_v1) :: ledger
 integer(c_int) :: rc
 integer(c_int64_t) :: ticket,epoch
 edges=[0.0_c_double,1.0_c_double];s=0;t=0
 call fusion_source_state_create(1_c_int,edges,s,t,0.0_c_double,123_c_int64_t,state,rc)
 if(rc/=0)stop 1
 call fusion_source_state_snapshot(state,s,t,ledger,time,epoch,rc)
 if(rc/=0)stop 2
 call fusion_source_state_begin(state,.1_c_double,ticket,rc)
 if(rc/=0)stop 3
 call fusion_source_state_stage(state,ticket,s,t,ledger,rc)
 if(rc/=0)stop 4
 call fusion_source_state_commit(state,ticket,rc)
 if(rc/=0)stop 5
 call fusion_source_state_destroy(state)
 if(c_associated(state))stop 6
end program
