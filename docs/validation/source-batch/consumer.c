#include <fusion_source_state.h>
int main(void){
 double edges[2]={0,2e-15},s[6]={0},t[6]={0},time;uint64_t tickets[2],epoch;
 fusion_source_state_v1 *states[2]={0,0};fusion_source_ledger_v1 ledger={0};
 for(int i=0;i<2;++i){
  if(fusion_c_source_state_create(1,edges,s,t,0,(uint64_t)i,&states[i]))return 1;
  if(fusion_c_source_state_begin(states[i],.1,&tickets[i]))return 2;
  if(fusion_c_source_state_stage(states[i],tickets[i],s,t,&ledger))return 3;
 }
 if(fusion_c_source_state_commit_many(2,states,tickets))return 4;
 for(int i=0;i<2;++i){
  if(fusion_c_source_state_snapshot(states[i],s,t,&ledger,&time,&epoch))return 5;
  if(time!=.1||epoch!=1)return 6;
  fusion_c_source_state_destroy(states[i]);
 }
 return 0;
}
