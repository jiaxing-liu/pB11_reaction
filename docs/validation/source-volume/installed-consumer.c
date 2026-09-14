#include <fusion_source_state.h>
#include <stdlib.h>
#include <stdio.h>
int main(void) {
 const double edges[2]={0,2},old[6]={4,0,0,0,0,0},zero[6]={0},trial[6]={8,0,0,0,0,0};
 double s[6],t[6],inert[6],ref=0,vol=0,time=0;
 fusion_source_state_v1* p=NULL;
 fusion_source_ledger_v1 source={0};fusion_transport_ledger_v1 move={0};uint64_t ticket=0,epoch=0;
 if(sizeof(move)!=336)return 1;
 if(fusion_c_source_state_create_volume(1,edges,old,zero,2,0,91,&p))return 2;
 if(fusion_c_source_state_begin(p,.5,&ticket))return 3;
 if(fusion_c_source_state_stage_volume(p,ticket,trial,zero,1,1,&source,zero,&move))return 4;
 if(fusion_c_source_state_commit(p,ticket))return 5;
 if(fusion_c_source_state_snapshot_volume(p,s,t,&source,inert,&move,&ref,&vol,&time,&epoch))return 6;
 if(s[0]!=8||ref!=2||vol!=1||time!=.5||epoch!=1)return 7;
 fusion_c_source_state_destroy(p);puts("PASS installed C11 volume-state consumer");return 0;
}
