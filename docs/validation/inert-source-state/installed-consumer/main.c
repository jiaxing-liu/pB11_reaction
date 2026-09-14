#include "fusion_source_state.h"
#include <stdlib.h>
#include <stdio.h>
int main(void){
 double edges[3]={0,2,4},s[12]={0,20},t[12]={0},sn[12]={5,15},heat[6]={10};
 fusion_source_state_v1 *state=NULL,*restored=NULL;fusion_source_ledger_v1 ledger={0};
 uint64_t ticket=0,epoch=0;double time=0,outheat[6]={0};size_t bytes=0,written=0;
 if(fusion_c_source_state_create(2,edges,s,t,0,987,&state))return 1;
 if(fusion_c_source_state_begin(state,.5,&ticket))return 2;
 if(fusion_c_source_state_stage_inert(state,ticket,sn,t,&ledger,heat))return 3;
 if(fusion_c_source_state_commit(state,ticket))return 4;
 if(fusion_c_source_state_pack_size(state,&bytes))return 5;
 unsigned char *buffer=malloc(bytes);if(!buffer)return 6;
 if(fusion_c_source_state_pack(state,buffer,bytes,&written)||written!=bytes)return 7;
 if(fusion_c_source_state_unpack(buffer,bytes,987,&restored))return 8;
 if(fusion_c_source_state_snapshot_inert(restored,s,t,&ledger,outheat,&time,&epoch))return 9;
 if(s[0]!=5||s[1]!=15||outheat[0]!=10||time!=.5||epoch!=1)return 10;
 if(fusion_c_source_state_snapshot(restored,s,t,&ledger,&time,&epoch)==0)return 11;
 fusion_c_source_state_destroy(state);fusion_c_source_state_destroy(restored);free(buffer);
 puts("Installed C11 inert source-state v2 restart passed");return 0;
}
