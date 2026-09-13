#include <fusion_source_state.h>
#include <stdlib.h>
int main(void){double edges[2]={0,1},pop[6]={0};fusion_source_state_v1 *p=NULL,*q=NULL;
 uint64_t ticket=0;size_t length=0,written=0;fusion_source_ledger_v1 step={0};int cells=0;
 if(fusion_c_source_state_create(1,edges,pop,pop,0,123,&p))return 1;
 if(fusion_c_source_state_begin(p,.1,&ticket))return 2;
 if(fusion_c_source_state_stage(p,ticket,pop,pop,&step))return 3;
 if(fusion_c_source_state_commit(p,ticket))return 4;
 if(fusion_c_source_state_pack_size(p,&length))return 5;
 unsigned char *bytes=malloc(length);if(!bytes)return 6;
 if(fusion_c_source_state_pack(p,bytes,length,&written))return 7;
 if(fusion_c_source_state_unpack(bytes,written,123,&q))return 8;
 if(fusion_c_source_state_cells(q,&cells)||cells!=1)return 9;
 fusion_c_source_state_destroy(p);fusion_c_source_state_destroy(q);free(bytes);return 0;
}
