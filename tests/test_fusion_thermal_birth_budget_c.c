/* Real C translation unit: C linkage and documented error clearing only. */
#include "fusion_thermal_birth_budget.h"
#include <stddef.h>
#include <string.h>
_Static_assert(sizeof(((fusion_thermal_birth_budget_report_v1*)0)->joint_nodes)==sizeof(uint64_t),"counts are uint64_t");
_Static_assert(offsetof(fusion_thermal_birth_budget_policy_v1,number_L1_tolerance)>=4*sizeof(int),"policy layout");
int main(void){
 fusion_thermal_birth_budget_policy_v1 p={2,4,4,0,.01,.01};
 fusion_thermal_birth_v1 out;
 fusion_thermal_birth_budget_report_v1 report;
 double birth[7],edges[2]={0,1};
 unsigned char zero_out[sizeof(out)]={0},zero_report[sizeof(report)]={0};
 memset(&out,0x5a,sizeof(out));memset(&report,0x5a,sizeof(report));
 for(int i=0;i<7;++i)birth[i]=99;
 if(fusion_c_pb_thermal_birth_budget_grid(1,0,&p,1,edges,birth,&out,&report)==0)return 1;
 if(memcmp(&out,zero_out,sizeof(out))||memcmp(&report,zero_report,sizeof(report)))return 2;
 for(int i=0;i<7;++i)if(birth[i]!=0)return 3;
 return 0;
}
