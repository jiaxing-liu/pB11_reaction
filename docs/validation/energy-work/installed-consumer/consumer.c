#include <fusion_energy_work.h>
#include <math.h>
#include <stdio.h>
_Static_assert(sizeof(fusion_energy_work_ledger_v1)==88,"ledger ABI");
int main(void){double e[]={1,3},old[]={4},out[1],face[2];fusion_energy_work_ledger_v1 b;
 int s=fusion_c_energy_work_trial(1,.5,1,e,old,out,face,&b);
 if(s || fabs(out[0]-3.2)>1e-13 || fabs(b.work_on_particles_J_m3+.8)>1e-13)return 1;
 puts("PASS installed C11 work consumer");return 0;}
