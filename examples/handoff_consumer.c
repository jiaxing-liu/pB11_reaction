/* Installed C11 interface example: trace candidate near a 10-keV Maxwellian.
 * Sources are trial AMOUNTS; commit them once only after the bath accepts
 * the signed correction. This example does not advance a plasma background. */
#include <fusion_handoff.h>
#include <math.h>
#include <stdio.h>
_Static_assert(sizeof(fusion_handoff_ledger_v1)==13*sizeof(double),"handoff ABI");
_Static_assert(sizeof(fusion_maxwellian_grid_v1)==4*sizeof(double),"grid ABI");
int main(void){
 enum {n=1000};const double T=1.602176634e-15;
 double edges[n+1],q[n],old[n],trial[n];
 for(int i=0;i<=n;i++)edges[i]=20*T*i/n;
 fusion_maxwellian_grid_v1 grid={0};
 int status=fusion_c_maxwellian_energy_grid(n,T,edges,q,&grid);
 if(status)return status;
 for(int i=0;i<n;i++)old[i]=1e12*q[i];
 int projected=0;fusion_handoff_ledger_v1 out={0};
 status=fusion_c_maxwellian_handoff_trial(n,T,.001,.001,edges,old,trial,&projected,&out);
 if(status||!projected)return 10+status;
 if(fabs((out.fluid_energy_J_m3+out.bath_energy_correction_J_m3)/out.initial_energy_J_m3-1)>1e-12)return 20;
 printf("fluid N %.9g m^-3, fluid U %.9g J/m^3, bath correction %+.9g J/m^3\n",out.fluid_number_m3,out.fluid_energy_J_m3,out.bath_energy_correction_J_m3);
 return 0;
}
