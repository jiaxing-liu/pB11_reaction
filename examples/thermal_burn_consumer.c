#include "fusion_thermal_burn.h"
#include <stdio.h>
#include <math.h>
int main(void) {
 /* Supplied-coefficient numerical example, not an EXL plasma prediction. */
 double n[6]={1e20,0,0,0,0,1e20},u[6]={1e5,0,0,0,0,1e5};
 double k[5]={1e-22,0,0,0,0},ea[5]={1e-15,0,0,0,0},eb[5]={1e-15,0,0,0,0};
 double nn[6],uu[6];fusion_thermal_burn_v1 ledger;
 if(fusion_c_thermal_burn_trial(1,n,u,k,ea,eb,nn,uu,&ledger)!=PB11_STATUS_OK)return 1;
 if(fabs(ledger.fast_product_birth_m3[FUSION_HELIUM4]-3*ledger.events_m3[0])>1e-12*n[0])return 2;
 printf("pB events %.9g m^-3, remaining proton %.9g m^-3, FAST alpha births %.9g m^-3\n",
  ledger.events_m3[0],nn[0],ledger.fast_product_birth_m3[4]);
 printf("Removed reactant energy %.9g J/m^3; product owner still adds consistent Q.\n",
  ledger.reactant_removed_energy_J_m3[0]+ledger.reactant_removed_energy_J_m3[5]);
 return 0;
}
