#include "fusion_coupled_thermal.h"
#include <stdio.h>
_Static_assert(sizeof(fusion_thermal_increment_v1)==8*sizeof(double),"increment ABI");
int main(void){
 fusion_source_ledger_v1 ledger={0};
 fusion_thermal_increment_v1 out={0};double inert[6]={0};
 ledger.thermal_consumed_number_m3[1]=7;
 ledger.heat_to_bath_J_m3[0]=1e-30;
 inert[3]=-0.5;
 if(fusion_c_coupled_thermal_increment(&ledger,inert,&out)!=0)return 1;
 if(out.thermal_number_m3[1]!=-7||out.electron_energy_J_m3!=1e-30||out.ion_energy_J_m3!=-0.5)return 2;
 puts("C11 thermal-increment ABI and signed weak-source call passed");return 0;
}
