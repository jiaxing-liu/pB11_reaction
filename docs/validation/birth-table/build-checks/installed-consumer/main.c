#include "fusion_coupled_thermal.h"
#include <math.h>
#include <stdio.h>
int main(void){
 enum{n=200};const double keV=1.602176634e-16,MeV=1.602176634e-13,ne=1e22;
 double edges[n+1]={0},zero[6*n]={0},s[6*n],t[6*n],logs[48],N[6]={0},nextN[6],z2[6]={1,1,1,4,4,25};
 int low=n/3;for(int j=1;j<=low;++j)edges[j]=1e-10*keV*pow(1e12,(double)(j-1)/(low-1));for(int j=low+1;j<=n;++j)edges[j]=(100.+24900.*(j-low)/(n-low))*keV;for(int j=0;j<48;++j)logs[j]=15;
 fusion_coupled_thermal_options_v1 o={0};o.birth=(fusion_thermal_birth_options_v1){2.5*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,4,4};o.channels[3]=1;o.max_source_rate_error=o.max_source_debit_error=1e-5;o.handoff_max_L1=o.handoff_max_mean_error=.001;
 fusion_birth_table_control_v1 ctl={.001,.001,.001,.001,1e-5,1e-5,64,512,12};fusion_birth_table_v1 *table=NULL;
 if(fusion_c_birth_table_create(3,19.5*keV,23*keV,&o.birth,&ctl,n,edges,&table))return 1;
 fusion_birth_table_info_v1 info={0};if(fusion_c_birth_table_info(table,&info)||info.cells!=n||info.knots<2)return 2;
 const fusion_birth_table_v1*tables[5]={0};tables[3]=table;fusion_inert_ion_v1 carbon={.001*ne,12*1.66053906892e-27,36};N[1]=N[2]=(ne-6*carbon.density_m3)/2;
 fusion_coupled_thermal_v1 result={0};
 int st=fusion_c_coupled_thermal_table_trial(.000125,&o,tables,n,edges,N,1.5*ne*5*keV,1.5*(N[1]+N[2]+carbon.density_m3)*20*keV,ne,z2,1,&carbon,logs,zero,zero,zero,zero,nextN,s,t,&result);
 fusion_c_birth_table_destroy(table);
 if(st||result.ledger.events_m3[3]<=0||result.ledger.neutron_number_m3<=0||result.inert_ion_heat_J_m3[4]<=0||nextN[1]>=N[1])return 3;
 puts("Installed C11 table-backed DT coupled trial passed");return 0;
}
