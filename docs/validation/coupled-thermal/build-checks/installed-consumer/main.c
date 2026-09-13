#include <fusion_coupled_thermal.h>
#include <math.h>
#include <stdio.h>
int main(void){
 const int n=256;const double keV=1.602176634e-16,MeV=1.602176634e-13;
 fusion_coupled_thermal_options_v1 o={0};
 o.birth=(fusion_thermal_birth_options_v1){2.5*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,8,8};
 o.channels[3]=1;o.max_source_rate_error=o.max_source_debit_error=1e-5;
 o.handoff_max_L1=o.handoff_max_mean_error=.001;
 double e[257],N[6]={0,5e19,5e19,0,0,0},Z2[6]={1,1,1,4,4,25},logs[42];
 double s[1536]={0},t[1536]={0},birth[1536]={0},escape[1536]={0},sn[1536],tn[1536],Nn[6];
 e[0]=0;for(int j=1;j<=64;++j)e[j]=1e-12*keV*pow(1e14,(j-1)/63.);
 for(int j=65;j<=n;++j)e[j]=(100.+24900.*(j-64)/192.)*keV;
 for(int j=0;j<42;++j)logs[j]=15;
 double Ue=1.5e20*5*keV,Ui=1.5e20*20*keV;fusion_coupled_thermal_v1 r;
 int status=fusion_c_coupled_thermal_trial(1e-4,&o,n,e,N,Ue,Ui,1e20,Z2,0,NULL,logs,s,t,birth,escape,Nn,sn,tn,&r);
 if(status||!(r.ledger.events_m3[3]>0)||!(Nn[1]<N[1])||!(r.electron_energy_J_m3>Ue))return 1;
 double A=0;for(int j=0;j<n;++j)A+=sn[4*n+j]+tn[4*n+j];
 if(fabs(A/r.ledger.neutron_number_m3-1)>1e-10)return 2;
 puts("Installed C11 coupled DT birth and thermal feedback passed");return 0;
}
