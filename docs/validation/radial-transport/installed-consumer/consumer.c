#include <fusion_radial_transport.h>
#include <math.h>
#include <stdio.h>
_Static_assert(sizeof(fusion_radial_ledger_v1)==5*sizeof(double),"C ledger ABI");
int main(void){
 double vo[1]={2},vn[1]={1},a[2]={0,0},k[2]={0,0};
 double old[1]={3},bc[2]={0,0},out[1],face[2];
 fusion_radial_ledger_v1 ledger;
 int s=fusion_c_radial_transport_trial(1,1,.1,vo,vn,a,k,old,bc,out,face,&ledger);
 if(s || out[0]!=6 || ledger.initial_number!=6 || ledger.final_number!=6) return 1;
 puts("PASS installed C11 radial consumer");return 0;
}
