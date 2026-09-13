#include <fusion_two_component.h>
#include <math.h>
_Static_assert(sizeof(fusion_two_component_ledger_v1)==14*sizeof(double),"C ABI ledger");
int main(void){
 double edges[2]={0,2},s[1]={4},t[1]={1},bs[1]={3},bt[1]={.5},esc[1]={.2},tr[1]={.4},sn[1],tn[1];
 fusion_two_component_ledger_v1 out={0};
 int status=fusion_c_two_component_trial(1,0,.5,edges,s,t,0,0,bs,bt,esc,tr,sn,tn,0,&out);
 double wantS=5.5/1.3,wantT=(1.25+.2*wantS)/1.1;
 if(status||fabs(sn[0]-wantS)>1e-14||fabs(tn[0]-wantT)>1e-14)return 1;
 if(out.total.thermalized_number_m3!=0||fabs(out.transferred_number_m3-.2*wantS)>1e-14)return 2;
 return 0;
}
