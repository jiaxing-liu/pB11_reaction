#include <fusion_kinetic_geometry.h>
#include <stdio.h>
int main(void){
 double edges[2]={0,2},va[1]={2},vb[1]={1},a[2]={0},c[1]={0},k[12]={0},bc[12]={0};
 double s[6]={0,0,0,0,3,0},t[6]={0,0,0,0,1,0},ns[6],nt[6];fusion_transport_ledger_v1 movement;
 int status=fusion_c_kinetic_geometry_trial(1,1,.2,edges,va,vb,a,c,k,bc,bc,s,t,ns,nt,&movement);
 if(status||ns[4]!=6||nt[4]!=2||movement.work_J[4]!=0)return 1;
 puts("PASS installed C11 kinetic geometry consumer");return 0;
}
