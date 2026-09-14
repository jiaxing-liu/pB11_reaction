#include <fusion_fast_moments.h>
#include <math.h>
int main(void){double e[2]={0,2e-15},s[6]={0,0,0,0,3,0},t[6]={0};fusion_fast_moments_v1 m;
if(sizeof(m)!=120||fusion_c_fast_moments(1,e,s,t,&m))return 1;
return m.charge_number_m3!=6||m.charge_squared_number_m3!=12||fabs(m.pressure_Pa-2e-15)>1e-28;}
