#include <fusion_thermal_birth.h>
#include <math.h>
#include <stdio.h>
int main(void) {
 const double MeV=1.602176634e-13;
 fusion_thermal_birth_options_v1 o={MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,8,8};
 double edges[201],birth[1400];fusion_thermal_birth_v1 r;
 for(int j=0;j<=200;++j)edges[j]=20*MeV*j/200;
 int status=fusion_c_thermal_birth_grid(3,.001*MeV,&o,200,edges,birth,&r);
 if(status||r.reactivity_m3_s<=0||fabs(r.relative_rate_discrepancy)>1e-6)return 1;
 for(int id=0;id<7;++id){double N=r.below_number_m3_s[id]+r.above_number_m3_s[id];
  for(int j=0;j<200;++j)N+=birth[id*200+j];
  if(fabs(N-((id==4||id==6)?r.reactivity_m3_s:0))>1e-12*r.reactivity_m3_s)return 2;
 }
 puts("Installed C11 thermal birth consumer passed");return 0;
}
