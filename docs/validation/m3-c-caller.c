#include "fusion_coulomb.h"
#include "fusion_kinetics.h"
#include <math.h>
int main(void) {
  double edge[2]={1e-16,3e-16},old[1]={4e20},birth[1]={3e20},escape[1]={.2},trial[1];
  fusion_kinetic_ledger_v1 ledger;
  int status=fusion_c_energy_fp_trial(1,0,.5,edge,old,0,0,birth,escape,.4,trial,0,&ledger);
  if(status || fabs(trial[0]/((4e20+.5*3e20)/1.3)-1)>1e-14) return 1;
  fusion_maxwellian_bath_v1 bath={1e20,3.3435837768e-27,1,1.602176634e-16,15};
  fusion_coulomb_energy_v1 coefficient;
  status=fusion_c_coulomb_energy(20*1.602176634e-16,6.644657345e-27,2,&bath,&coefficient);
  return status || !(coefficient.diffusion_J2_s>0 && coefficient.mean_energy_rate_J_s<0);
}
