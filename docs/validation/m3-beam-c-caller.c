#include "fusion_beam.h"
#include <math.h>
#include <stdio.h>
int main(void) {
    const double unit=1.66053906892e-27,k=1.602176634e-16;
    fusion_beam_window_v1 result;
    int status=fusion_c_thermal_pair_maxwellian_window(FUSION_DT_ALPHAN,
        2*unit,3*unit,5*k,12*k,&result);
    if(status || result.domain_incomplete!=1 || !(result.resolved_reactivity_m3_s>0)) return 1;
    double sum=result.projectile_energy_reactivity_J_m3_s+result.target_energy_reactivity_J_m3_s;
    if(fabs(sum-result.relative_energy_reactivity_J_m3_s-result.cm_energy_reactivity_J_m3_s)>1e-10*sum) return 2;
    printf("DT resolved-window K = %.12e m^3/s\n",result.resolved_reactivity_m3_s);
    printf("Reaction-conditioned Ea,Eb = %.9f, %.9f keV\n",
        result.projectile_energy_reactivity_J_m3_s/result.resolved_reactivity_m3_s/k,
        result.target_energy_reactivity_J_m3_s/result.resolved_reactivity_m3_s/k);
    printf("Unresolved relative-pair probability = %.12e; incomplete flag = %d\n",
        result.unresolved_pair_probability,result.domain_incomplete);
    return 0;
}
