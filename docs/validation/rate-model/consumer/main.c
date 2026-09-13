#include <fusion_rate_model.h>
#include <stddef.h>
#include <math.h>
_Static_assert(sizeof(fusion_rate_model_v1)==4*sizeof(fusion_beam_window_v1),"C ABI layout");
int main(void){
 fusion_rate_model_v1 r;
 int status=fusion_c_thermal_pair_maxwellian_model(FUSION_DT_ALPHAN,
  FUSION_ENDPOINT_S,FUSION_PB_LOW_TB,3.3435837768e-27,5.0073567512e-27,
  3*1.602176634e-16,3*1.602176634e-16,&r);
 return status || !(isfinite(r.total.resolved_reactivity_m3_s) &&
                    r.total.resolved_reactivity_m3_s>0 && r.total.domain_incomplete==0);
}
