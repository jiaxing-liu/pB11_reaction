#include "fusion_pb_population.h"
#include "fusion_nuclear_data.h"
#include <math.h>
int main(void){fusion_nuclear_mass_v1 p,b;fusion_pb_population_rates_v1 r;if(fusion_c_nuclear_mass(0,&p)||fusion_c_nuclear_mass(5,&b))return 1;int s=fusion_c_pb_thermal_population_rates(1,0,p.mass_kg,b.mass_kg,1.602176634e-15,1.1215236438e-15,.051,1,&r);if(s)return s;return fabs((r.alpha0_peak.resolved_reactivity_m3_s+r.narrow_remainder.resolved_reactivity_m3_s+r.other_remainder.resolved_reactivity_m3_s)/r.total.resolved_reactivity_m3_s-1)>1e-8;}
