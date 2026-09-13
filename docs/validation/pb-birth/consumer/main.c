#include "fusion_pb_birth.h"
#include <math.h>
int main(void){const double u=1.602176634e-13;double e[3]={0,4*u,8*u},b[2];fusion_pb_birth_v1 r;int s=fusion_c_pb_cm_source_grid(8.84*u,.09184*u,.05,.95,13,1,.001*u,.76,0,32,32,2,e,b,&r);if(s)return s;return fabs(r.mapped_number+r.below_number+r.above_number-3)>1e-12||fabs((r.mapped_energy_J+r.below_energy_J+r.above_energy_J)/(8.84*u)-1)>1e-12;}
