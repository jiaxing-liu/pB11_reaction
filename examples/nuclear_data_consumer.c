#include "fusion_nuclear_data.h"
#include <stdio.h>
#include <math.h>
int main(void){
 fusion_nuclear_channel_v1 r;fusion_particle_four_vector_v1 p[2];
 const double direction[3]={0,0,1},Ecm=1.602176634e-15;
 if(fusion_c_nuclear_channel(FUSION_DT_ALPHAN,&r))return 1;
 if(fusion_c_nuclear_two_body_cm(FUSION_DT_ALPHAN,Ecm,direction,p))return 2;
 if(fabs(p[0].kinetic_energy_J+p[1].kinetic_energy_J-Ecm-r.q_J)>1e-12*r.q_J)return 3;
 printf("DT Q %.12g MeV, alpha %.12g MeV, neutron %.12g MeV at Ecm=10keV\n",
 r.q_J/1.602176634e-13,p[0].kinetic_energy_J/1.602176634e-13,p[1].kinetic_energy_J/1.602176634e-13);
 return 0;
}
