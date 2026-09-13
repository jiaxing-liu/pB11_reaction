#include "fusion_reaction_event.h"
#include <math.h>
int main(void){double pa[3]={1e-21,0,0},pb[3]={0,0,0},d[3]={0,0,1};fusion_reaction_parent_v1 p;fusion_boost_ledger_v1 l;fusion_particle_four_vector_v1 o[3];int s=fusion_c_reaction_lab_event(0,0,pa,pb,d,91.84*1.602176634e-16,.2,1,o,&p,&l);if(s)return s;return fabs(l.energy_residual_J)>1e-10*l.expected_kinetic_energy_J;}
