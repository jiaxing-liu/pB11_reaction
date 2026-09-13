#include "fusion_nuclear_data.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
namespace {
void check(bool x){if(!x)throw std::runtime_error("nuclear birth four-momentum check");}
}
int main(){
 const double mev=1.602176634e-13,dir[]={.6,0,.8};
 const double velocities[2][3]={{0,0,0},{1e7,-2e7,3e7}};
 std::cout<<std::setprecision(17)<<"channel,Ecm_MeV,boosted,Q_MeV,lab_kinetic_gain_MeV,relative_Q_residual,relative_momentum_residual\n";
 for(int channel=1;channel<5;++channel)for(double E:{0.,.01,1.})for(int iv=0;iv<2;++iv){
  fusion_nuclear_channel_v1 r{};check(fusion_c_nuclear_channel(channel,&r)==0);
  fusion_nuclear_mass_v1 a{},b{};
  check(fusion_c_nuclear_mass(r.reactant_ids[0],&a)==0);
  check(fusion_c_nuclear_mass(r.reactant_ids[1],&b)==0);
  fusion_particle_four_vector_v1 reactants[2]{},products[2]{},lab_in[2]{},lab_out[2]{};
  check(fusion_c_two_body_cm(a.mass_kg,b.mass_kg,E*mev,dir,reactants)==0);
  check(fusion_c_nuclear_two_body_cm(channel,E*mev,dir,products)==0);
  fusion_boost_ledger_v1 bi{},bo{};
  check(fusion_c_boost_particles(2,velocities[iv],reactants,lab_in,&bi)==0);
  check(fusion_c_boost_particles(2,velocities[iv],products,lab_out,&bo)==0);
  long double gain=0,p2=0,pscale=0;
  for(int i=0;i<2;++i)gain+=(long double)lab_out[i].kinetic_energy_J-lab_in[i].kinetic_energy_J;
  for(int j=0;j<3;++j){
   long double dp=0;
   for(int i=0;i<2;++i){dp+=(long double)lab_out[i].momentum_kg_m_s[j]-lab_in[i].momentum_kg_m_s[j];
    pscale+=std::abs((long double)lab_out[i].momentum_kg_m_s[j])+std::abs((long double)lab_in[i].momentum_kg_m_s[j]);}
   p2+=dp*dp;
  }
  long double er=std::abs(gain-r.q_J)/r.q_J,pr=std::sqrt(p2)/pscale;
  check(er<2e-11L&&pr<2e-12L);
  std::cout<<channel<<','<<E<<','<<iv<<','<<r.q_J/mev<<','<<(double)(gain/mev)<<','<<(double)er<<','<<(double)pr<<'\n';
 }
}
