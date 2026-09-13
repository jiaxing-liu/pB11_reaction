#include "fusion_reaction_event.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
namespace{using R=long double;constexpr double u=1.602176634e-16;constexpr R c2=299792458.L*299792458.L;
void need(bool ok,const char*s){if(!ok)throw std::runtime_error(s);}
bool near(R a,R b,R eps=2e-11L){return std::abs(a-b)<=eps*std::max(std::abs(a),std::abs(b));}
}
int main(){try{
 std::mt19937_64 gen(271828);std::normal_distribution<double>N;std::uniform_real_distribution<double>U;
 for(int ch=0;ch<5;++ch){fusion_nuclear_channel_v1 reaction{};need(!fusion_c_nuclear_channel(ch,&reaction),"channel");fusion_nuclear_mass_v1 masses[2]{};for(int i=0;i<2;++i)need(!fusion_c_nuclear_mass(reaction.reactant_ids[i],&masses[i]),"mass");
  for(int z=0;z<40;++z){double pa[3]{},pb[3]{},d[3]{};double norm=0;for(int j=0;j<3;++j){pa[j]=z?N(gen)*std::sqrt(masses[0].mass_kg*(1+1000*U(gen))*u):0;pb[j]=z?N(gen)*std::sqrt(masses[1].mass_kg*(1+1000*U(gen))*u):0;d[j]=N(gen);norm+=d[j]*d[j];}for(double&x:d)x/=std::sqrt(norm);
   double cos=2*U(gen)-1,phi=6.28*U(gen);fusion_reaction_parent_v1 parents[2]{};
   for(int mode=0;mode<2;++mode){fusion_particle_four_vector_v1 out[3]{};fusion_boost_ledger_v1 l{};need(!fusion_c_reaction_lab_event(ch,mode,pa,pb,d,ch==0?91.84*u:0,ch==0?cos:0,ch==0?phi:0,out,&parents[mode],&l),"event status");const auto&p=parents[mode];need(p.product_count==reaction.product_count&&p.convention==mode,"metadata");
    R E=0,psum[3]{},pscale=0;for(int i=0;i<p.product_count;++i){fusion_nuclear_mass_v1 m{};need(!fusion_c_nuclear_mass(reaction.product_ids[i],&m)&&out[i].mass_kg==m.mass_kg,"canonical products");R P2=0;for(int j=0;j<3;++j){P2+=R(out[i].momentum_kg_m_s[j])*out[i].momentum_kg_m_s[j];psum[j]+=out[i].momentum_kg_m_s[j];pscale+=std::abs(R(out[i].momentum_kg_m_s[j]));}
     R K=out[i].kinetic_energy_J;need(K>=0&&near(P2*c2,K*(K+2*R(m.mass_kg)*c2)),"product mass shell");E+=K;}
    R input=p.reactant_kinetic_J[0]+R(p.reactant_kinetic_J[1])+reaction.q_J;need(near(E,input),"lab energy closure");for(int j=0;j<3;++j)need(std::abs(psum[j]-R(pa[j])-pb[j])<=2e-12L*pscale,"momentum closure");need(near(l.expected_kinetic_energy_J,input),"ledger expected");if(ch)need(out[2].mass_kg==0&&out[2].kinetic_energy_J==0,"unused output clears");if(z==0)need(near(p.available_cm_energy_J,reaction.q_J),"stationary Q");
   }
   R delta=0;for(int i=0;i<2;++i){need(parents[0].classical_kinetic_J[i]>=parents[0].on_shell_kinetic_J[i],"classical upper energy");delta+=R(parents[0].reactant_kinetic_J[i])-parents[1].reactant_kinetic_J[i];}
   if(z)need(near(delta,parents[0].classical_minus_on_shell_J,2e-9L),"convention difference");
  }
 }
 {double zero[3]{},bad[3]{1,1,1},d[3]{0,0,1};fusion_reaction_parent_v1 p{};fusion_boost_ledger_v1 l{};fusion_particle_four_vector_v1 o[3]{};need(fusion_c_reaction_lab_event(3,0,zero,zero,d,u,0,0,o,&p,&l)==PB11_STATUS_OUT_OF_RANGE,"unused q rejected");p.available_cm_energy_J=1;o[0].mass_kg=1;l.output_kinetic_energy_J=1;need(fusion_c_reaction_lab_event(0,0,zero,zero,bad,u,0,0,o,&p,&l)==PB11_STATUS_OUT_OF_RANGE&&p.available_cm_energy_J==0&&o[0].mass_kg==0&&l.output_kinetic_energy_J==0,"bad direction clears");need(fusion_c_reaction_parent(0,2,zero,zero,&p)==PB11_STATUS_INVALID_ARGUMENT,"bad mode");bad[0]=std::numeric_limits<double>::quiet_NaN();need(fusion_c_reaction_parent(0,0,bad,zero,&p)==PB11_STATUS_INVALID_ARGUMENT,"NaN");need(fusion_c_reaction_lab_event(0,0,zero,zero,d,u,1.01,0,o,&p,&l)==PB11_STATUS_OUT_OF_RANGE,"cosine guard");}
 std::cout<<"400 laboratory events preserve energy momentum and product mass shells\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
