#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
int main(){
 constexpr double kev=1.602176634e-16;
 std::cout<<std::setprecision(17)<<"channel,T_keV,K_fit,K_low,K_high,E_fit,E_low,E_high\n";
 for(int ch=0;ch<5;++ch){
  fusion_nuclear_channel_v1 c{};fusion_nuclear_mass_v1 a{},b{};
  if(fusion_c_nuclear_channel(ch,&c)||fusion_c_nuclear_mass(c.reactant_ids[0],&a)||fusion_c_nuclear_mass(c.reactant_ids[1],&b))throw std::runtime_error("mass");
  for(double T:{.01,.03,.05,.1,.2,.5,1.,3.,10.,30.,100.,190.,500.}){
   fusion_rate_model_v1 r{};int status=fusion_c_thermal_pair_maxwellian_model(ch,1,0,a.mass_kg,b.mass_kg,T*kev,T*kev,&r);
   if(status){std::cerr<<ch<<' '<<T<<" status "<<status<<'\n';return 1;}
   std::cout<<ch<<','<<T<<','<<r.fit.resolved_reactivity_m3_s<<','<<r.below.resolved_reactivity_m3_s<<','<<r.above.resolved_reactivity_m3_s;
   for(const auto*p:{&r.fit,&r.below,&r.above})std::cout<<','<<(p->resolved_reactivity_m3_s>0?p->relative_energy_reactivity_J_m3_s/p->resolved_reactivity_m3_s/kev:0);
   std::cout<<'\n';
  }
 }
}
