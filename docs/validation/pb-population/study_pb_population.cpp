#include "fusion_pb_population.h"
#include "fusion_nuclear_data.h"
#include <iostream>
#include <iomanip>
int main(){constexpr double u=1.602176634e-16;fusion_nuclear_mass_v1 p{},b{};fusion_c_nuclear_mass(0,&p);fusion_c_nuclear_mass(5,&b);std::cout<<std::setprecision(17)<<"kind,Ta_or_beam_keV,Tb_keV,component,K,Ea, Eb,Erel,Ecm,probability\n";
for(int kind:{0,1})for(double T:{2.,10.,50.,150.,300.,1000.}){double Tb=kind==0?T*.7:3.;fusion_pb_population_rates_v1 r{};int s=kind==0?fusion_c_pb_thermal_population_rates(1,0,p.mass_kg,b.mass_kg,T*u,Tb*u,.051,1,&r):fusion_c_pb_beam_population_rates(1,0,p.mass_kg,b.mass_kg,T*u,Tb*u,.051,1,&r);if(s){std::cerr<<"status "<<s<<" kind "<<kind<<" T "<<T<<'\n';return 1;}fusion_beam_window_v1* v[]={&r.total,&r.alpha0_peak,&r.narrow_remainder,&r.other_remainder,&r.continuum_extrapolated};for(int i=0;i<5;++i)std::cout<<kind<<','<<T<<','<<Tb<<','<<i<<','<<v[i]->resolved_reactivity_m3_s<<','<<v[i]->projectile_energy_reactivity_J_m3_s<<','<<v[i]->target_energy_reactivity_J_m3_s<<','<<v[i]->relative_energy_reactivity_J_m3_s<<','<<v[i]->cm_energy_reactivity_J_m3_s<<','<<v[i]->resolved_pair_probability<<'\n';}
}
