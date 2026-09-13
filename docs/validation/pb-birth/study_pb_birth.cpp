#include "fusion_pb_birth.h"
#include <iostream>
#include <iomanip>
#include <vector>
int main(){constexpr double mev=1.602176634e-13;std::vector<double> e(101),b(100);for(int i=0;i<=100;++i)e[i]=i*.06*mev;std::cout<<std::setprecision(17)<<"alpha0_fraction,E_MeV,N_per_event\n";for(double f:{0.,.05,1.}){fusion_pb_birth_v1 r{};int s=fusion_c_pb_cm_source_grid(8.84*mev,.09184*mev,f,1-f,13,1,.001*mev,.76,0,1024,1024,100,e.data(),b.data(),&r);if(s){std::cerr<<"status "<<s<<'\n';return 1;}for(int i=0;i<100;++i)std::cout<<f<<','<<.03+.06*i<<','<<b[i]<<'\n';std::cerr<<std::setprecision(17)<<f<<' '<<r.mapped_number+r.below_number+r.above_number<<' '<<(r.mapped_energy_J+r.below_energy_J+r.above_energy_J)/mev<<' '<<r.below_number<<' '<<r.above_number<<' '<<r.primary_alpha0_energy_J/mev<<' '<<r.secondary_min_J/mev<<' '<<r.secondary_max_J/mev<<'\n';}}
