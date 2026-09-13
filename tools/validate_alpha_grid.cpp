// Standalone numerical sensitivity diagnostic; stdout is CSV.
#include "fusion_alpha_spectrum.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
int main(){
 constexpr double mev=1.602176634e-13,phase=.67*2*3.141592653589793;
 constexpr int cells=32;
 std::vector<double> edges(cells+1),birth(cells);
 for(int i=0;i<=cells;++i)edges[i]=8.*mev*i/cells;
 std::cout<<std::setprecision(17)<<"mode,A_MeV,nq,ncos,cutoff_MeV,normalization_J2,l1_norm_J2,l3_norm_J2,below_N,below_E_MeV,above_N,above_E_MeV,N_residual,E_residual_MeV,pruned_events";
 for(int i=0;i<cells;++i)std::cout<<",N_bin"<<i;
 std::cout<<'\n';
 for(double A:{8.829,9.3,12.})for(int mode:{1,2,3,13})for(int n:{32,64,128,256})for(double cutoff:{.001,.002,.004,.01}){
  // Resolve all shapes at finest cut; cutoff sweep only at finest quadrature.
  if(n!=256 && cutoff!=.001)continue;
  fusion_alpha_spectrum_v1 o{};
  int status=fusion_c_alpha_spectrum_grid(mode,A*mev,cutoff*mev,.76,phase,n,n,cells,edges.data(),birth.data(),&o);
  if(status){std::cerr<<"status="<<status<<" A="<<A<<" mode="<<mode<<" n="<<n<<" cut="<<cutoff<<'\n';return 1;}
  std::cout<<mode<<','<<A<<','<<n<<','<<n<<','<<cutoff<<','<<o.normalization_J2<<','<<o.l1_normalization_J2<<','<<o.l3_normalization_J2<<','<<o.below_number<<','<<o.below_energy_J/mev<<','<<o.above_number<<','<<o.above_energy_J/mev<<','<<o.number_residual<<','<<o.energy_residual_J/mev<<','<<o.pruned_events;
  for(double x:birth)std::cout<<','<<x;
  std::cout<<'\n'<<std::flush;
 }
}
