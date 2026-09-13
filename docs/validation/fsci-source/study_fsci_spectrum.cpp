#include "fusion_alpha_spectrum.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
int main(int argc,char**argv){
 const int n=argc>1?std::stoi(argv[1]):256;constexpr double mev=1.602176634e-13;
 const double cutoff=(argc>2?std::stod(argv[2]):1)*.001*mev;
 std::vector<double> edges(101),birth(100);for(int i=0;i<=100;++i)edges[i]=.06*i*mev;
 std::cout<<std::setprecision(17)<<"mode,policy,A_MeV,energy_keV,number_per_event\n";
 const double k=argc>3?std::stod(argv[3]):.76;
 const std::vector<std::pair<int,int>> cases=argc>3?std::vector<std::pair<int,int>>{{13,0},{13,1}}:std::vector<std::pair<int,int>>{{2,0},{2,1},{13,1}};
 for(auto item:cases){
  int mode=item.first,policy=item.second;double A=mode==2?8.84:9.3;fusion_alpha_spectrum_v1 r{};
  int status=fusion_c_alpha_spectrum_model_grid(mode,policy,A*mev,cutoff,k,.67*2*std::acos(-1.),n,n,100,edges.data(),birth.data(),&r);
  if(status){std::cerr<<"API status "<<status<<'\n';return 1;}
  std::cerr<<std::setprecision(17)<<mode<<' '<<policy<<' '<<A<<' '<<r.mapped_number<<' '<<r.below_number<<' '<<r.above_number<<' '<<r.mapped_energy_J/mev<<' '<<r.below_energy_J/mev<<' '<<r.above_energy_J/mev<<' '<<r.pruned_events<<' '<<r.l1_normalization_J2<<' '<<r.l3_normalization_J2<<'\n';
  for(int i=0;i<100;++i)std::cout<<mode<<','<<policy<<','<<A<<','<<30+60*i<<','<<birth[i]<<'\n';
 }
}
