#include "fusion_alpha_spectrum.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
int main(int argc,char**argv){
 int n=argc>1?std::stoi(argv[1]):128;
 constexpr double mev=1.602176634e-13;
 std::vector<double>edges(101),birth(100);
 for(int i=0;i<=100;++i)edges[i]=.06*i*mev;
 std::cout<<std::setprecision(17)<<"mode,energy_keV,number_per_event\n";
 const double kval=argc>2?std::stod(argv[2]):.76;
 const std::vector<int> modes=argc>2?std::vector<int>{13}:std::vector<int>{1,3,13};
 for(int mode:modes){
  fusion_alpha_spectrum_v1 r{};int s=fusion_c_alpha_spectrum_grid(mode,9.3*mev,.001*mev,kval,.67*2*std::acos(-1.),n,n,100,edges.data(),birth.data(),&r);
  if(s){std::cerr<<"status "<<s<<'\n';return 1;}
  std::cerr<<mode<<" mapped="<<r.mapped_number<<" spill="<<r.below_number+r.above_number<<" e="<<r.mapped_energy_J/mev<<" n1="<<r.l1_normalization_J2<<" n3="<<r.l3_normalization_J2<<'\n';
  for(int i=0;i<100;++i)std::cout<<mode<<','<<30+60*i<<','<<birth[i]<<'\n';
 }
}
