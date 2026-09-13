#include "fusion_thermal_birth.h"
#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace {
using R=long double;
constexpr R c2=299792458.L*299792458.L, pi=3.1415926535897932384626433832795L;
constexpr double MeV=1.602176634e-13;
void need(bool p,const char*s){if(!p)throw std::runtime_error(s);}
// Independent composite midpoint integration in E/T and |V| sqrt(M/2T).
// Does not call the production CM quadrature, parent, boost or mapper.
std::array<R,3> reference(int ch,double T,int resolution){
 fusion_nuclear_channel_v1 r{}; need(!fusion_c_nuclear_channel(ch,&r),"channel");
 fusion_nuclear_mass_v1 a{},b{},p[2]{};
 fusion_c_nuclear_mass(r.reactant_ids[0],&a);fusion_c_nuclear_mass(r.reactant_ids[1],&b);
 for(int j=0;j<2;++j)fusion_c_nuclear_mass(r.product_ids[j],&p[j]);
 R M=R(a.mass_kg)+b.mass_kg,mu=R(a.mass_kg)*b.mass_kg/M;
 R ro[2]={R(p[0].mass_kg)*c2,R(p[1].mass_kg)*c2},rest=ro[0]+ro[1];
 R dx=50.L/resolution,dy=std::sqrt(40.L)/resolution;
 R pref=std::sqrt(8*R(T)/(pi*mu));std::array<R,3> out{};
 for(int i=0;i<resolution;++i){R x=(i+.5L)*dx,E=x*T;double sigma=0;
  need(!fusion_c_cross_section_model(ch,1,0,double(E),&sigma),"cross section");
  R wr=pref*x*std::exp(-x)*sigma*dx;
  for(int k=0;k<resolution;++k){R y=(k+.5L)*dy,C=y*y*T;
   R w=wr*4/std::sqrt(pi)*y*y*std::exp(-y*y)*dy;
   R lab=rest+r.q_J+E+C,pc2=2*M*c2*C;
   R W=std::sqrt(lab*lab-pc2),A=W-rest;
   R gamma=lab/W,gb2=pc2/(W*W);
   out[0]+=w;
   for(int j=0;j<2;++j){R K=A*(A+2*ro[1-j])/(2*W);
    R mid=gamma*(ro[j]+K)-ro[j];
    out[j+1]+=w*(mid*mid+gb2*K*(K+2*ro[j])/3);
   }
  }
 }
 return out;
}
}
int main(){try{
 const double T=.02*MeV;
 fusion_thermal_birth_options_v1 o{1*MeV,40, .09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,24,12,8,8};
 for(int ch=1;ch<5;++ch){
  const int n=4000;std::vector<double> edges(n+1),birth(7*n);
  for(int j=0;j<=n;++j)edges[j]=20*MeV*j/n;
  fusion_thermal_birth_v1 result{};
  int status=fusion_c_thermal_birth_grid(ch,T,&o,n,edges.data(),birth.data(),&result); if(status)std::cerr<<"channel "<<ch<<" status "<<status<<"\n"; need(!status,"source status");
  auto a=reference(ch,T,512),b=reference(ch,T,1024);
  need(std::abs(a[0]/b[0]-1)<2e-7L,"independent rate convergence");
  need(std::abs(R(result.reactivity_m3_s)/b[0]-1)<2e-7L,"independent rate parity");
  fusion_nuclear_channel_v1 r{};fusion_c_nuclear_channel(ch,&r);
  for(int k=0;k<2;++k){int id=r.product_ids[k];R s=0;
   need(result.below_number_m3_s[id]==0&&result.above_number_m3_s[id]==0,"complete grid support");
   for(int j=0;j<n;++j){R center=(R(edges[j])+edges[j+1])/2;s+=birth[id*n+j]*center*center;}
   R reference_error=std::abs(a[k+1]-b[k+1]);
   R spacing=R(edges[1])-edges[0],hat_bound=result.reactivity_m3_s*spacing*spacing/4;
   R uncertainty=4*reference_error+2e-9L*b[k+1];
   // Positive linear interpolation increases E^2 by at most h^2/4 per particle.
   need(s-b[k+1]>=-uncertainty&&s-b[k+1]<=hat_bound+uncertainty,"independent spectral second moment");
   std::cout<<ch<<' '<<id<<" relative second-moment difference "<<double(s/b[k+1]-1)<<" bound "<<double((hat_bound+uncertainty)/b[k+1])<<'\n';
  }
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
