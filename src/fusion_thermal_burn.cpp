#include "fusion_thermal_burn.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
namespace {
using R=long double;
constexpr int ns=FUSION_SPECIES_COUNT,nc=FUSION_CHANNEL_COUNT;
constexpr int a[nc]={0,1,1,1,1},b[nc]={5,1,1,2,3};
bool nonnegative(double x){return std::isfinite(x)&&x>=0;}
bool put(R x,double& v){
 if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;
 v=static_cast<double>(x);return true; // rounding a positive subnormal to zero is allowed
}
// Monotone normalized D equation. Evaluate losses without subtracting nearly
// equal populations; its unique root is in (0,1] for finite coefficients.
bool solve_D(R D,R T,R H,R dt,const double*K,R& x){
 x=1;if(D==0)return true;
 R dd=dt*(R(K[1])+K[2])*D,at=dt*K[3]*D,ah=dt*K[4]*D;
 R bt=T/D,bh=H/D,stiff=dd+bt*at+bh*ah;
 if(!std::isfinite(stiff))return false;
 if(stiff==0)return true;
 R lo=1/(1+stiff),hi=2/(1+std::sqrt(1+4*dd));
 auto evaluate=[&](R y){
  R zt=at*y,zh=ah*y;
  return y+dd*y*y+bt*zt/(1+zt)+bh*zh/(1+zh)-1;
 };
 x=std::sqrt(lo)*std::sqrt(hi);
 for(int i=0;i<256;++i){
  R f=evaluate(x);
  if(!std::isfinite(f))return false;
  if(std::abs(f)<=8*std::numeric_limits<R>::epsilon())return true;
  if(f>0)hi=x;else lo=x;
  R next;
  if(hi>4*lo)next=std::sqrt(lo)*std::sqrt(hi);
  else{
   R zt=1+at*x,zh=1+ah*x;
   R df=1+2*dd*x+bt*at/(zt*zt)+bh*ah/(zh*zh);
   next=x-f/df;
   if(!(next>lo&&next<hi)||!std::isfinite(next))next=(lo+hi)/2;
  }
  if(next==x)return std::abs(f)<=2e-15L;
  x=next;
 }
 return false;
}
}
extern "C" int fusion_c_thermal_burn_trial(double dt,const double*oldN,
 const double*oldU,const double*K,const double*Ea,const double*Eb,
 double*newN,double*newU,fusion_thermal_burn_v1*out){
 if(newN)std::fill(newN,newN+ns,0.);
 if(newU)std::fill(newU,newU+ns,0.);
 if(out)*out={};
 if(!newN||!newU||!out)return PB11_STATUS_NULL_OUTPUT;
 if(!oldN||!oldU||!K||!Ea||!Eb||!std::isfinite(dt)||dt<=0)return PB11_STATUS_INVALID_ARGUMENT;
 for(int s=0;s<ns;++s)if(!nonnegative(oldN[s])||!nonnegative(oldU[s])||
    (oldN[s]==0&&oldU[s]!=0))return PB11_STATUS_INVALID_ARGUMENT;
 for(int c=0;c<nc;++c)if(!nonnegative(K[c])||!nonnegative(Ea[c])||!nonnegative(Eb[c])||
    (K[c]==0&&(Ea[c]!=0||Eb[c]!=0)))return PB11_STATUS_INVALID_ARGUMENT;
 try{
  std::array<R,ns>N{},loss{},debit{};for(int s=0;s<ns;++s)N[s]=oldN[s];
  fusion_thermal_burn_v1 r{};
  // Distinct p/B pair: cancellation-free positive quadratic root of the
  // surviving minor population, including equal populations and K=0.
  R m=std::min(N[0],N[5]),delta=std::abs(N[0]-N[5]),A=R(dt)*K[0];
  if(m>0&&A>0){
   R v=1+A*delta;
   R minor=2*m/(v+std::sqrt(v*v+4*A*m));
   if(!std::isfinite(minor)||minor<0)return PB11_STATUS_NUMERICAL_FAILURE;
   if(N[0]<=N[5]){N[0]=minor;N[5]=delta+minor;}
   else{N[5]=minor;N[0]=delta+minor;}
  }
  R x;
  if(!solve_D(N[1],N[2],N[3],dt,K,x))return PB11_STATUS_NUMERICAL_FAILURE;
  N[1]*=x;
  N[2]/=1+R(dt)*K[3]*N[1];N[3]/=1+R(dt)*K[4]*N[1];
  // Compute events from the implicit populations, not old-new cancellation.
  // Round event amounts once; all product/energy ledgers use those amounts.
  for(int c=0;c<nc;++c){
   R events=R(dt)*K[c]*N[a[c]]*N[b[c]]/(a[c]==b[c]?2:1);
   if(!put(events,r.events_m3[c]))return PB11_STATUS_NUMERICAL_FAILURE;
   R e=r.events_m3[c];loss[a[c]]+=e;loss[b[c]]+=e;
   debit[a[c]]+=e*Ea[c];debit[b[c]]+=e*Eb[c];
  }
  fusion_particle_sources_v1 counts{};
  int status=fusion_c_particle_sources(r.events_m3,&counts);if(status)return status;
  for(int s=0;s<ns;++s)r.fast_product_birth_m3[s]=counts.product_birth[s];
  r.neutron_birth_m3=counts.neutron_birth;
  std::array<double,ns>numbers{},energies{};
  for(int s=0;s<ns;++s){
   R U=R(oldU[s])-debit[s];
   if(U<0||N[s]<0||!put(N[s],numbers[s])||!put(U,energies[s])||
     !put(loss[s],r.reactant_removed_m3[s])||!put(debit[s],r.reactant_removed_energy_J_m3[s])||
     (numbers[s]==0&&energies[s]!=0))
    return PB11_STATUS_NUMERICAL_FAILURE;
   R nr=R(oldN[s])-numbers[s]-r.reactant_removed_m3[s];
   R er=R(oldU[s])-energies[s]-r.reactant_removed_energy_J_m3[s];
   if(!put(nr,r.number_residual_m3[s])||!put(er,r.energy_residual_J_m3[s])||
     std::abs(nr)>2e-12L*std::max(R(oldN[s]),loss[s])||
     std::abs(er)>2e-12L*std::max(R(oldU[s]),debit[s]))return PB11_STATUS_NUMERICAL_FAILURE;
  }
  std::copy(numbers.begin(),numbers.end(),newN);std::copy(energies.begin(),energies.end(),newU);
  *out=r;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
