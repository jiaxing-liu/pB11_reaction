#include "fusion_handoff.h"
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using R=long double;
bool put(R v,double& out,bool permit_underflow=false) {
 if(!std::isfinite(v)||std::abs(v)>std::numeric_limits<double>::max())return false;
 out=static_cast<double>(v);return permit_underflow||v==0||out!=0;
}
R cdf(R x){return std::isinf(x)?1:boost::math::gamma_p(1.5L,x);}
R survival(R x){return std::isinf(x)?0:boost::math::gamma_q(1.5L,x);}
int valid_grid(int n,double T,const double* e){
 if(n<1 || n==std::numeric_limits<int>::max())return PB11_STATUS_INVALID_ARGUMENT;
 if(!e || !std::isfinite(T))return PB11_STATUS_INVALID_ARGUMENT;
 if(T<=0)return PB11_STATUS_OUT_OF_RANGE;
 for(int i=0;i<=n;i++){
  if(!std::isfinite(e[i])||e[i]<0)return PB11_STATUS_INVALID_ARGUMENT;
  if(i&&e[i]<=e[i-1])return PB11_STATUS_OUT_OF_RANGE;
 }
 return PB11_STATUS_OK;
}
}
extern "C" int fusion_c_maxwellian_energy_grid(int n,double T,const double* e,
 double* q,fusion_maxwellian_grid_v1* out){
 if(out)*out={};
 if(n<1||n==std::numeric_limits<int>::max())return PB11_STATUS_INVALID_ARGUMENT;
 if(q)std::fill(q,q+n,0.);
 if(!q||!out)return PB11_STATUS_NULL_OUTPUT;
 int s=valid_grid(n,T,e);if(s)return s;
 try{
  std::vector<double> p(n);fusion_maxwellian_grid_v1 r{};
  const R lo=R(e[0])/T, hi=R(e[n])/T;
  if(!put(cdf(lo),r.below_probability,true)||!put(survival(hi),r.above_probability,true))
   return PB11_STATUS_NUMERICAL_FAILURE;
  R total=R(r.below_probability)+r.above_probability,mean=0;
  for(int i=0;i<n;i++){
   const R a=R(e[i])/T,b=R(e[i+1])/T;
   // Use the survival function in the upper tail to avoid subtracting two
   // CDFs indistinguishable from one. Nonnegative underflow is permitted.
   const R prob=a>=1.5L?survival(a)-survival(b):cdf(b)-cdf(a);
   if(prob<0||!put(prob,p[i],true))return PB11_STATUS_NUMERICAL_FAILURE;
   total+=p[i];mean+=R(p[i])*(R(e[i])+e[i+1])/2;
  }
  const R error=total-1;
  if(std::abs(error)>1e-12L || !put(mean,r.represented_mean_energy_J,true) ||
     !put(error,r.probability_balance_error,true))return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(p.begin(),p.end(),q);*out=r;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
extern "C" int fusion_c_maxwellian_handoff_trial(int n,double T,double maxL1,
 double maxMean,const double* e,const double* old,double* trial,int* projected,
 fusion_handoff_ledger_v1* out){
 if(out)*out={};
 if(projected)*projected=0;
 if(n<1||n==std::numeric_limits<int>::max())return PB11_STATUS_INVALID_ARGUMENT;
 if(trial)std::fill(trial,trial+n,0.);
 if(!out||!trial||!projected)return PB11_STATUS_NULL_OUTPUT;
 if(!old || !std::isfinite(maxL1)||!std::isfinite(maxMean))return PB11_STATUS_INVALID_ARGUMENT;
 if(maxL1<0||maxL1>2||maxMean<0||maxMean>1)return PB11_STATUS_OUT_OF_RANGE;
 int s=valid_grid(n,T,e);if(s)return s;
 try{
  R N=0,U=0;
  for(int i=0;i<n;i++){
   if(!std::isfinite(old[i])||old[i]<0)return PB11_STATUS_INVALID_ARGUMENT;
   N+=old[i];U+=R(old[i])*(R(e[i])+e[i+1])/2;
  }
  fusion_handoff_ledger_v1 r{};
  if(!put(N,r.initial_number_m3)||!put(U,r.initial_energy_J_m3))return PB11_STATUS_NUMERICAL_FAILURE;
  std::vector<double> q(n);fusion_maxwellian_grid_v1 grid{};
  s=fusion_c_maxwellian_energy_grid(n,T,e,q.data(),&grid);if(s)return s;
  const R outside=R(grid.below_probability)+grid.above_probability;
  if(!put(outside,r.outside_grid_probability,true))return PB11_STATUS_NUMERICAL_FAILURE;
  r.represented_Maxwellian_mean_energy_J=grid.represented_mean_energy_J;
  R l1=0,mean_error=0;
  if(N>0){
   l1=outside;for(int i=0;i<n;i++)l1+=std::abs(R(old[i])/N-q[i]);
   mean_error=std::abs(U/(1.5L*T*N)-1);
  }
  if(!put(l1,r.distribution_L1,true)||!put(mean_error,r.relative_mean_energy_error,true))
   return PB11_STATUS_NUMERICAL_FAILURE;
  const bool accept=N>0 && l1<=maxL1 && mean_error<=maxMean;
  if(accept){
   r.fluid_number_m3=r.initial_number_m3;
   if(!put(1.5L*T*r.fluid_number_m3,r.fluid_energy_J_m3)||
      !put(U-r.fluid_energy_J_m3,r.bath_energy_correction_J_m3,true))return PB11_STATUS_NUMERICAL_FAILURE;
   // Balance the quantities actually returned at double precision.
   const R nr=R(r.fluid_number_m3)-N;
   const R er=R(r.fluid_energy_J_m3)+r.bath_energy_correction_J_m3-U;
   if(std::abs(nr)>1e-12L*(N+r.fluid_number_m3)||
      std::abs(er)>1e-12L*(U+r.fluid_energy_J_m3+std::abs(R(r.bath_energy_correction_J_m3)))||
      !put(nr,r.particle_balance_error_m3,true)||!put(er,r.energy_balance_error_J_m3,true))
    return PB11_STATUS_NUMERICAL_FAILURE;
  }else{
   r.remaining_number_m3=r.initial_number_m3;r.remaining_energy_J_m3=r.initial_energy_J_m3;
  }
  if(!accept)std::copy(old,old+n,trial);
  *projected=accept?1:0;*out=r;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
