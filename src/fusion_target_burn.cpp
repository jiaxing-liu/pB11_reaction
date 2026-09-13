#include "fusion_target_burn.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using L=long double;
bool nonnegative(double x){return std::isfinite(x) && x>=0;}
bool representable(L x){return std::isfinite(x) && std::abs(x)<=std::numeric_limits<double>::max();}
void clear(int n,double*trial,double*loss,fusion_target_burn_v1*out){
 if(out)*out={};
 if(n>0){if(trial)std::fill(trial,trial+n,0.);if(loss)std::fill(loss,loss+n,0.);}
}
bool grid_ok(int n,const double*e,const double*N){
 if(n<1 || !e || !N)return false;
 for(int i=0;i<=n;++i)if(!nonnegative(e[i]) || (i>0 && e[i]<=e[i-1]))return false;
 for(int i=0;i<n;++i)if(!nonnegative(N[i]))return false;
 return true;
}
}
extern "C" int fusion_c_target_burn_trial(int n,double dt,const double*edges,
 const double*old,double B,double U,const double*K,const double*M,
 double*trial,double*loss,fusion_target_burn_v1*out){
 clear(n,trial,loss,out);
 if(!trial || !loss || !out)return PB11_STATUS_NULL_OUTPUT;
 if(!grid_ok(n,edges,old) || !K || !M || !std::isfinite(dt) || dt<=0 ||
  !nonnegative(B) || !nonnegative(U) || (B==0 && U!=0))return PB11_STATUS_INVALID_ARGUMENT;
 for(int i=0;i<n;++i)if(!nonnegative(K[i]) || !nonnegative(M[i]) || (K[i]==0 && M[i]!=0))return PB11_STATUS_INVALID_ARGUMENT;
 try{
  std::vector<L> a(n);L stiffness=0,fast_initial=0;
  for(int i=0;i<n;++i){a[i]=L(dt)*K[i]*B;stiffness+=L(dt)*K[i]*old[i];fast_initial+=old[i];}
  L x=1;
  if(B>0 && stiffness>0){
   // A positive lower bound prevents loss of a very small surviving target
   // population. Safeguarded Newton uses geometric bisection across scales.
   L lo=std::max(1/(1+stiffness),std::max(L(0),1-fast_initial/B)),hi=1;
   x=lo;bool converged=false;
   for(int iteration=0;iteration<256;++iteration){
    L f=x-1,df=1;
    for(int i=0;i<n;++i){L z=a[i]*x,d=1+z,ratio=L(old[i])/B;f+=ratio*z/d;df+=ratio*a[i]/(d*d);}
    if(std::abs(f)<2e-16L){converged=true;break;}
    if(f>0)hi=x;else lo=x;
    L next=x-f/df;
    if(!(next>lo && next<hi) || !std::isfinite(next))next=std::sqrt(lo)*std::sqrt(hi);
    if(next==x){converged=std::abs(f)<2e-13L;break;}
    x=next;
   }
   if(!converged)return PB11_STATUS_NUMERICAL_FAILURE;
  }
  std::vector<double> numbers(n),removed(n);
  L fast_final=0,E0=0,E1=0,R=0,removedE=0,removedT=0;
  for(int i=0;i<n;++i){
   L z=a[i]*x,Ni=L(old[i])/(1+z),Ri=L(old[i])*z/(1+z),E=(L(edges[i])+edges[i+1])/2;
   L debit=K[i]>0?Ri*L(M[i])/K[i]:0;
   if(!representable(Ni) || !representable(Ri) || !representable(debit))return PB11_STATUS_NUMERICAL_FAILURE;
   numbers[i]=double(Ni);removed[i]=double(Ri);
   fast_final+=numbers[i];E0+=L(old[i])*E;E1+=L(numbers[i])*E;
   R+=removed[i];removedE+=L(removed[i])*E;removedT+=L(removed[i])*(K[i]>0?L(M[i])/K[i]:0);
  }
  L B1=L(B)*x,U1=L(U)-removedT;
  if(U1<0 || !representable(B1) || !representable(U1))return PB11_STATUS_NUMERICAL_FAILURE;
  // Residuals use the actual double-precision state returned to the caller.
  B1=double(B1);U1=double(U1);
  L values[]={fast_initial,fast_final,E0,E1,B,B1,U,U1,R,removedE,removedT,
   fast_initial-fast_final-R,L(B)-B1-R,E0-E1-removedE+L(U)-U1-removedT};
  for(L v:values)if(!representable(v))return PB11_STATUS_NUMERICAL_FAILURE;
  fusion_target_burn_v1 result{double(values[0]),double(values[1]),double(values[2]),double(values[3]),
   double(values[4]),double(values[5]),double(values[6]),double(values[7]),double(values[8]),
   double(values[9]),double(values[10]),double(values[11]),double(values[12]),double(values[13])};
  if(std::abs(values[11])>2e-12L*std::max(fast_initial,R) ||
     std::abs(values[12])>2e-12L*std::max(L(B),R) ||
     std::abs(values[13])>2e-12L*std::max(E0+U,removedE+removedT))return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(numbers.begin(),numbers.end(),trial);std::copy(removed.begin(),removed.end(),loss);*out=result;
  return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
extern "C" int fusion_c_beam_target_burn_window_trial(int channel,int n,double dt,
 double ma,double mb,const double*edges,const double*old,double B,double T,
 double*trial,double*loss,fusion_beam_window_v1*windows,fusion_target_burn_v1*out){
 clear(n,trial,loss,out);
 if(n>0 && windows)std::fill(windows,windows+n,fusion_beam_window_v1{});
 if(!trial || !loss || !windows || !out)return PB11_STATUS_NULL_OUTPUT;
 if(!grid_ok(n,edges,old) || !nonnegative(B) || !nonnegative(T) ||
  !std::isfinite(dt) || dt<=0)return PB11_STATUS_INVALID_ARGUMENT;
 L U=1.5L*B*T;if(!representable(U))return PB11_STATUS_NUMERICAL_FAILURE;
 try{
  std::vector<double>K(n),M(n);std::vector<fusion_beam_window_v1>w(n);
  for(int i=0;i<n;++i){
   double E=double((L(edges[i])+edges[i+1])/2);
   int status=fusion_c_beam_maxwellian_window(channel,ma,mb,E,T,&w[i]);if(status)return status;
   K[i]=w[i].resolved_reactivity_m3_s;M[i]=w[i].target_energy_reactivity_J_m3_s;
  }
  int status=fusion_c_target_burn_trial(n,dt,edges,old,B,double(U),K.data(),M.data(),trial,loss,out);
  if(status)return status;
  std::copy(w.begin(),w.end(),windows);return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
