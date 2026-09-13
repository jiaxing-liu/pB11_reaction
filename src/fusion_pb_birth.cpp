#include "fusion_pb_birth.h"
#include "fusion_alpha_spectrum.h"
#include "fusion_nuclear_data.h"
#include "fusion_products.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using R=long double;
constexpr double mev=1.602176634e-13;
bool finite_value(double x){return std::isfinite(x);}
int validate(double A,double q,int n,const double* edges,double* birth,fusion_pb_birth_v1* out){
 if(out)*out={};
 if(n>0 && birth)std::fill(birth,birth+n,0.);
 if(!birth || !out)return PB11_STATUS_NULL_OUTPUT;
 if(n<1 || !edges || !finite_value(A) || !finite_value(q))return PB11_STATUS_INVALID_ARGUMENT;
 if(A<=0 || A>12*mev || q<0 || q>A)return PB11_STATUS_OUT_OF_RANGE;
 for(int i=0;i<=n;++i){
  if(!finite_value(edges[i]) || edges[i]<0)return PB11_STATUS_INVALID_ARGUMENT;
  if(i && edges[i]<=edges[i-1])return PB11_STATUS_OUT_OF_RANGE;
 }
 return PB11_STATUS_OK;
}
struct Accum {
 std::vector<R> c,b;
 R bn=0,be=0,an=0,ae=0;
 explicit Accum(int n,const double* e):c(n),b(n){for(int i=0;i<n;++i)c[i]=(R(e[i])+e[i+1])/2;}
 void delta(R e,R w){
  if(e<c.front()){bn+=w;be+=w*e;return;}
  if(e>c.back()){an+=w;ae+=w*e;return;}
  auto it=std::lower_bound(c.begin(),c.end(),e);auto j=it-c.begin();
  if(*it==e)b[j]+=w;
  else {R f=(e-c[j-1])/(c[j]-c[j-1]);b[j]+=w*f;b[j-1]+=w*(1-f);}
 }
 void box(R lo,R hi){
  if(hi==lo){delta(lo,2);return;}
  const R width=hi-lo;
  auto part=[&](R l,R u,R& count,R& energy){
   l=std::max(l,lo);u=std::min(u,hi);if(u<=l)return;
   R w=2*(u-l)/width;count+=w;energy+=w*(l+(u-l)/2);
  };
  part(lo,c.front(),bn,be);part(c.back(),hi,an,ae);
  for(std::size_t j=1;j<c.size();++j){
   R l=std::max(lo,c[j-1]),u=std::min(hi,c[j]);if(u<=l)continue;
   R w=2*(u-l)/width,mid=l+(u-l)/2,f=(mid-c[j-1])/(c[j]-c[j-1]);
   b[j]+=w*f;b[j-1]+=w*(1-f);
  }
 }
 int finish(double A,double* birth,fusion_pb_birth_v1& result){
  std::vector<double> v(b.size());R mn=0,me=0;
  for(std::size_t i=0;i<b.size();++i){
   if(!std::isfinite(b[i]) || b[i]<0 || b[i]>3.00000000001L)return PB11_STATUS_NUMERICAL_FAILURE;
   v[i]=double(b[i]);mn+=v[i];me+=c[i]*v[i];
  }
  result.mapped_number=double(mn);result.mapped_energy_J=double(me);
  result.below_number=double(bn);result.below_energy_J=double(be);
  result.above_number=double(an);result.above_energy_J=double(ae);
  R nr=3-mn-R(result.below_number)-result.above_number;
  R er=R(A)-me-R(result.below_energy_J)-result.above_energy_J;
  result.number_residual=double(nr);result.energy_residual_J=double(er);
  if(!std::isfinite(nr) || !std::isfinite(er) || std::abs(nr)>3e-12L || std::abs(er)>1e-12L*A)return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(v.begin(),v.end(),birth);return PB11_STATUS_OK;
 }
};
}
extern "C" int fusion_c_pb_alpha0_grid(double A,double q,int n,const double* edges,
 double* birth,fusion_pb_birth_v1* out){
 int status=validate(A,q,n,edges,birth,out);if(status)return status;
 try {
  fusion_nuclear_mass_v1 mass{};status=fusion_c_nuclear_mass(FUSION_HELIUM4,&mass);if(status)return status;
  fusion_three_body_cm_v1 e{};status=fusion_c_three_equal_sequential_cm(mass.mass_kg,A,q,1,&e);if(status)return status;
  fusion_pb_birth_v1 result{};result.alpha0_fraction=1;
  result.primary_alpha0_energy_J=e.kinetic_energy_J[0];
  result.secondary_min_J=std::min(e.kinetic_energy_J[1],e.kinetic_energy_J[2]);
  result.secondary_max_J=std::max(e.kinetic_energy_J[1],e.kinetic_energy_J[2]);
  Accum acc(n,edges);acc.delta(result.primary_alpha0_energy_J,1);acc.box(result.secondary_min_J,result.secondary_max_J);
  status=acc.finish(A,birth,result);if(status)return status;
  *out=result;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
extern "C" int fusion_c_pb_cm_source_grid(double A,double q,double f0,double flow,
 int mode,int policy,double cutoff,double k,double phase,int nq,int nc,int n,
 const double* edges,double* birth,fusion_pb_birth_v1* out){
 int status=validate(A,q,n,edges,birth,out);if(status)return status;
 if(!finite_value(f0)||!finite_value(flow)||!finite_value(cutoff)||!finite_value(k)||!finite_value(phase) ||
    (mode!=1 && mode!=3 && mode!=13)||(policy!=0 && policy!=1)||nq<4||nq>1024||nc<4||nc>1024)
  return PB11_STATUS_INVALID_ARGUMENT;
 if(f0<0||flow<0||f0+flow>1||cutoff<.001*mev||cutoff>.01*mev||k<0||k>1)return PB11_STATUS_OUT_OF_RANGE;
 try {
  const double fb=1-(f0+flow);std::vector<double> part(n);fusion_pb_birth_v1 result{};
  status=fusion_c_pb_alpha0_grid(A,q,n,edges,part.data(),&result);if(status)return status;
  Accum acc(n,edges);for(int i=0;i<n;++i)acc.b[i]=R(f0)*part[i];
  acc.bn=R(f0)*result.below_number;acc.be=R(f0)*result.below_energy_J;
  acc.an=R(f0)*result.above_number;acc.ae=R(f0)*result.above_energy_J;
  const double weights[2]={flow,fb};const int modes[2]={2,mode};
  for(int j=0;j<2;++j)if(weights[j]>0){
   fusion_alpha_spectrum_v1 s{};
   status=fusion_c_alpha_spectrum_model_grid(modes[j],policy,A,cutoff,k,phase,nq,nc,n,edges,part.data(),&s);if(status)return status;
   const R w=weights[j];for(int i=0;i<n;++i)acc.b[i]+=w*part[i];
   acc.bn+=w*s.below_number;acc.be+=w*s.below_energy_J;acc.an+=w*s.above_number;acc.ae+=w*s.above_energy_J;
  }
  result.alpha0_fraction=f0;result.low_alpha1_fraction=flow;result.broad_alpha1_fraction=fb;
  status=acc.finish(A,birth,result);if(status)return status;
  *out=result;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
