#include "fusion_radial_transport.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using R=long double;
bool nonneg(double x){return std::isfinite(x)&&x>=0;}
bool put(R x,double&v){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;v=double(x);return true;}
}
extern "C" int fusion_c_radial_transport_trial(int z,int c,double dt,
 const double*vo,const double*vn,const double*a,const double*k,const double*old,
 const double*boundary,double*trial,double*face,fusion_radial_ledger_v1*ledger){
 const long long entries=static_cast<long long>(c)*(2LL*z+1);
 if(z<1||c<1||entries>50000000)return PB11_STATUS_INVALID_ARGUMENT;
 const auto nz=static_cast<std::size_t>(z),nc=static_cast<std::size_t>(c),nf=nz+1;
 if(trial)std::fill(trial,trial+nz*nc,0.);
 if(face)std::fill(face,face+nf*nc,0.);
 if(ledger)std::fill(ledger,ledger+nc,fusion_radial_ledger_v1{});
 if(!trial||!face||!ledger)return PB11_STATUS_NULL_OUTPUT;
 if(!vo||!vn||!a||!k||!old||!boundary||!std::isfinite(dt)||dt<=0)return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<z;++j)if(!std::isfinite(vo[j])||vo[j]<=0||!std::isfinite(vn[j])||vn[j]<=0)return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<=z;++j)if(!std::isfinite(a[j]))return PB11_STATUS_INVALID_ARGUMENT;
 for(std::size_t j=0;j<nz*nc;++j)if(!nonneg(old[j]))return PB11_STATUS_INVALID_ARGUMENT;
 for(std::size_t j=0;j<nf*nc;++j)if(!nonneg(k[j]))return PB11_STATUS_INVALID_ARGUMENT;
 for(std::size_t j=0;j<2*nc;++j)if(!nonneg(boundary[j]))return PB11_STATUS_INVALID_ARGUMENT;
 try{
  std::vector<double> result(nz*nc),flows(nf*nc);
  std::vector<fusion_radial_ledger_v1> balances(nc);
  std::vector<R> left(nf),right(nf),diag(nz),lower(nz),upper(nz),rhs(nz),n(nz),round(nz),fr(nf);
  for(std::size_t p=0;p<nc;++p){
   const auto off=p*nz,foff=p*nf;
   for(std::size_t j=0;j<nf;++j){left[j]=std::max(R(a[j]),R(0))+k[foff+j];right[j]=std::min(R(a[j]),R(0))-k[foff+j];}
   for(std::size_t j=0;j<nz;++j){
    diag[j]=R(vn[j])+R(dt)*(left[j+1]-right[j]);
    lower[j]=j?-R(dt)*left[j]:0;upper[j]=j+1<nz?R(dt)*right[j+1]:0;
    rhs[j]=R(vo[j])*old[off+j];
   }
   rhs[0]+=R(dt)*left[0]*boundary[2*p];rhs[nz-1]-=R(dt)*right[nz]*boundary[2*p+1];
   for(std::size_t j=1;j<nz;++j){
    if(!std::isfinite(diag[j-1])||diag[j-1]<=0)return PB11_STATUS_NUMERICAL_FAILURE;
    R factor=lower[j]/diag[j-1];diag[j]-=factor*upper[j-1];rhs[j]-=factor*rhs[j-1];
   }
   for(std::size_t jj=nz;jj>0;--jj){const auto j=jj-1;
    if(!std::isfinite(diag[j])||diag[j]<=0)return PB11_STATUS_NUMERICAL_FAILURE;
    n[j]=(rhs[j]-(j+1<nz?upper[j]*n[j+1]:0))/diag[j];
    if(!std::isfinite(n[j])||n[j]<0||!put(n[j],result[off+j]))return PB11_STATUS_NUMERICAL_FAILURE;
    round[j]=std::abs(n[j]-R(result[off+j]));
   }
   for(std::size_t j=0;j<nf;++j){
    R nl=j?result[off+j-1]:boundary[2*p],nr=j<nz?result[off+j]:boundary[2*p+1];
    R amount=R(dt)*(left[j]*nl+right[j]*nr);
    if(!put(amount,flows[foff+j]))return PB11_STATUS_NUMERICAL_FAILURE;
    fr[j]=std::abs(amount-R(flows[foff+j]))+R(dt)*(left[j]*(j?round[j-1]:0)-right[j]*(j<nz?round[j]:0));
   }
   R initial=0,final=0,round_total=0;
   for(std::size_t j=0;j<nz;++j){
    R ni=R(vo[j])*old[off+j],nn=R(vn[j])*result[off+j];
    R err=nn-ni+R(flows[foff+j+1])-flows[foff+j];
    R scale=ni+nn+std::abs(R(flows[foff+j]))+std::abs(R(flows[foff+j+1]));
    R allowance=R(vn[j])*round[j]+fr[j]+fr[j+1];
    if(!std::isfinite(err)||!std::isfinite(allowance)||std::abs(err)>1e-10L*scale+allowance)return PB11_STATUS_NUMERICAL_FAILURE;
    initial+=ni;final+=nn;round_total+=allowance;
   }
   R err=final-initial-R(flows[foff])+flows[foff+nz];
   R scale=initial+final+std::abs(R(flows[foff]))+std::abs(R(flows[foff+nz]));
   if(!std::isfinite(err)||std::abs(err)>1e-10L*scale+round_total)return PB11_STATUS_NUMERICAL_FAILURE;
   auto&b=balances[p];
   if(!put(initial,b.initial_number)||!put(final,b.final_number)||!put(err,b.balance_error))return PB11_STATUS_NUMERICAL_FAILURE;
   b.inner_inward_number=flows[foff];b.outer_outward_number=flows[foff+nz];
  }
  std::copy(result.begin(),result.end(),trial);std::copy(flows.begin(),flows.end(),face);std::copy(balances.begin(),balances.end(),ledger);
  return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
