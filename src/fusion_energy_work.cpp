#include "fusion_energy_work.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using R=long double;
bool put(R x,double&v){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;v=double(x);return true;}
}
extern "C" int fusion_c_energy_work_trial(int n,double dt,double c,
 const double*edge,const double*old,double*out,double*face,fusion_energy_work_ledger_v1*ledger){
 if(n<1||n>1000000)return PB11_STATUS_INVALID_ARGUMENT;
 if(out)std::fill(out,out+n,0.);
 if(face)std::fill(face,face+n+1,0.);
 if(ledger)*ledger={};
 if(!out||!face||!ledger)return PB11_STATUS_NULL_OUTPUT;
 if(!edge||!old||!std::isfinite(dt)||dt<=0||!std::isfinite(c))return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<=n;++j)if(!std::isfinite(edge[j])||edge[j]<0||(j&&edge[j]<=edge[j-1]))return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<n;++j)if(!std::isfinite(old[j])||old[j]<0)return PB11_STATUS_INVALID_ARGUMENT;
 try{
  std::vector<R> width(n),energy(n),rate(n+1),value(n),round(n),flow_round(n+1);
  std::vector<double> result(n),flows(n+1);
  for(int j=0;j<n;++j){width[j]=R(edge[j+1])-edge[j];energy[j]=(R(edge[j+1])+edge[j])/2;}
  for(int j=0;j<=n;++j){
   const int donor=c>0?j:j-1;
   rate[j]=(donor<0||donor>=n)?0:std::abs(R(c))*edge[j]/width[donor];
  }
  for(int i=0;i<n;++i){
   const int j=c>0?n-1-i:i;
   const R outgoing=rate[c>0?j:j+1],incoming=rate[c>0?j+1:j];
   const int donor=c>0?j+1:j-1;
   const R rhs=R(old[j])+((donor<0||donor>=n)?0:R(dt)*incoming*value[donor]);
   value[j]=rhs/(1+R(dt)*outgoing);
   if(value[j]<0||!put(value[j],result[j]))return PB11_STATUS_NUMERICAL_FAILURE;
   round[j]=std::abs(value[j]-R(result[j]));
  }
  for(int j=0;j<=n;++j){
   const int donor=c>0?j:j-1;
   const R sign=c>0?-1:1;
   const R amount=(donor<0||donor>=n)?0:sign*R(dt)*rate[j]*result[donor];
   if(!put(amount,flows[j]))return PB11_STATUS_NUMERICAL_FAILURE;
   flow_round[j]=std::abs(amount-R(flows[j]))+((donor<0||donor>=n)?0:R(dt)*rate[j]*round[donor]);
  }
  R ni=0,nf=0,ui=0,uf=0,nround=0,eround=0;
  for(int j=0;j<n;++j){
   const R err=R(result[j])-old[j]+R(flows[j+1])-flows[j];
   const R scale=R(result[j])+old[j]+std::abs(R(flows[j]))+std::abs(R(flows[j+1]));
   const R allowance=round[j]+flow_round[j]+flow_round[j+1];
   if(!std::isfinite(err)||!std::isfinite(allowance)||std::abs(err)>1e-10L*scale+allowance)return PB11_STATUS_NUMERICAL_FAILURE;
   ni+=old[j];nf+=result[j];ui+=energy[j]*old[j];uf+=energy[j]*result[j];
   nround+=allowance;eround+=energy[j]*allowance;
  }
  const R low=-R(flows[0]),high=flows[n],elo=low*edge[0],ehi=high*edge[n];
  R work=(energy[0]-edge[0])*flows[0]+(R(edge[n])-energy[n-1])*flows[n];
  for(int j=1;j<n;++j)work+=(energy[j]-energy[j-1])*flows[j];
  const R ne=nf+low+high-ni,ee=uf+elo+ehi-ui-work;
  if(std::abs(ne)>1e-10L*(ni+nf+low+high)+nround||std::abs(ee)>1e-10L*(ui+uf+elo+ehi+std::abs(work))+eround)return PB11_STATUS_NUMERICAL_FAILURE;
  fusion_energy_work_ledger_v1 b{};
  if(!put(ni,b.initial_number_m3)||!put(nf,b.final_number_m3)||!put(ui,b.initial_energy_J_m3)||!put(uf,b.final_energy_J_m3)||!put(low,b.lower_number_m3)||!put(elo,b.lower_energy_J_m3)||!put(high,b.upper_number_m3)||!put(ehi,b.upper_energy_J_m3)||!put(work,b.work_on_particles_J_m3)||!put(ne,b.particle_balance_error_m3)||!put(ee,b.energy_balance_error_J_m3))return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(result.begin(),result.end(),out);std::copy(flows.begin(),flows.end(),face);*ledger=b;
  return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
