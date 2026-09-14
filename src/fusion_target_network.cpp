#include "fusion_target_network.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using R=long double;
struct Sum {R s=0,c=0;void add(R x){R t=s+x;if(std::abs(s)>=std::abs(x))c+=(s-t)+x;else c+=(x-t)+s;s=t;}R value()const{return s+c;}};
bool nn(double x){return std::isfinite(x)&&x>=0;}
bool put(R x,double& d){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;d=double(x);return true;}
bool fastdim(int n){return n>=1&&n<=600000;}
bool targetdim(int n){return n>=1&&n<=6;}
bool edgedim(int n){return n>=0&&n<=2000000;}
}
extern "C" int fusion_c_target_network_trial(int nf,int nt,int ne,double dt,
 const double*energy,const double*old,const double*B,const double*U,
 const fusion_target_network_edge_v1*edges,double*next,double*nextB,double*nextU,
 double*events,fusion_target_network_v1*out){
 if(out)*out={};
 if(fastdim(nf)&&next)std::fill(next,next+nf,0.);
 if(targetdim(nt)){if(nextB)std::fill(nextB,nextB+nt,0.);if(nextU)std::fill(nextU,nextU+nt,0.);}
 if(edgedim(ne)&&events)std::fill(events,events+ne,0.);
 if(!out||!next||!nextB||!nextU||(ne&&!events))return PB11_STATUS_NULL_OUTPUT;
 if(!fastdim(nf)||!targetdim(nt)||!edgedim(ne)||!energy||!old||!B||!U||(ne&&!edges)||!std::isfinite(dt)||dt<=0)return PB11_STATUS_INVALID_ARGUMENT;
 for(int i=0;i<nf;++i)if(!nn(energy[i])||!nn(old[i]))return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<nt;++j)if(!nn(B[j])||!nn(U[j])||(B[j]==0&&U[j]!=0))return PB11_STATUS_INVALID_ARGUMENT;
 for(int e=0;e<ne;++e){const auto&v=edges[e];if(v.fast_index<0||v.fast_index>=nf||v.target_index<0||v.target_index>=nt||!nn(v.reactivity_m3_s)||!nn(v.target_energy_reactivity_J_m3_s)||(v.reactivity_m3_s==0&&v.target_energy_reactivity_J_m3_s!=0))return PB11_STATUS_INVALID_ARGUMENT;}
 try{
  std::vector<R>a(size_t(nf)*nt,0),f(nf),den(nf);
  for(int e=0;e<ne;++e)a[size_t(edges[e].fast_index)*nt+edges[e].target_index]+=R(dt)*edges[e].reactivity_m3_s;
  std::array<R,6>b{},g{},step{};std::array<int,6>active{};int na=0;
  for(int j=0;j<nt;++j){R stiff=0;for(int i=0;i<nf;++i)stiff+=a[size_t(i)*nt+j]*old[i];b[j]=R(B[j])/(1+stiff);if(B[j]>0)active[na++]=j;}
  auto evaluate=[&](const std::array<R,6>&x,bool derivatives,std::array<R,6>&gradient,R h[6][6]){
   std::array<Sum,6> sums{};
   for(int j=0;j<nt;++j){sums[j].add(x[j]);sums[j].add(-R(B[j]));if(derivatives){for(int k=0;k<nt;++k)h[j][k]=0;h[j][j]=x[j];}}
   for(int i=0;i<nf;++i){R d=1;for(int j=0;j<nt;++j)d+=a[size_t(i)*nt+j]*x[j];den[i]=d;f[i]=R(old[i])/d;
    std::array<R,6>p{};for(int j=0;j<nt;++j){p[j]=a[size_t(i)*nt+j]*x[j]/d;if(p[j]>.5L){sums[j].add(old[i]);R complement=1;for(int k=0;k<nt;++k)if(k!=j)complement+=a[size_t(i)*nt+k]*x[k];sums[j].add(-f[i]*complement);}else sums[j].add(f[i]*a[size_t(i)*nt+j]*x[j]);}
    if(derivatives)for(int j=0;j<nt;++j){R other=1/d;for(int k=0;k<nt;++k)if(k!=j)other+=p[k];h[j][j]+=R(old[i])*p[j]*other;for(int k=0;k<j;++k){R v=R(old[i])*p[j]*p[k];h[j][k]-=v;h[k][j]-=v;}}
   }
   for(int j=0;j<nt;++j)gradient[j]=sums[j].value();
   R residual=0;for(int z=0;z<na;++z){int j=active[z];residual=std::max(residual,std::abs(gradient[j])/B[j]);}return residual;
  };
  bool converged=false;int iteration=0;
  for(;iteration<512;++iteration){
   R h[6][6]{};R norm=evaluate(b,true,g,h);if(!std::isfinite(norm))return PB11_STATUS_NUMERICAL_FAILURE;
   R local=0;for(int z=0;z<na;++z){int j=active[z];local=std::max(local,std::abs(g[j])/h[j][j]);}
   if(norm<=8*std::numeric_limits<double>::epsilon()&&local<=2e-14L){converged=true;break;}
   // Hessian of sum B'-B0*log(B') + sum F0*log(1+sum dt*K*B').
   // Row scaling by initial target inventory keeps heterogeneous baths usable.
   R mat[6][7]{};for(int z=0;z<na;++z){int j=active[z];for(int q=0;q<na;++q)mat[z][q]=h[j][active[q]]/B[j];mat[z][na]=-g[j]/B[j];}
   bool solved=true;
   for(int k=0;k<na;++k){int pivot=k;for(int z=k+1;z<na;++z)if(std::abs(mat[z][k])>std::abs(mat[pivot][k]))pivot=z;
    if(mat[pivot][k]==0||!std::isfinite(mat[pivot][k])){solved=false;break;}
    for(int q=k;q<=na;++q)std::swap(mat[k][q],mat[pivot][q]);
    for(int z=k+1;z<na;++z){R factor=mat[z][k]/mat[k][k];for(int q=k;q<=na;++q)mat[z][q]-=factor*mat[k][q];}
   }
   step.fill(0);if(solved)for(int z=na-1;z>=0;--z){R v=mat[z][na];for(int q=z+1;q<na;++q)v-=mat[z][q]*step[active[q]];step[active[z]]=v/mat[z][z];if(!std::isfinite(step[active[z]]))solved=false;}
   R descent=0,maxstep=0;if(solved)for(int z=0;z<na;++z){int j=active[z];descent+=g[j]*step[j];maxstep=std::max(maxstep,std::abs(step[j]));}
   if(!solved||descent>=0){descent=0;maxstep=0;for(int z=0;z<na;++z){int j=active[z];step[j]=-g[j]/h[j][j];descent+=g[j]*step[j];maxstep=std::max(maxstep,std::abs(step[j]));}}
   if(!std::isfinite(maxstep)||maxstep==0)return PB11_STATUS_NUMERICAL_FAILURE;
   R alpha=std::min(R(1),8/maxstep);bool accepted=false;
   for(int search=0;search<80;++search){std::array<R,6>candidate=b,delta{},vstep{};Sum objective;
    for(int z=0;z<na;++z){int j=active[z];R v=alpha*step[j];vstep[j]=v;delta[j]=std::expm1(v);candidate[j]=b[j]*std::exp(v);objective.add(b[j]*delta[j]);objective.add(-R(B[j])*v);}
    for(int i=0;i<nf;++i){int dominant=-1;for(int j=0;j<nt;++j)if(a[size_t(i)*nt+j]*b[j]/den[i]>.5L)dominant=j;
     R relative=0;if(dominant>=0){R shift=vstep[dominant];objective.add(R(old[i])*shift);relative=std::expm1(-shift)/den[i];for(int j=0;j<nt;++j)if(j!=dominant)relative+=a[size_t(i)*nt+j]*b[j]/den[i]*std::expm1(vstep[j]-shift);}
     else {for(int j=0;j<nt;++j)relative+=a[size_t(i)*nt+j]*b[j]/den[i]*delta[j];}
     objective.add(R(old[i])*std::log1p(relative));}
    R change=objective.value();
    if(std::isfinite(change)&&change<=1e-4L*alpha*descent){b=candidate;accepted=true;break;}alpha*=.5L;
   }
   if(!accepted)return PB11_STATUS_NUMERICAL_FAILURE;
  }
  if(!converged)return PB11_STATUS_NUMERICAL_FAILURE;
  std::vector<double>newf(nf),loss(ne);std::vector<R>removed(nf,0);std::array<R,6>targetloss{},targetdebit{};
  R total=0,fastE=0,targetE=0;
  for(int i=0;i<nf;++i)if(!put(f[i],newf[i]))return PB11_STATUS_NUMERICAL_FAILURE;
  for(int e=0;e<ne;++e){const auto&v=edges[e];R amount=R(dt)*v.reactivity_m3_s*f[v.fast_index]*b[v.target_index];if(!put(amount,loss[e]))return PB11_STATUS_NUMERICAL_FAILURE;
   // All ledger sums are based on the actual returned event double.
   amount=loss[e];R debit=v.reactivity_m3_s>0?amount*R(v.target_energy_reactivity_J_m3_s)/v.reactivity_m3_s:0;
   removed[v.fast_index]+=amount;targetloss[v.target_index]+=amount;targetdebit[v.target_index]+=debit;total+=amount;fastE+=amount*energy[v.fast_index];targetE+=debit;
  }
  fusion_target_network_v1 result{};std::array<double,6>bn{},un{};
  for(int j=0;j<nt;++j){R remaining=R(U[j])-targetdebit[j];if(remaining<0||!put(b[j],bn[j])||!put(remaining,un[j]))return PB11_STATUS_NUMERICAL_FAILURE;
   R scale=std::max(R(B[j]),targetloss[j]);R r=scale>0?std::abs(R(B[j])-bn[j]-targetloss[j])/scale:0;
   result.max_target_number_relative_residual=std::max(result.max_target_number_relative_residual,double(r));
   scale=std::max(R(U[j]),targetdebit[j]);r=scale>0?std::abs(R(U[j])-un[j]-targetdebit[j])/scale:0;
   result.max_target_energy_relative_residual=std::max(result.max_target_energy_relative_residual,double(r));
  }
  for(int i=0;i<nf;++i){R scale=std::max(R(old[i]),removed[i]);R r=scale>0?std::abs(R(old[i])-newf[i]-removed[i])/scale:0;result.max_fast_number_relative_residual=std::max(result.max_fast_number_relative_residual,double(r));}
  if(result.max_fast_number_relative_residual>2e-12||result.max_target_number_relative_residual>2e-12||result.max_target_energy_relative_residual>2e-12||!put(total,result.reactions_m3)||!put(fastE,result.removed_fast_energy_J_m3)||!put(targetE,result.removed_target_energy_J_m3))return PB11_STATUS_NUMERICAL_FAILURE;
  result.iterations=iteration;std::copy(newf.begin(),newf.end(),next);std::copy(bn.begin(),bn.begin()+nt,nextB);std::copy(un.begin(),un.begin()+nt,nextU);if(ne)std::copy(loss.begin(),loss.end(),events);*out=result;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
