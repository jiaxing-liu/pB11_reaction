#include "fusion_thermal_birth_budget.h"
#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
namespace {
constexpr double kev=1.602176634e-16,mev=1.602176634e-13;
constexpr int cells=8;
using Birth=std::array<double,7*cells>;
const std::array<double,9> edges{{0,.1*mev,.2*mev,.4*mev,.8*mev,1.2*mev,1.6*mev,2*mev,2.4*mev}};
void need(bool b,const char*m){if(!b)throw std::runtime_error(m);}
template<class T> bool zero(const T&v){const unsigned char*p=reinterpret_cast<const unsigned char*>(&v);return std::all_of(p,p+sizeof(v),[](unsigned char x){return x==0;});}
fusion_thermal_birth_options_v1 options(int remainder=0){
 fusion_thermal_birth_options_v1 o{};o.relative_max_J=2*mev;o.cm_max_kT=40;o.ground_state_q_J=91.84*kev;o.cutoff_J=.001*mev;o.l1_fraction=.76;o.narrow_peak_fraction=.051;o.continuum_peak_scale=1;o.continuation=FUSION_ENDPOINT_S;o.pb_low=FUSION_PB_LOW_TB;o.remainder_policy=remainder;o.broad_mode=13;o.relative_order=4;o.cm_order=4;o.nq=8;o.ncos=8;return o;
}
fusion_thermal_birth_budget_policy_v1 policy(double n,double e){return {1,8,4,0,n,e};}
struct Result{Birth b{};fusion_thermal_birth_v1 out{};fusion_thermal_birth_budget_report_v1 report{};};
Result evaluate(double T,const fusion_thermal_birth_options_v1&o,const fusion_thermal_birth_budget_policy_v1&p){Result r;need(fusion_c_pb_thermal_birth_budget_grid(T,&o,&p,cells,edges.data(),r.b.data(),&r.out,&r.report)==0,"budget evaluation failed");return r;}
void counts(const Result&r,const fusion_thermal_birth_budget_policy_v1&p){
 auto&q=r.report;need(q.version==1&&q.full_nq==8&&q.full_ncos==8&&q.coarse_nq==p.coarse_nq&&q.coarse_ncos==p.coarse_ncos,"report quadrature metadata");
 need(q.high_pairs+q.low_pairs==q.joint_nodes&&q.joint_nodes==q.relative_nodes*q.cm_nodes,"all positive parents retained");
 need(q.mixed==(q.low_pairs>0),"mixed flag tracks low parents");
 need(q.high_relative_nodes<=q.relative_nodes&&q.low_relative_nodes<=q.relative_nodes,"group accounting");
 need(q.number_L1_bound>=0&&q.energy_L1_bound>=0&&q.number_L1_bound<=p.number_L1_tolerance&&q.energy_L1_bound<=p.energy_L1_tolerance,"separate budgets");
 need(q.number_roundoff_allowance==1e-12&&q.energy_roundoff_allowance==1e-12,"declared engineering allowance");
 if(p.number_L1_tolerance==0||p.energy_L1_tolerance==0)need(q.low_pairs==0&&q.high_pairs==q.joint_nodes&&q.number_L1_bound==0&&q.energy_L1_bound==0,"either zero forces all high");
}
void compare(const Result&full,const Result&mix,const fusion_thermal_birth_budget_policy_v1&p){
 counts(mix,p);auto&a=full.out;auto&b=mix.out;
 need(a.reactivity_m3_s==b.reactivity_m3_s&&a.reference_reactivity_m3_s==b.reference_reactivity_m3_s,"rate unchanged");
 for(int j=0;j<2;++j)need(a.reactant_energy_moment_J_m3_s[j]==b.reactant_energy_moment_J_m3_s[j]&&a.reference_reactant_energy_moment_J_m3_s[j]==b.reference_reactant_energy_moment_J_m3_s[j],"debits unchanged");
 long double N=0,U=0,dN=0,dU=0,mixN=0,mixU=0;
 for(int s=0;s<7;++s){
  for(int i=0;i<cells;++i){int k=s*cells+i;need(std::isfinite(mix.b[k])&&mix.b[k]>=0,"positive finite mapped source");long double center=(static_cast<long double>(edges[i])+edges[i+1])/2,delta=std::abs(static_cast<long double>(full.b[k])-mix.b[k]);N+=full.b[k];U+=full.b[k]*center;mixN+=mix.b[k];mixU+=mix.b[k]*center;dN+=delta;dU+=delta*center;}
  for(int side=0;side<2;++side){double an=side?a.above_number_m3_s[s]:a.below_number_m3_s[s],bn=side?b.above_number_m3_s[s]:b.below_number_m3_s[s],au=side?a.above_energy_J_m3_s[s]:a.below_energy_J_m3_s[s],bu=side?b.above_energy_J_m3_s[s]:b.below_energy_J_m3_s[s];need(std::isfinite(bn)&&std::isfinite(bu)&&bn>=0&&bu>=0,"positive finite spills");N+=an;U+=au;mixN+=bn;mixU+=bu;dN+=std::abs(static_cast<long double>(an)-bn);dU+=std::abs(static_cast<long double>(au)-bu);}
 }
 need(N>0&&U>0&&a.above_number_m3_s[4]>0,"nonzero reference with real above-grid tail");
 fusion_nuclear_channel_v1 reaction{};need(fusion_c_nuclear_channel(0,&reaction)==0,"public reaction ledger");
 long double expectedN=3*static_cast<long double>(a.reactivity_m3_s),expectedU=static_cast<long double>(a.reactivity_m3_s)*reaction.q_J+a.reactant_energy_moment_J_m3_s[0]+a.reactant_energy_moment_J_m3_s[1];
 need(std::abs(N/expectedN-1)<1e-12L&&std::abs(U/expectedU-1)<1e-12L,"independent stoichiometry and Q plus debit closure");
 need(std::abs(mix.report.reaction_weight_m3_s/a.reactivity_m3_s-1)<1e-12&&std::abs(mix.report.product_energy_weight_J_m3_s/static_cast<double>(expectedU)-1)<1e-12,"report parent weights match physical ledger");
 need(dN/N<=mix.report.number_L1_bound+1e-12L&&dU/U<=mix.report.energy_L1_bound+1e-12L,"cells plus spills obey declared mixture bound");
 need(std::abs(mixN/N-1)<1e-12L&&std::abs(mixU/U-1)<1e-12L,"number and energy preserved");
}
void real_cases(){
 for(auto setup: {std::pair<double,int>{.5*kev,0},{10*kev,1},{10*kev,2}}){
  auto o=options(setup.second);auto p=policy(0,0);auto full=evaluate(setup.first,o,p);Birth old{};fusion_thermal_birth_v1 out{};
  need(fusion_c_thermal_birth_grid(0,setup.first,&o,cells,edges.data(),old.data(),&out)==0,"old API regression");
  need(std::memcmp(old.data(),full.b.data(),sizeof(old))==0&&std::memcmp(&out,&full.out,sizeof(out))==0,"zero budget whole output byte identity");counts(full,p);
  for(auto budget:{policy(.01,.001),policy(.001,.01),policy(2,2)}){auto mixed=evaluate(setup.first,o,budget);compare(full,mixed,budget);if(budget.number_L1_tolerance==2)need(mixed.report.high_pairs==0&&mixed.report.low_pairs>0,"maximal budget takes all low");}
  for(auto budget:{policy(0,.01),policy(.01,0)}){auto same=evaluate(setup.first,o,budget);counts(same,budget);need(std::memcmp(same.b.data(),full.b.data(),sizeof(Birth))==0&&std::memcmp(&same.out,&full.out,sizeof(out))==0,"one zero budget byte identity");}
 }
}
void rejection(){
 auto base=options();auto good=policy(.01,.01);double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
 auto reject=[&](double T,const fusion_thermal_birth_options_v1*o,const fusion_thermal_birth_budget_policy_v1*p,int n,const double*e,bool nullbirth=false,bool nullout=false,bool nullreport=false){
  Result r;std::memset(&r,0x5a,sizeof(r));int status=fusion_c_pb_thermal_birth_budget_grid(T,o,p,n,e,nullbirth?nullptr:r.b.data(),nullout?nullptr:&r.out,nullreport?nullptr:&r.report);
  need(status!=0,"invalid input accepted");if(!nullout)need(zero(r.out),"failure out cleared");if(!nullreport)need(zero(r.report),"failure report cleared");if(!nullbirth&&n>=1&&n<=100000)need(zero(r.b),"failure birth cleared");
 };
 for(auto T:{nan,inf,-kev,0.})reject(T,&base,&good,cells,edges.data());
 reject(kev,nullptr,&good,cells,edges.data());reject(kev,&base,nullptr,cells,edges.data());reject(kev,&base,&good,cells,nullptr);
 reject(kev,&base,&good,cells,edges.data(),true);reject(kev,&base,&good,cells,edges.data(),false,true);reject(kev,&base,&good,cells,edges.data(),false,false,true);
 reject(kev,&base,&good,0,edges.data());reject(kev,&base,&good,100001,edges.data());
 for(int mode=0;mode<10;++mode){auto p=good;switch(mode){case 0:p.version=2;break;case 1:p.reserved=1;break;case 2:p.coarse_nq=3;break;case 3:p.coarse_ncos=9;break;case 4:p.number_L1_tolerance=-.1;break;case 5:p.energy_L1_tolerance=2.1;break;case 6:p.number_L1_tolerance=nan;break;case 7:p.energy_L1_tolerance=inf;break;case 8:p.coarse_nq=9;break;case 9:p.coarse_ncos=3;break;}reject(kev,&base,&p,cells,edges.data());}
 for(int mode=0;mode<6;++mode){auto e=edges;if(mode==0)e[2]=e[1];if(mode==1)e[2]=nan;if(mode==2)e[0]=-1;if(mode==3)e[8]=inf;if(mode==4)e[2]=e[1]/2;if(mode==5)e[0]=nan;reject(kev,&base,&good,cells,e.data());}
 for(int mode=0;mode<5;++mode){auto o=base;if(mode==0)o.relative_max_J=nan;if(mode==1)o.relative_order=3;if(mode==2)o.nq=3;if(mode==3)o.cm_max_kT=7;if(mode==4)o.relative_max_J=20*mev;reject(kev,&o,&good,cells,edges.data());}
}
}
int main(){try{real_cases();rejection();std::cout<<"PASS budget API discrete-source/positive-tail/zero-budget/error contracts; fixed engineering allowance1e-12, not interval or physical convergence certificate\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
