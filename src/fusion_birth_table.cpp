#include "fusion_birth_table.h"
#include "fusion_birth_table_internal.h"
#include "fusion_nuclear_data.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
namespace {
using R=long double;
constexpr int BAD=PB11_STATUS_INVALID_ARGUMENT,NUM=PB11_STATUS_NUMERICAL_FAILURE;
struct Failure {int status;};
void require(bool p,int status=NUM){if(!p)throw Failure{status};}
template<class C,class F>void fields(C&c,F f){
 f(c.reactivity_m3_s);
 for(auto&v:c.reactant_energy_moment_J_m3_s)f(v);
 for(auto&v:c.below_number_m3_s)f(v);
 for(auto&v:c.below_energy_J_m3_s)f(v);
 for(auto&v:c.above_number_m3_s)f(v);
 for(auto&v:c.above_energy_J_m3_s)f(v);
}
double stored(R v){require(std::isfinite(v)&&v>=0&&v<=std::numeric_limits<double>::max());return double(v);}
struct Node {double T=0;std::vector<double> grid;fusion_birth_coefficients_v1 c{};};
using Ptr=std::shared_ptr<Node>;
using Errors=std::array<double,4>;
R weight(double T,double a,double b){
 const R denominator=std::log(R(b)/a);require(std::isfinite(denominator)&&denominator>0);
 const R w=std::log(R(T)/a)/denominator;require(std::isfinite(w)&&w>=0&&w<=1);return w;
}
double geometric(double a,double b,R fraction){
 const R value=std::exp((1-fraction)*std::log(R(a))+fraction*std::log(R(b)));
 const double out=stored(value);require(out>a&&out<b);return out;
}
fusion_birth_coefficients_v1 mix(const Node&a,const Node&b,R w){
 std::array<double,31> av{},bv{};size_t i=0;fields(a.c,[&](double v){av[i++]=v;});
 i=0;fields(b.c,[&](double v){bv[i++]=v;});
 fusion_birth_coefficients_v1 c{};i=0;fields(c,[&](double&v){v=stored((1-w)*av[i]+w*bv[i]);++i;});return c;
}
R relative(R difference,R reference){
 if(reference==0)return difference==0?0:std::numeric_limits<R>::infinity();
 return difference/reference;
}
bool balanced(R a,R b){return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=1e-10L*(std::abs(a)+std::abs(b));}
}
struct fusion_birth_table_v1 {
 fusion_birth_table_info_v1 info{};
 std::vector<double> edges;
 std::vector<Ptr> knots;
 fusion_nuclear_channel_v1 channel{};
};
namespace {
void conservative(const fusion_birth_table_v1&t,const double*grid,const fusion_birth_coefficients_v1&c){
 bool valid=true;fields(c,[&](double v){if(!std::isfinite(v)||v<0)valid=false;});require(valid);
 const int n=t.info.cells;R totalE=0;
 for(int species=0;species<7;++species){
  R N=R(c.below_number_m3_s[species])+c.above_number_m3_s[species];
  R E=R(c.below_energy_J_m3_s[species])+c.above_energy_J_m3_s[species];
  for(int j=0;j<n;++j){const double g=grid[species*n+j];require(std::isfinite(g)&&g>=0);N+=g;E+=R(g)*(R(t.edges[j])+t.edges[j+1])/2;}
  int multiplicity=0;for(int j=0;j<t.channel.product_count;++j){int id=t.channel.product_ids[j];if(id==FUSION_MASS_NEUTRON)id=6;if(id==species)++multiplicity;}
  require(balanced(N,R(multiplicity)*c.reactivity_m3_s));totalE+=E;
 }
 require(balanced(totalE,R(t.channel.q_J)*c.reactivity_m3_s+c.reactant_energy_moment_J_m3_s[0]+c.reactant_energy_moment_J_m3_s[1]));
}
Ptr direct(fusion_birth_table_v1&t,double T){
 require(t.info.direct_evaluations<t.info.control.max_evaluations);
 ++t.info.direct_evaluations;
 auto p=std::make_shared<Node>();p->T=T;p->grid.resize(7*t.info.cells);
 fusion_thermal_birth_v1 r{};int st=fusion_c_thermal_birth_grid(t.info.channel,T,&t.info.source,t.info.cells,t.edges.data(),p->grid.data(),&r);
 require(st==PB11_STATUS_OK,st);
 require(std::abs(r.relative_rate_discrepancy)<=t.info.control.max_direct_rate_discrepancy&&r.relative_reactant_energy_discrepancy<=t.info.control.max_direct_debit_discrepancy);
 t.info.max_sampled_direct_rate_discrepancy=std::max(t.info.max_sampled_direct_rate_discrepancy,std::abs(r.relative_rate_discrepancy));
 t.info.max_sampled_direct_debit_discrepancy=std::max(t.info.max_sampled_direct_debit_discrepancy,r.relative_reactant_energy_discrepancy);
 auto&c=p->c;c.reactivity_m3_s=r.reactivity_m3_s;
 for(int i=0;i<2;++i)c.reactant_energy_moment_J_m3_s[i]=r.reactant_energy_moment_J_m3_s[i];
 for(int i=0;i<7;++i){c.below_number_m3_s[i]=r.below_number_m3_s[i];c.below_energy_J_m3_s[i]=r.below_energy_J_m3_s[i];c.above_number_m3_s[i]=r.above_number_m3_s[i];c.above_energy_J_m3_s[i]=r.above_energy_J_m3_s[i];}
 conservative(t,p->grid.data(),c);return p;
}
Errors discrepancy(const fusion_birth_table_v1&t,const Node&a,const Node&b,const Node&reference){
 const R w=weight(reference.T,a.T,b.T);const auto c=mix(a,b,w);const auto&r=reference.c;
 Errors e{};
 e[0]=double(relative(std::abs(R(c.reactivity_m3_s)-r.reactivity_m3_s),r.reactivity_m3_s));
 for(int j=0;j<2;++j)e[1]=std::max(e[1],double(relative(std::abs(R(c.reactant_energy_moment_J_m3_s[j])-r.reactant_energy_moment_J_m3_s[j]),r.reactant_energy_moment_J_m3_s[j])));
 const int n=t.info.cells;
 for(int species=0;species<7;++species){
  R dN=0,dE=0,N=0,E=0;
  dN+=std::abs(R(c.below_number_m3_s[species])-r.below_number_m3_s[species])+std::abs(R(c.above_number_m3_s[species])-r.above_number_m3_s[species]);
  dE+=std::abs(R(c.below_energy_J_m3_s[species])-r.below_energy_J_m3_s[species])+std::abs(R(c.above_energy_J_m3_s[species])-r.above_energy_J_m3_s[species]);
  N+=R(r.below_number_m3_s[species])+r.above_number_m3_s[species];E+=R(r.below_energy_J_m3_s[species])+r.above_energy_J_m3_s[species];
  for(int j=0;j<n;++j){int k=species*n+j;R center=(R(t.edges[j])+t.edges[j+1])/2;
   R v=stored((1-w)*a.grid[k]+w*b.grid[k]),diff=std::abs(v-reference.grid[k]);dN+=diff;dE+=center*diff;N+=reference.grid[k];E+=center*reference.grid[k];}
  e[2]=std::max(e[2],double(relative(dN,N)));e[3]=std::max(e[3],double(relative(dE,E)));
 }
 return e;
}
void refine(fusion_birth_table_v1&t,const Ptr&a,const Ptr&b,int depth){
 const auto&q=t.info.control;Errors worst{};Ptr midpoint;
 for(int j=1;j<=3;++j){auto ref=direct(t,geometric(a->T,b->T,R(j)/4));const auto e=discrepancy(t,*a,*b,*ref);for(int i=0;i<4;++i)worst[i]=std::max(worst[i],e[i]);if(j==2)midpoint=ref;}
 const bool pass=worst[0]<=q.max_rate_error&&worst[1]<=q.max_debit_error&&worst[2]<=q.max_number_L1&&worst[3]<=q.max_energy_L1;
 if(pass){require(t.knots.size()+2<=size_t(q.max_knots));t.knots.push_back(a);
  t.info.max_validated_rate_error=std::max(t.info.max_validated_rate_error,worst[0]);
  t.info.max_validated_debit_error=std::max(t.info.max_validated_debit_error,worst[1]);
  t.info.max_validated_number_L1=std::max(t.info.max_validated_number_L1,worst[2]);
  t.info.max_validated_energy_L1=std::max(t.info.max_validated_energy_L1,worst[3]);return;
 }
 require(depth<q.max_depth);refine(t,a,midpoint,depth+1);refine(t,midpoint,b,depth+1);
}
}
extern "C" int fusion_c_birth_table_create(int channel,double lower,double upper,
 const fusion_thermal_birth_options_v1*source,const fusion_birth_table_control_v1*control,
 int n,const double*edges,fusion_birth_table_v1**out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out=nullptr;
 if(!source||!control||!edges||channel<0||channel>4||n<1||n>100000)return BAD;
 if(!std::isfinite(lower)||!std::isfinite(upper)||lower<=0||upper<=lower)return BAD;
 const auto&q=*control;
 for(double v:{q.max_rate_error,q.max_debit_error,q.max_number_L1,q.max_energy_L1,q.max_direct_rate_discrepancy,q.max_direct_debit_discrepancy})if(!std::isfinite(v)||v<0||v>1)return BAD;
 if(q.max_knots<2||q.max_knots>100000||q.max_evaluations<5||q.max_evaluations>1000000||q.max_depth<0||q.max_depth>24)return BAD;
 if(7ULL*n*(q.max_knots+3*q.max_depth+6)>50000000ULL)return BAD;
 for(int j=0;j<=n;++j)if(!std::isfinite(edges[j])||edges[j]<0||(j&&edges[j]<=edges[j-1]))return BAD;
 try{auto t=std::make_unique<fusion_birth_table_v1>();auto&i=t->info;i.channel=channel;i.cells=n;i.lower_kT_J=lower;i.upper_kT_J=upper;i.source=*source;i.control=q;t->edges.assign(edges,edges+n+1);
  int st=fusion_c_nuclear_channel(channel,&t->channel);require(st==PB11_STATUS_OK,st);
  auto a=direct(*t,lower),b=direct(*t,upper);refine(*t,a,b,0);t->knots.push_back(b);i.knots=int(t->knots.size());*out=t.release();return PB11_STATUS_OK;
 }catch(const Failure&f){return f.status;}catch(...){return PB11_STATUS_EXCEPTION;}
}
extern "C" void fusion_c_birth_table_destroy(fusion_birth_table_v1*t){delete t;}
extern "C" int fusion_c_birth_table_info(const fusion_birth_table_v1*t,fusion_birth_table_info_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};if(!t)return BAD;*out=t->info;return PB11_STATUS_OK;
}
extern "C" int fusion_c_birth_table_evaluate(const fusion_birth_table_v1*t,double T,int n,double*grid,fusion_birth_coefficients_v1*out){
 if(out)*out={};
 if(grid&&n>0&&n<=100000)std::fill(grid,grid+7*n,0.);
 if(!out||!grid)return PB11_STATUS_NULL_OUTPUT;
 if(!t||n!=t->info.cells||!std::isfinite(T))return BAD;
 if(T<t->info.lower_kT_J||T>t->info.upper_kT_J)return PB11_STATUS_OUT_OF_RANGE;
 try{auto it=std::lower_bound(t->knots.begin(),t->knots.end(),T,[](const Ptr&p,double x){return p->T<x;});require(it!=t->knots.end());
  std::vector<double> values;fusion_birth_coefficients_v1 c{};
  if((*it)->T==T){values=(*it)->grid;c=(*it)->c;}
  else{require(it!=t->knots.begin());const auto&a=**(it-1);const auto&b=**it;R w=weight(T,a.T,b.T);c=mix(a,b,w);values.resize(7*n);for(int j=0;j<7*n;++j)values[j]=stored((1-w)*a.grid[j]+w*b.grid[j]);}
  conservative(*t,values.data(),c);std::copy(values.begin(),values.end(),grid);*out=c;return PB11_STATUS_OK;
 }catch(const Failure&f){return f.status;}catch(...){return PB11_STATUS_EXCEPTION;}
}

namespace fusion_detail {
bool birth_table_matches(const fusion_birth_table_v1*t,int channel,
 const fusion_thermal_birth_options_v1&s,int n,const double*edges) noexcept {
 if(!t||!edges||channel!=t->info.channel||n!=t->info.cells)return false;
 const auto&a=t->info.source;
 if(a.relative_max_J!=s.relative_max_J||a.cm_max_kT!=s.cm_max_kT||
    a.ground_state_q_J!=s.ground_state_q_J||a.cutoff_J!=s.cutoff_J||
    a.l1_fraction!=s.l1_fraction||a.relative_phase!=s.relative_phase||
    a.narrow_peak_fraction!=s.narrow_peak_fraction||a.continuum_peak_scale!=s.continuum_peak_scale||
    a.continuation!=s.continuation||a.pb_low!=s.pb_low||a.remainder_policy!=s.remainder_policy||
    a.broad_mode!=s.broad_mode||a.fsci_policy!=s.fsci_policy||a.relative_order!=s.relative_order||
    a.cm_order!=s.cm_order||a.nq!=s.nq||a.ncos!=s.ncos)return false;
 return std::equal(t->edges.begin(),t->edges.end(),edges);
}
}
