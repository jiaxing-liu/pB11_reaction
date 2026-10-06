#include "fusion_birth_table.h"
#include "fusion_birth_table_internal.h"
#include "fusion_thermal_parent_internal.h"
#include "fusion_cache_identity_internal.h"
#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <type_traits>
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
void refine(fusion_birth_table_v1&t,const Ptr&a,const Ptr&b,int depth,Ptr midpoint_hint={}){
 const auto&q=t.info.control;Errors worst{};Ptr midpoint,quarter1,quarter3;
 for(int j=1;j<=3;++j){const double T=geometric(a->T,b->T,R(j)/4);auto ref=(j==2&&midpoint_hint&&T==midpoint_hint->T)?midpoint_hint:direct(t,T);const auto e=discrepancy(t,*a,*b,*ref);for(int i=0;i<4;++i)worst[i]=std::max(worst[i],e[i]);if(j==1)quarter1=ref;else if(j==2)midpoint=ref;else quarter3=ref;}
 const bool pass=worst[0]<=q.max_rate_error&&worst[1]<=q.max_debit_error&&worst[2]<=q.max_number_L1&&worst[3]<=q.max_energy_L1;
 if(pass){require(t.knots.size()+2<=size_t(q.max_knots));t.knots.push_back(a);
  t.info.max_validated_rate_error=std::max(t.info.max_validated_rate_error,worst[0]);
  t.info.max_validated_debit_error=std::max(t.info.max_validated_debit_error,worst[1]);
  t.info.max_validated_number_L1=std::max(t.info.max_validated_number_L1,worst[2]);
  t.info.max_validated_energy_L1=std::max(t.info.max_validated_energy_L1,worst[3]);return;
 }
 require(depth<q.max_depth);refine(t,a,midpoint,depth+1,quarter1);refine(t,midpoint,b,depth+1,quarter3);
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

// Thermal-table persistence deliberately shares the generated source/data
// identity with beam tables but uses its own magic and a dense knot payload.
namespace {
constexpr uint64_t thermal_magic=0x3154544e4f495346ULL;
constexpr uint64_t thermal_format=1;
constexpr size_t thermal_byte_cap=512ULL*1024*1024;
constexpr uint64_t thermal_double_budget=50000000ULL;
constexpr size_t thermal_header_bytes=24+64;
constexpr size_t thermal_info_words=38;
constexpr double thermal_mev=1.602176634e-13;
static_assert(sizeof(double)==8&&std::numeric_limits<double>::is_iec559,
              "binary64 required");

template<class I,class F,class G>void table_fields(I&i,F real,G integer){
 real(i.lower_kT_J);real(i.upper_kT_J);
 real(i.max_validated_rate_error);real(i.max_validated_debit_error);
 real(i.max_validated_number_L1);real(i.max_validated_energy_L1);
 real(i.max_sampled_direct_rate_discrepancy);
 real(i.max_sampled_direct_debit_discrepancy);
 integer(i.channel);integer(i.cells);integer(i.knots);integer(i.direct_evaluations);
 auto&s=i.source;
 real(s.relative_max_J);real(s.cm_max_kT);real(s.ground_state_q_J);real(s.cutoff_J);
 real(s.l1_fraction);real(s.relative_phase);real(s.narrow_peak_fraction);
 real(s.continuum_peak_scale);
 integer(s.continuation);integer(s.pb_low);integer(s.remainder_policy);
 integer(s.broad_mode);integer(s.fsci_policy);integer(s.relative_order);
 integer(s.cm_order);integer(s.nq);integer(s.ncos);
 auto&c=i.control;
 real(c.max_rate_error);real(c.max_debit_error);real(c.max_number_L1);
 real(c.max_energy_L1);real(c.max_direct_rate_discrepancy);
 real(c.max_direct_debit_discrepancy);
 integer(c.max_knots);integer(c.max_evaluations);integer(c.max_depth);
}

uint64_t thermal_checksum(const unsigned char*p,size_t n){
 uint64_t h=14695981039346656037ULL;
 for(size_t i=0;i<n;++i){h^=p[i];h*=1099511628211ULL;}
 return h;
}
struct ThermalWriter{
 unsigned char*p;size_t n,pos=0;
 void bytes(const void*src,size_t count){
  require(pos<=n&&count<=n-pos,BAD);std::memcpy(p+pos,src,count);pos+=count;
 }
 void integer(uint64_t v){
  unsigned char b[8];for(int i=0;i<8;++i)b[i]=static_cast<unsigned char>(v>>(8*i));
  bytes(b,8);
 }
 void real(double v){
  uint64_t bits=0;std::memcpy(&bits,&v,8);integer(bits);
 }
};
struct ThermalReader{
 const unsigned char*p;size_t n,pos=0;
 void bytes(void*dst,size_t count){
  require(pos<=n&&count<=n-pos,BAD);std::memcpy(dst,p+pos,count);pos+=count;
 }
 uint64_t integer(){
  unsigned char b[8];bytes(b,8);uint64_t v=0;
  for(int i=0;i<8;++i)v|=uint64_t(b[i])<<(8*i);
  return v;
 }
 double real(){
  uint64_t bits=integer();double v=0;std::memcpy(&v,&bits,8);return v;
 }
};

size_t checked_bytes(size_t base,size_t count,size_t unit){
 require(base<=thermal_byte_cap,BAD);
 require(unit==0||count<=(thermal_byte_cap-base)/unit,BAD);
 return base+count*unit;
}
size_t thermal_wire_size(int cells,int knots){
 require(cells>=1&&cells<=100000&&knots>=2&&knots<=100000,BAD);
 const size_t grid_words=32+7*size_t(cells);
 const size_t node_bytes=checked_bytes(0,grid_words,8);
 size_t size=thermal_header_bytes+thermal_info_words*8+8;
 size=checked_bytes(size,size_t(cells)+1,8);
 size=checked_bytes(size,size_t(knots),node_bytes);
 return size;
}

void validate_source(const fusion_birth_table_info_v1&i,
                     const fusion_nuclear_channel_v1&reaction){
 const auto&s=i.source;
 bool finite=true;
 for(double v:{s.relative_max_J,s.cm_max_kT,s.ground_state_q_J,s.cutoff_J,
               s.l1_fraction,s.relative_phase,s.narrow_peak_fraction,
               s.continuum_peak_scale})
  finite=finite&&std::isfinite(v);
 require(finite,BAD);
 require(s.relative_max_J>0&&s.cm_max_kT>=8&&s.cm_max_kT<=80&&
         s.ground_state_q_J>=0&&s.cutoff_J>=.001*thermal_mev&&
         s.cutoff_J<=.01*thermal_mev&&s.l1_fraction>=0&&s.l1_fraction<=1&&
         s.narrow_peak_fraction>=0&&s.narrow_peak_fraction<=1&&
         s.continuum_peak_scale>=0&&s.continuum_peak_scale<=2,BAD);
 require(s.continuation>=1&&s.continuation<=2&&s.pb_low>=0&&s.pb_low<=3&&
         (i.channel==FUSION_PB11_3ALPHA||s.pb_low==FUSION_PB_LOW_TB),BAD);
 require(s.remainder_policy>=0&&s.remainder_policy<=2&&
         (s.broad_mode==1||s.broad_mode==3||s.broad_mode==13)&&
         s.fsci_policy>=0&&s.fsci_policy<=1,BAD);
 require(s.relative_order>=4&&s.relative_order<=64&&
         s.cm_order>=4&&s.cm_order<=32&&s.nq>=4&&s.nq<=1024&&
         s.ncos>=4&&s.ncos<=1024,BAD);
 require(i.channel!=FUSION_PB11_3ALPHA||
         s.ground_state_q_J<=reaction.q_J,BAD);
}

void validate_table_metadata(const fusion_birth_table_info_v1&i,
                             const fusion_nuclear_channel_v1&reaction){
 bool finite=true;
 table_fields(i,
              [&](double v){finite=finite&&std::isfinite(v);},
              [](auto){});
 require(finite,BAD);
 require(i.channel>=0&&i.channel<FUSION_CHANNEL_COUNT&&
         i.cells>=1&&i.cells<=100000&&i.knots>=2&&i.knots<=100000&&
         i.lower_kT_J>0&&i.upper_kT_J>i.lower_kT_J,BAD);
 const auto&q=i.control;
 for(double v:{q.max_rate_error,q.max_debit_error,q.max_number_L1,
               q.max_energy_L1,q.max_direct_rate_discrepancy,
               q.max_direct_debit_discrepancy})
  require(std::isfinite(v)&&v>=0&&v<=1,BAD);
 require(q.max_knots>=2&&q.max_knots<=100000&&
         q.max_evaluations>=5&&q.max_evaluations<=1000000&&
         q.max_depth>=0&&q.max_depth<=24,BAD);
 require(i.knots<=q.max_knots&&
         i.direct_evaluations>=std::max(5,i.knots)&&
         i.direct_evaluations<=q.max_evaluations,BAD);
 const uint64_t term=uint64_t(q.max_knots)+3*uint64_t(q.max_depth)+6;
 const uint64_t slots=7*uint64_t(i.cells);
 require(term>0&&slots<=thermal_double_budget/term,BAD);
 const double measured[]={i.max_validated_rate_error,
  i.max_validated_debit_error,i.max_validated_number_L1,
  i.max_validated_energy_L1,i.max_sampled_direct_rate_discrepancy,
  i.max_sampled_direct_debit_discrepancy};
 const double gates[]={q.max_rate_error,q.max_debit_error,
  q.max_number_L1,q.max_energy_L1,q.max_direct_rate_discrepancy,
  q.max_direct_debit_discrepancy};
 for(int j=0;j<6;++j)require(std::isfinite(measured[j])&&
                              measured[j]>=0&&measured[j]<=gates[j],BAD);
 validate_source(i,reaction);
}

void validate_edges(const fusion_birth_table_v1&t){
 require(t.edges.size()==size_t(t.info.cells)+1,BAD);
 for(size_t j=0;j<t.edges.size();++j)
  require(std::isfinite(t.edges[j])&&t.edges[j]>=0&&
          (j==0||t.edges[j]>t.edges[j-1]),BAD);
}
void validate_node(const fusion_birth_table_v1&t,const Node&node){
 fusion_nuclear_mass_v1 ma{},mb{};
 require(fusion_c_nuclear_mass(t.channel.reactant_ids[0],&ma)==PB11_STATUS_OK,BAD);
 require(fusion_c_nuclear_mass(t.channel.reactant_ids[1],&mb)==PB11_STATUS_OK,BAD);
 require(fusion_detail::thermal_parent_preflight(t.info.channel,t.info.source.relative_max_J,
   R(node.T)*t.info.source.cm_max_kT,ma.mass_kg,mb.mass_kg,t.channel.q_J,
   t.info.source.ground_state_q_J)==PB11_STATUS_OK,BAD);
 const size_t count=7*size_t(t.info.cells);
 require(node.grid.size()==count,BAD);
 bool finite=true;
 fields(node.c,[&](double v){finite=finite&&std::isfinite(v)&&v>=0;});
 require(finite,BAD);
 for(double v:node.grid)require(std::isfinite(v)&&v>=0,BAD);
 conservative(t,node.grid.data(),node.c);
}
void validate_table_structure(const fusion_birth_table_v1&t){
 fusion_nuclear_channel_v1 reaction{};
 require(fusion_c_nuclear_channel(t.info.channel,&reaction)==PB11_STATUS_OK,BAD);
 validate_table_metadata(t.info,reaction);
 require(t.knots.size()==size_t(t.info.knots),BAD);
 require(t.channel.product_count==reaction.product_count&&
         t.channel.q_J==reaction.q_J,BAD);
 validate_edges(t);
 double previous=0;
 for(size_t i=0;i<t.knots.size();++i){
  const auto&node=t.knots[i];require(bool(node),BAD);
  require(std::isfinite(node->T)&&node->T>=t.info.lower_kT_J&&
          node->T<=t.info.upper_kT_J&&(i==0||node->T>previous),BAD);
  validate_node(t,*node);previous=node->T;
 }
 require(!t.knots.empty()&&t.knots.front()->T==t.info.lower_kT_J&&
         t.knots.back()->T==t.info.upper_kT_J,BAD);
}
size_t thermal_packed_size(const fusion_birth_table_v1&t){
 validate_table_structure(t);
 return thermal_wire_size(t.info.cells,t.info.knots);
}
}

extern "C" int fusion_c_birth_table_matches_request(
 const fusion_birth_table_v1*t,int channel,double lower,double upper,
 const fusion_thermal_birth_options_v1*source,
 const fusion_birth_table_control_v1*control,int cells,const double*edges,
 int*matches){
 if(!matches)return PB11_STATUS_NULL_OUTPUT;
 *matches=0;
 if(!t||!source||!control||!edges||cells<1||cells>100000)return BAD;
 if(!fusion_detail::birth_table_matches(t,channel,*source,cells,edges))
  return PB11_STATUS_OK;
 const auto&a=t->info.control;const auto&b=*control;
 *matches=t->info.lower_kT_J==lower&&t->info.upper_kT_J==upper&&
  a.max_rate_error==b.max_rate_error&&a.max_debit_error==b.max_debit_error&&
  a.max_number_L1==b.max_number_L1&&a.max_energy_L1==b.max_energy_L1&&
  a.max_direct_rate_discrepancy==b.max_direct_rate_discrepancy&&
  a.max_direct_debit_discrepancy==b.max_direct_debit_discrepancy&&
  a.max_knots==b.max_knots&&a.max_evaluations==b.max_evaluations&&
  a.max_depth==b.max_depth;
 return PB11_STATUS_OK;
}
extern "C" const char*fusion_c_birth_table_kernel_identity(void){
 return fusion_detail::beam_cache_kernel_identity;
}
extern "C" int fusion_c_birth_table_pack_size(
 const fusion_birth_table_v1*t,size_t*required){
 if(!required)return PB11_STATUS_NULL_OUTPUT;
 *required=0;if(!t)return BAD;
 try{*required=thermal_packed_size(*t);return PB11_STATUS_OK;}
 catch(const Failure&f){return f.status;}catch(...){return PB11_STATUS_EXCEPTION;}
}
extern "C" int fusion_c_birth_table_pack(
 const fusion_birth_table_v1*t,void*buffer,size_t capacity,size_t*written){
 if(!written)return PB11_STATUS_NULL_OUTPUT;
 *written=0;if(!t||!buffer)return BAD;
 try{
  const size_t size=thermal_packed_size(*t);require(capacity>=size,BAD);
  ThermalWriter w{static_cast<unsigned char*>(buffer),size};
  w.integer(thermal_magic);w.integer(thermal_format);w.integer(size);
  w.bytes(fusion_detail::beam_cache_kernel_identity,64);
  table_fields(t->info,[&](double v){w.real(v);},
               [&](auto v){w.integer(uint64_t(v));});
  for(double v:t->edges)w.real(v);
  for(const auto&node:t->knots){
   w.real(node->T);fields(node->c,[&](double v){w.real(v);});
   for(double v:node->grid)w.real(v);
  }
  require(w.pos==size-8,BAD);
  w.integer(thermal_checksum(w.p,w.pos));*written=size;
  return PB11_STATUS_OK;
 }catch(const Failure&f){return f.status;}catch(...){return PB11_STATUS_EXCEPTION;}
}
extern "C" int fusion_c_birth_table_unpack(
 const void*buffer,size_t length,fusion_birth_table_v1**out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out=nullptr;
 if(!buffer||length<thermal_header_bytes+thermal_info_words*8+8||
    length>thermal_byte_cap)return BAD;
 try{
  const auto*data=static_cast<const unsigned char*>(buffer);
  ThermalReader footer{data+length-8,8};
  require(footer.integer()==thermal_checksum(data,length-8),BAD);
  ThermalReader r{data,length-8};
  require(r.integer()==thermal_magic&&r.integer()==thermal_format&&
          r.integer()==length,BAD);
  char identity[64];r.bytes(identity,64);
  require(std::memcmp(identity,fusion_detail::beam_cache_kernel_identity,64)==0,BAD);
  auto t=std::make_unique<fusion_birth_table_v1>();
  table_fields(t->info,[&](double&v){v=r.real();},
   [&](auto&v){uint64_t x=r.integer();using V=std::decay_t<decltype(v)>;
    require(x<=uint64_t(std::numeric_limits<V>::max()),BAD);
    v=static_cast<V>(x);});
  fusion_nuclear_channel_v1 reaction{};
  require(fusion_c_nuclear_channel(t->info.channel,&reaction)==PB11_STATUS_OK,BAD);
  validate_table_metadata(t->info,reaction);
  require(thermal_wire_size(t->info.cells,t->info.knots)==length,BAD);
  t->channel=reaction;
  t->edges.resize(size_t(t->info.cells)+1);
  for(size_t j=0;j<t->edges.size();++j){
   const double v=r.real();
   require(std::isfinite(v)&&v>=0&&(j==0||v>t->edges[j-1]),BAD);
   t->edges[j]=v;
  }
  t->knots.reserve(size_t(t->info.knots));double previous=0;
  const size_t count=7*size_t(t->info.cells);
  for(int i=0;i<t->info.knots;++i){
   auto node=std::make_shared<Node>();node->T=r.real();
   require(std::isfinite(node->T)&&node->T>=t->info.lower_kT_J&&
           node->T<=t->info.upper_kT_J&&(i==0||node->T>previous),BAD);
   fields(node->c,[&](double&v){v=r.real();});
   node->grid.resize(count);
   for(double&v:node->grid)v=r.real();
   validate_node(*t,*node);t->knots.push_back(node);previous=node->T;
  }
  require(r.pos==r.n&&t->knots.front()->T==t->info.lower_kT_J&&
          t->knots.back()->T==t->info.upper_kT_J,BAD);
  *out=t.release();return PB11_STATUS_OK;
 }catch(const Failure&f){return f.status;}catch(...){return PB11_STATUS_EXCEPTION;}
}
