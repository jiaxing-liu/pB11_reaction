#include "fusion_thermal_birth.h"
#include "fusion_beam_birth.h"
#include "fusion_alpha_events.h"
#include "fusion_rate_model.h"
#include "fusion_pb_population.h"
#include "fusion_pb_population_data.h"
#include "fusion_nuclear_data.h"
#include "fusion_reaction_event.h"
#include "fusion_products.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
namespace {
using R=long double;
constexpr R c=299792458.L,c2=c*c,pi=3.141592653589793238462643383279502884L;
constexpr double kev=1.602176634e-16,mev=1.602176634e-13;
bool finite_value(double x){return std::isfinite(x);}
struct Mapper {
 int n;std::vector<R> center,birth;
 std::array<R,7> bn{},be{},an{},ae{},number{},energy{};
 Mapper(int cells,const double*e):n(cells),center(n),birth(7*n){
  for(int j=0;j<n;++j)center[j]=(R(e[j])+e[j+1])/2;
 }
 void delta(int id,R e,R w){
  if(w==0)return;
  number[id]+=w;energy[id]+=w*e;
  if(e<center.front()){bn[id]+=w;be[id]+=w*e;return;}
  if(e>center.back()){an[id]+=w;ae[id]+=w*e;return;}
  auto p=std::lower_bound(center.begin(),center.end(),e);int j=int(p-center.begin());
  if(*p==e)birth[id*n+j]+=w;
  else {R f=(e-center[j-1])/(center[j]-center[j-1]);birth[id*n+j]+=w*f;birth[id*n+j-1]+=w*(1-f);}
 }
 void box(int id,R lo,R hi,R w){
  if(w==0)return;
  if(lo==hi){delta(id,lo,w);return;}
  if(!(lo>=0&&hi>lo&&w>0))throw std::runtime_error("Invalid birth box");
  number[id]+=w;energy[id]+=w*(lo+(hi-lo)/2);
  auto part=[&](R l,R u,R& num,R& en){l=std::max(l,lo);u=std::min(u,hi);
   if(u>l){R z=w*(u-l)/(hi-lo);num+=z;en+=z*(l+(u-l)/2);}};
  part(lo,center.front(),bn[id],be[id]);part(center.back(),hi,an[id],ae[id]);
  int start=std::max(1,int(std::upper_bound(center.begin(),center.end(),lo)-center.begin()));
  for(int j=start;j<n&&center[j-1]<hi;++j){
   R l=std::max(lo,center[j-1]),u=std::min(hi,center[j]);if(u<=l)continue;
   R z=w*(u-l)/(hi-lo),mid=l+(u-l)/2,f=(mid-center[j-1])/(center[j]-center[j-1]);
   birth[id*n+j]+=z*f;birth[id*n+j-1]+=z*(1-f);
  }
 }
 struct Boost {R beta,gamma,d;explicit Boost(R b):beta(b){R root=std::sqrt(1-beta*beta);gamma=1/root;d=beta*beta/(root*(1+root));}};
 void isotropic(int id,R mass,R K,const Boost&boost,R w){
  const R beta=boost.beta,gamma=boost.gamma,d=boost.d,rest=mass*c2;
  R mid=gamma*K+d*rest,half=gamma*beta*std::sqrt(K*(K+2*rest));
  R hi=mid+half;
  // (mid-half)*(mid+half)=(K-(gamma-1)*mc^2)^2 avoids
  // catastrophic cancellation near a product brought to rest in the lab.
  R lo=hi==0?0:(K-d*rest)*(K-d*rest)/hi;
  box(id,lo,hi,w);
 }
};
struct CMNode {R kinetic,weight;};
std::vector<CMNode> cm_nodes(double T,const fusion_thermal_birth_options_v1&o){
 const auto g=fusion_detail::gauss_legendre(o.cm_order);
 std::vector<R> knots{0};R upper=std::sqrt(R(o.cm_max_kT));
 for(R v:{.25L,.5L,1.L,2.L,3.L,4.L,6.L,8.L})if(v<upper)knots.push_back(v);
 knots.push_back(upper);std::vector<CMNode> result;
 for(std::size_t i=1;i<knots.size();++i){R mid=(knots[i]+knots[i-1])/2,h=(knots[i]-knots[i-1])/2;
  for(auto q:g){R y=mid+h*q.x;result.push_back({R(T)*y*y,h*q.w*4/std::sqrt(pi)*y*y*std::exp(-y*y)});}}
 return result;
}
std::vector<double> relative_knots(double T,double limit,int ch){
 std::vector<double> k{0,limit};
 for(R e=R(T)/128;e<limit;e*=2){if(e>0)k.push_back(double(e));else break;}
 double lo=0,hi=0;int st=fusion_c_cross_section_domain(ch,&lo,&hi);if(st)throw std::runtime_error("Channel domain");
 for(double e:{lo,hi})if(e>0&&e<limit)k.push_back(e);
 // Split at internal piecewise cross-section fit boundaries, as the
 // independent rate reference does; a single rule must not straddle them.
 if(ch==FUSION_DT_ALPHAN&&530*kev<limit)k.push_back(530*kev);
 if(ch==FUSION_DHE3_ALPHAP&&900*kev<limit)k.push_back(900*kev);
 if(ch==0){
  for(double x:{101.,124.5,138.6,143.3,145.65,148.,150.35,152.7,157.4,171.5,195.,400.,668.,1211.,2340.,3294.,5700.})
   if(x*kev<limit)k.push_back(x*kev);
  // Population-table interpolation knots, converted from equivalent p lab
  // energy to relative energy using the same canonical masses as its API.
  fusion_nuclear_mass_v1 p{},b{};fusion_c_nuclear_mass(0,&p);fusion_c_nuclear_mass(5,&b);
  for(double x:pb_population_data::lab_keV){
   double e=x*kev*b.mass_kg/(p.mass_kg+b.mass_kg);if(e<limit)k.push_back(e);
  }
 }
 std::sort(k.begin(),k.end());k.erase(std::unique(k.begin(),k.end()),k.end());return k;
}
// The NR event has zero total momentum. One common p^2 multiplier preserves
// that identity while putting all alphas on shell at the requested total A.
std::array<R,3> on_shell(const fusion_detail::AlphaEvent&e,R A,R rest,R& shift){
 R old=R(e.energy_J[0])+e.energy_J[1]+e.energy_J[2];
 R lo=A/old,hi=lo*(1+A/(2*rest)),s=(lo+hi)/2;
 std::array<R,3> K{};bool done=false;
 for(int it=0;it<24;++it){R sum=0,der=0;
  for(int j=0;j<3;++j){R x=s*e.energy_J[j],root=std::sqrt(rest*rest+2*rest*x);
   K[j]=2*rest*x/(root+rest);sum+=K[j];der+=rest*e.energy_J[j]/root;}
  R error=sum-A;if(std::abs(error)<4e-17L*A){done=true;break;}
  if(error>0)hi=s;else lo=s;R next=s-error/der;
  s=(next>lo&&next<hi)?next:(lo+hi)/2;
 }
 if(!done)throw std::runtime_error("On-shell CM remapping did not converge");
 for(int j=0;j<3;++j)shift=std::max(shift,std::abs(K[j]-R(e.energy_J[j])*A/old)/A);
 return K;
}
struct SourceWeights {R f0=0,flow=0,fbroad=0;fusion_detail::AlphaEvents low{},broad{};};
int prepare_weights(int ch,R E,const fusion_thermal_birth_options_v1&o,
 const fusion_nuclear_channel_v1&reaction,SourceWeights&weights){
 int st=0;
    const double A0=double(R(reaction.q_J)+E);
    if(ch==0){fusion_pb_population_v1 pop{};
     st=fusion_c_pb_population(o.continuation,o.pb_low,double(E),o.narrow_peak_fraction,o.continuum_peak_scale,&pop);if(st)return st;
     weights.f0=pop.alpha0_peak_fraction;
     if(o.remainder_policy==0){weights.flow=pop.narrow_remainder_fraction;weights.fbroad=pop.other_remainder_fraction;}
     else if(o.remainder_policy==1)weights.flow=R(pop.narrow_remainder_fraction)+pop.other_remainder_fraction;
     else weights.fbroad=R(pop.narrow_remainder_fraction)+pop.other_remainder_fraction;
     if(weights.flow>0){st=fusion_detail::alpha_events(2,o.fsci_policy,A0,o.cutoff_J,o.l1_fraction,o.relative_phase,o.nq,o.ncos,weights.low);if(st)return st;}
     if(weights.fbroad>0){st=fusion_detail::alpha_events(o.broad_mode,o.fsci_policy,A0,o.cutoff_J,o.l1_fraction,o.relative_phase,o.nq,o.ncos,weights.broad);if(st)return st;}
    }
 return PB11_STATUS_OK;
}
int emit_products(int ch,R A,R beta,R w,const fusion_thermal_birth_options_v1&o,
 const fusion_nuclear_channel_v1&reaction,const fusion_nuclear_mass_v1*products,
 const std::vector<fusion_detail::QuadNode>&angular,const SourceWeights&weights,Mapper&map,R&shell_shift){
 const Mapper::Boost boost(beta);
 int st=0;
     if(ch!=0){fusion_particle_four_vector_v1 pair[2]{};double direction[3]={0,0,1};
      st=fusion_c_two_body_cm(products[0].mass_kg,products[1].mass_kg,double(A),direction,pair);if(st)return st;
      for(int j=0;j<2;++j)map.isotropic(reaction.product_ids[j],products[j].mass_kg,pair[j].kinetic_energy_J,boost,w);
     }else{
      if(weights.f0>0)for(auto angle:angular){fusion_three_body_cm_v1 event{};
       st=fusion_c_three_equal_sequential_cm(products[0].mass_kg,double(A),o.ground_state_q_J,angle.x,&event);if(st)return st;
       for(double K:event.kinetic_energy_J)map.isotropic(4,products[0].mass_kg,K,boost,w*weights.f0*angle.w/2);
      }
      auto alpha1=[&](const fusion_detail::AlphaEvents&events,R fraction){
       if(fraction==0)return;
       for(const auto&event:events.events){if(event.weight==0)continue;
        auto K=on_shell(event,A,R(products[0].mass_kg)*c2,shell_shift);
        for(R k:K)map.isotropic(4,products[0].mass_kg,k,boost,w*fraction*event.weight);
       }
      };
      alpha1(weights.low,weights.flow);alpha1(weights.broad,weights.fbroad);
     }
 return PB11_STATUS_OK;
}
int finish_source(Mapper&map,int n,const fusion_nuclear_channel_v1&reaction,
 R rate,R ea,R eb,fusion_thermal_birth_v1&result,double*birth,fusion_thermal_birth_v1*out){
  if(result.reference_reactivity_m3_s>0)result.relative_rate_discrepancy=double(rate/result.reference_reactivity_m3_s-1);
  else if(rate>0)return PB11_STATUS_NUMERICAL_FAILURE;
  for(int j=0;j<2;++j){double ref=result.reference_reactant_energy_moment_J_m3_s[j];
   if(ref>0)result.relative_reactant_energy_discrepancy=std::max(result.relative_reactant_energy_discrepancy,std::abs(result.reactant_energy_moment_J_m3_s[j]/ref-1));
   else if(result.reactant_energy_moment_J_m3_s[j]>0)return PB11_STATUS_NUMERICAL_FAILURE;
  }
  std::vector<double> mapped(7*n);R number_error=0,energy_total=0;
  for(int id=0;id<7;++id){R N=0,U=0;
   for(int j=0;j<n;++j){R value=map.birth[id*n+j];if(!std::isfinite(value)||value<0||value>std::numeric_limits<double>::max())return PB11_STATUS_NUMERICAL_FAILURE;
    mapped[id*n+j]=double(value);N+=mapped[id*n+j];U+=R(mapped[id*n+j])*map.center[j];}
   result.below_number_m3_s[id]=double(map.bn[id]);result.below_energy_J_m3_s[id]=double(map.be[id]);
   result.above_number_m3_s[id]=double(map.an[id]);result.above_energy_J_m3_s[id]=double(map.ae[id]);
   N+=R(result.below_number_m3_s[id])+result.above_number_m3_s[id];
   U+=R(result.below_energy_J_m3_s[id])+result.above_energy_J_m3_s[id];
   int multiplicity=0;for(int j=0;j<reaction.product_count;++j)if(reaction.product_ids[j]==id)++multiplicity;
   number_error=std::max(number_error,std::abs(N-multiplicity*rate));energy_total+=U;
  }
  R expected=R(reaction.q_J)*rate+ea+eb,energy_error=energy_total-expected;
  result.product_energy_moment_J_m3_s=double(energy_total);result.number_residual_m3_s=double(number_error);
  result.energy_residual_J_m3_s=double(energy_error);
  if(!std::isfinite(expected)||!std::isfinite(energy_total)||!std::isfinite(rate)||
     number_error>1e-10L*rate||std::abs(energy_error)>1e-10L*expected)return PB11_STATUS_NUMERICAL_FAILURE;
  // Long-double accumulation can remain finite while a public double
  // coefficient overflows. Validate the ABI representation before publishing.
  for(double x:{result.reactivity_m3_s,result.reference_reactivity_m3_s,
      result.product_energy_moment_J_m3_s,result.number_residual_m3_s,
      result.energy_residual_J_m3_s,result.relative_rate_discrepancy,
      result.relative_reactant_energy_discrepancy,result.cm_retained_probability,
      result.cm_tail_probability,result.cm_tail_energy_moment_J,
      result.max_cm_energy_shift_fraction,result.max_shell_remap_fraction})
    if(!std::isfinite(x))return PB11_STATUS_NUMERICAL_FAILURE;
  for(int j=0;j<2;++j)
    if(!std::isfinite(result.reactant_energy_moment_J_m3_s[j])||
       !std::isfinite(result.reference_reactant_energy_moment_J_m3_s[j]))
      return PB11_STATUS_NUMERICAL_FAILURE;
  for(int j=0;j<7;++j)
    if(!std::isfinite(result.below_number_m3_s[j])||!std::isfinite(result.below_energy_J_m3_s[j])||
       !std::isfinite(result.above_number_m3_s[j])||!std::isfinite(result.above_energy_J_m3_s[j]))
      return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(mapped.begin(),mapped.end(),birth);*out=result;return PB11_STATUS_OK;
}
int parent(int ch,R E,R C,double ma,double mb,fusion_reaction_parent_v1&out){
 R M=R(ma)+mb,mu=R(ma)*mb/M,v=std::sqrt(2*C/M),p=std::sqrt(2*mu*E);
 double pa[3]={double(p),0,double(ma*v)},pb[3]={-double(p),0,double(mb*v)};
 return fusion_c_reaction_parent(ch,FUSION_REACTANT_CLASSICAL_BUDGET,pa,pb,&out);
}
}
static int pair_birth_impl(int ch,double Ta,double Tb,int correlation_order,const fusion_thermal_birth_options_v1*op,
 int n,const double*edges,double*birth,fusion_thermal_birth_v1*out){
 double T=Ta;
 if(out)*out={};
 if(n>0&&n<=100000&&birth)std::fill(birth,birth+7*n,0.);
 if(!out||!birth)return PB11_STATUS_NULL_OUTPUT;
 if(!op||!edges||n<1||n>100000||!finite_value(Ta)||!finite_value(Tb))return PB11_STATUS_INVALID_ARGUMENT;
 if(correlation_order<4||correlation_order>32)return PB11_STATUS_INVALID_ARGUMENT;
 const auto&o=*op;
 for(double x:{o.relative_max_J,o.cm_max_kT,o.ground_state_q_J,o.cutoff_J,o.l1_fraction,o.relative_phase,o.narrow_peak_fraction,o.continuum_peak_scale})
  if(!finite_value(x))return PB11_STATUS_INVALID_ARGUMENT;
 if(Ta<=0||Tb<=0||o.relative_max_J<=0||o.cm_max_kT<8||o.cm_max_kT>80||o.ground_state_q_J<0||
    o.cutoff_J<.001*mev||o.cutoff_J>.01*mev||o.l1_fraction<0||o.l1_fraction>1||
    o.narrow_peak_fraction<0||o.narrow_peak_fraction>1||o.continuum_peak_scale<0||o.continuum_peak_scale>2)
  return PB11_STATUS_OUT_OF_RANGE;
 if(o.relative_order<4||o.relative_order>64||o.cm_order<4||o.cm_order>32||o.nq<4||o.nq>1024||o.ncos<4||o.ncos>1024||
   o.remainder_policy<0||o.remainder_policy>2||(o.broad_mode!=1&&o.broad_mode!=3&&o.broad_mode!=13)||o.fsci_policy<0||o.fsci_policy>1)
  return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<=n;++j)if(!finite_value(edges[j])||edges[j]<0||(j&&edges[j]<=edges[j-1]))return PB11_STATUS_INVALID_ARGUMENT;
 try{
  fusion_nuclear_channel_v1 reaction{};int st=fusion_c_nuclear_channel(ch,&reaction);if(st)return st;
  fusion_nuclear_mass_v1 ma{},mb{},products[3]{};
  fusion_c_nuclear_mass(reaction.reactant_ids[0],&ma);fusion_c_nuclear_mass(reaction.reactant_ids[1],&mb);
  for(int j=0;j<reaction.product_count;++j)fusion_c_nuclear_mass(reaction.product_ids[j],&products[j]);
  const bool correlated=Ta!=Tb;
  const R mass_sum=R(ma.mass_kg)+mb.mass_kg;
  const R reduced=R(ma.mass_kg)*mb.mass_kg/mass_sum;
  const R velocity_variance=R(Ta)/ma.mass_kg+R(Tb)/mb.mass_kg;
  const R drift=(R(Ta)-Tb)/(mass_sum*velocity_variance);
  double cm_temperature=Ta;
  if(correlated){
   T=double(reduced*velocity_variance);
   cm_temperature=double(mass_sum*R(Ta)*Tb/(R(ma.mass_kg)*Tb+R(mb.mass_kg)*Ta));
   if(!finite_value(T)||!finite_value(cm_temperature)||T<=0||cm_temperature<=0)
    return PB11_STATUS_NUMERICAL_FAILURE;
  }
  fusion_rate_model_v1 reference{};
  st=fusion_c_thermal_pair_maxwellian_model(ch,o.continuation,o.pb_low,ma.mass_kg,mb.mass_kg,Ta,Tb,&reference);if(st)return st;
  fusion_reaction_parent_v1 maximum{};
  R max_cm=R(T)*o.cm_max_kT;
  if(correlated){
   R speed=std::abs(drift)*std::sqrt(2*R(o.relative_max_J)/reduced)+
     std::sqrt(2*R(cm_temperature)*o.cm_max_kT/mass_sum);
   max_cm=mass_sum*speed*speed/2;
  }
  st=parent(ch,o.relative_max_J,max_cm,ma.mass_kg,mb.mass_kg,maximum);if(st)return st;
  if(ch==0&&(maximum.available_cm_energy_J>12*mev||o.ground_state_q_J>reaction.q_J))return PB11_STATUS_OUT_OF_RANGE;
  auto cm=cm_nodes(cm_temperature,o);auto nodes=fusion_detail::gauss_legendre(o.relative_order);
  auto angular=fusion_detail::gauss_legendre(o.ncos);auto knots=relative_knots(T,o.relative_max_J,ch);
  auto directions=correlated?fusion_detail::gauss_legendre(correlation_order):
    std::vector<fusion_detail::QuadNode>{{0,2}};
  Mapper map(n,edges);fusion_thermal_birth_v1 result{};
  R rate=0,ea=0,eb=0,max_shift=0,shell_shift=0;
  R M=R(ma.mass_kg)+mb.mass_kg,mu=R(ma.mass_kg)*mb.mass_kg/M;
  R prefactor=std::sqrt(8/(pi*mu))/(R(T)*std::sqrt(R(T)));
  R cmprob=0;for(auto q:cm)cmprob+=q.weight;
  for(std::size_t interval=1;interval<knots.size();++interval){
   R mid=(R(knots[interval])+knots[interval-1])/2,h=(R(knots[interval])-knots[interval-1])/2;
   for(auto q:nodes){R E=mid+h*q.x;double sigma=0;
    st=fusion_c_cross_section_model(ch,o.continuation,o.pb_low,double(E),&sigma);if(st)return st;
    R wr=h*q.w*prefactor*E*sigma*std::exp(-E/T);if(wr==0)continue;
    // Distributions far below representable coefficient precision contribute
    // no stored source, just as the scalar model; do not spend on spectra.
    if(double(wr)==0)continue;
    SourceWeights weights;
    st=prepare_weights(ch,E,o,reaction,weights);if(st)return st;
    for(auto cq:cm)for(auto direction:directions){R w=wr*cq.weight;
     if(correlated)w*=R(direction.w)/2;
     if(w==0)continue;
     fusion_reaction_parent_v1 par{};
     if(correlated){
      // u lies along x; the residual CM Gaussian has polar cosine direction.x.
      R u=std::sqrt(2*E/mu),v=std::sqrt(2*cq.kinetic/M);
      R vx=drift*u+v*direction.x,vy=v*std::sqrt(1-R(direction.x)*direction.x);
      double pa[3]={double(ma.mass_kg*(vx+R(mb.mass_kg)/M*u)),double(ma.mass_kg*vy),0};
      double pb[3]={double(mb.mass_kg*(vx-R(ma.mass_kg)/M*u)),double(mb.mass_kg*vy),0};
      st=fusion_c_reaction_parent(ch,FUSION_REACTANT_CLASSICAL_BUDGET,pa,pb,&par);
     }else st=parent(ch,E,cq.kinetic,ma.mass_kg,mb.mass_kg,par);
     if(st)return st;
     R A=par.available_cm_energy_J,beta=std::abs(R(par.boost_velocity_m_s[2]))/c;
     if(correlated){R v2=0;for(double v:par.boost_velocity_m_s)v2+=R(v)*v;beta=std::sqrt(v2)/c;}
     rate+=w;
     if(correlated){ea+=w*par.classical_kinetic_J[0];eb+=w*par.classical_kinetic_J[1];}
     else {ea+=w*(R(ma.mass_kg)/M*cq.kinetic+R(mb.mass_kg)/M*E);
      eb+=w*(R(mb.mass_kg)/M*cq.kinetic+R(ma.mass_kg)/M*E);}
     max_shift=std::max(max_shift,std::abs(A-(R(reaction.q_J)+E))/(R(reaction.q_J)+E));
     st=emit_products(ch,A,beta,w,o,reaction,products,angular,weights,map,shell_shift);if(st)return st;
    }
   }
  }
  result.reactivity_m3_s=double(rate);result.reference_reactivity_m3_s=reference.total.resolved_reactivity_m3_s;
  result.reactant_energy_moment_J_m3_s[0]=double(ea);result.reactant_energy_moment_J_m3_s[1]=double(eb);
  result.reference_reactant_energy_moment_J_m3_s[0]=reference.total.projectile_energy_reactivity_J_m3_s;
  result.reference_reactant_energy_moment_J_m3_s[1]=reference.total.target_energy_reactivity_J_m3_s;
  result.cm_retained_probability=double(cmprob);
  R X=o.cm_max_kT,tail=std::erfc(std::sqrt(X))+2/std::sqrt(pi)*std::sqrt(X)*std::exp(-X);
  result.cm_tail_probability=double(tail);
  result.cm_tail_energy_moment_J=double(R(cm_temperature)*(1.5L*tail+2/std::sqrt(pi)*X*std::sqrt(X)*std::exp(-X)));
  result.max_cm_energy_shift_fraction=double(max_shift);result.max_shell_remap_fraction=double(shell_shift);
  return finish_source(map,n,reaction,rate,ea,eb,result,birth,out);
 }catch(...){return PB11_STATUS_EXCEPTION;}
}

extern "C" int fusion_c_thermal_birth_grid(int ch,double T,const fusion_thermal_birth_options_v1*op,
 int n,const double*edges,double*birth,fusion_thermal_birth_v1*out){
 return pair_birth_impl(ch,T,T,4,op,n,edges,birth,out);
}
extern "C" int fusion_c_thermal_pair_birth_grid(int ch,double Ta,double Tb,int correlation_order,
 const fusion_thermal_birth_options_v1*op,int n,const double*edges,double*birth,fusion_thermal_birth_v1*out){
 return pair_birth_impl(ch,Ta,Tb,correlation_order,op,n,edges,birth,out);
}

extern "C" int fusion_c_beam_birth_grid(int ch,int slot,double projectile_E,double target_T,
 const fusion_beam_birth_options_v1*op,int n,const double*edges,double*birth,fusion_beam_birth_v1*out){
 if(out)*out={};
 if(n>0&&n<=100000&&birth)std::fill(birth,birth+7*n,0.);
 if(!out||!birth)return PB11_STATUS_NULL_OUTPUT;
 if(!op||!edges||n<1||n>100000||slot<0||slot>1||!finite_value(projectile_E)||!finite_value(target_T))
  return PB11_STATUS_INVALID_ARGUMENT;
 if(projectile_E<0||target_T<0)return PB11_STATUS_OUT_OF_RANGE;
 const auto&b=*op;
 const fusion_thermal_birth_options_v1 o{b.relative_max_J,b.angular_max_exponent,b.ground_state_q_J,b.cutoff_J,
  b.l1_fraction,b.relative_phase,b.narrow_peak_fraction,b.continuum_peak_scale,b.continuation,b.pb_low,
  b.remainder_policy,b.broad_mode,b.fsci_policy,b.relative_order,b.angular_order,b.nq,b.ncos};
 for(double x:{o.relative_max_J,o.cm_max_kT,o.ground_state_q_J,o.cutoff_J,o.l1_fraction,o.relative_phase,o.narrow_peak_fraction,o.continuum_peak_scale})
  if(!finite_value(x))return PB11_STATUS_INVALID_ARGUMENT;
 if(o.relative_max_J<=0||o.cm_max_kT<8||o.cm_max_kT>80||o.ground_state_q_J<0||
    o.cutoff_J<.001*mev||o.cutoff_J>.01*mev||o.l1_fraction<0||o.l1_fraction>1||
    o.narrow_peak_fraction<0||o.narrow_peak_fraction>1||o.continuum_peak_scale<0||o.continuum_peak_scale>2)
  return PB11_STATUS_OUT_OF_RANGE;
 if(o.relative_order<4||o.relative_order>64||o.cm_order<4||o.cm_order>32||o.nq<4||o.nq>1024||o.ncos<4||o.ncos>1024||
    o.remainder_policy<0||o.remainder_policy>2||(o.broad_mode!=1&&o.broad_mode!=3&&o.broad_mode!=13)||o.fsci_policy<0||o.fsci_policy>1)
  return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<=n;++j)if(!finite_value(edges[j])||edges[j]<0||(j&&edges[j]<=edges[j-1]))return PB11_STATUS_INVALID_ARGUMENT;
 try{
  fusion_nuclear_channel_v1 reaction{};int st=fusion_c_nuclear_channel(ch,&reaction);if(st)return st;
  fusion_nuclear_mass_v1 reactants[2]{},products[3]{};
  for(int j=0;j<2;++j){st=fusion_c_nuclear_mass(reaction.reactant_ids[j],&reactants[j]);if(st)return st;}
  for(int j=0;j<reaction.product_count;++j){st=fusion_c_nuclear_mass(reaction.product_ids[j],&products[j]);if(st)return st;}
  R mp=reactants[slot].mass_kg,mt=reactants[1-slot].mass_kg,M=mp+mt,mu=mp*mt/M;
  R v=std::sqrt(2*R(projectile_E)/mp);
  fusion_rate_model_v1 reference{};
  st=fusion_c_beam_maxwellian_model(ch,o.continuation,o.pb_low,double(mp),double(mt),projectile_E,target_T,&reference);if(st)return st;
  if(ch==0&&o.ground_state_q_J>reaction.q_J)return PB11_STATUS_OUT_OF_RANGE;
  auto parent_at=[&](R w,R cosine,fusion_reaction_parent_v1&par){
   R transverse=w*std::sqrt(1-cosine*cosine);
   double pp[3]={0,0,double(mp*v)},pt[3]={double(-mt*transverse),0,double(mt*(v-w*cosine))};
   return fusion_c_reaction_parent(ch,FUSION_REACTANT_CLASSICAL_BUDGET,slot?pt:pp,slot?pp:pt,&par);
  };
  Mapper map(n,edges);fusion_beam_birth_v1 result{};
  auto angular=fusion_detail::gauss_legendre(o.ncos);
  R rate=0,ea=0,eb=0,max_shift=0,shell_shift=0,relative_probability=0,retained_probability=0,omitted_probability=0;
  auto contribute=[&](R E,R w,R cosine,R weight,const SourceWeights&weights)->int{
   if(weight==0)return PB11_STATUS_OK;
   fusion_reaction_parent_v1 par{};int status=parent_at(w,cosine,par);if(status)return status;
   R A=par.available_cm_energy_J;
   if(ch==0&&A>12*mev)return PB11_STATUS_OUT_OF_RANGE;
   R speed2=0;for(double x:par.boost_velocity_m_s)speed2+=R(x)*x;
   rate+=weight;
   // Match the exact caller bin energy for the fast debit. Incoming momentum
   // conversion differs only by rounding, handled by the unchanged closure gate.
   R e0=slot?R(par.classical_kinetic_J[0]):R(projectile_E);
   R e1=slot?R(projectile_E):R(par.classical_kinetic_J[1]);
   ea+=weight*e0;eb+=weight*e1;
   max_shift=std::max(max_shift,std::abs(A-(R(reaction.q_J)+E))/(R(reaction.q_J)+E));
   return emit_products(ch,A,std::sqrt(speed2)/c,weight,o,reaction,products,angular,weights,map,shell_shift);
  };
  if(target_T==0){
   R E=mu*v*v/2;
   if(E<=o.relative_max_J){
    relative_probability=retained_probability=1;
    double sigma=0;st=fusion_c_cross_section_model(ch,o.continuation,o.pb_low,double(E),&sigma);if(st)return st;
    R weight=R(sigma)*v;
    if(weight>0){SourceWeights weights;st=prepare_weights(ch,E,o,reaction,weights);if(st)return st;
     st=contribute(E,v,1,weight,weights);if(st)return st;}
   }
  }else{
   R u=std::sqrt(2*R(target_T)/mt),s=v/u;
   if(!std::isfinite(s)||s>1e6L)return PB11_STATUS_NUMERICAL_FAILURE;
   double relative_T=double(mu/mt*target_T);
   if(!finite_value(relative_T)||relative_T<=0)return PB11_STATUS_NUMERICAL_FAILURE;
   auto knots=relative_knots(relative_T,o.relative_max_J,ch);
   for(R delta:{-40.L,-16.L,-8.L,-4.L,-2.L,-1.L,0.L,1.L,2.L,4.L,8.L,16.L,40.L}){
    R w=v+delta*u;if(w<=0)continue;R E=mu*w*w/2;
    if(E>0&&E<o.relative_max_J)knots.push_back(double(E));
   }
   std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
   auto nodes=fusion_detail::gauss_legendre(o.relative_order),polar=fusion_detail::gauss_legendre(o.cm_order);
   // At fixed w, the classical parent invariant grows with decreasing cosine;
   // wmax and cosine=-1 bound the retained domain, including angular tails.
   fusion_reaction_parent_v1 maximum{};
   st=parent_at(std::sqrt(2*R(o.relative_max_J)/mu),-1,maximum);if(st)return st;
   if(ch==0&&maximum.available_cm_energy_J>12*mev)return PB11_STATUS_OUT_OF_RANGE;
   for(std::size_t interval=1;interval<knots.size();++interval){
    R mid=(R(knots[interval])+knots[interval-1])/2,h=(R(knots[interval])-knots[interval-1])/2;
    for(auto q:nodes){
     R E=mid+h*q.x,w=std::sqrt(2*E/mu),x=w/u,kappa=2*x*s;
     R factor=kappa==0?1:-std::expm1(-2*kappa)/(2*kappa);
     R radial=4*x*x/std::sqrt(pi)/u*std::exp(-(x-s)*(x-s))*factor;
     R probability=h*q.w*radial/(mu*w);if(probability==0)continue;
     relative_probability+=probability;
     double sigma=0;st=fusion_c_cross_section_model(ch,o.continuation,o.pb_low,double(E),&sigma);if(st)return st;
     R wr=probability*w*sigma;
     // Match the coefficient representability convention of thermal birth.
     if(double(wr)==0)wr=0;
     SourceWeights weights;if(wr>0){st=prepare_weights(ch,E,o,reaction,weights);if(st)return st;}
     auto angle_at=[&](R cosine,R angle_weight){
      retained_probability+=probability*angle_weight;
      return contribute(E,w,cosine,wr*angle_weight,weights);
     };
     if(kappa<1){
      R norm=kappa==0?.5L:kappa/(-std::expm1(-2*kappa));
      for(auto a:polar){st=angle_at(a.x,R(a.w)*norm*std::exp(kappa*(R(a.x)-1)));if(st)return st;}
     }else{
      R upper=std::min(2*kappa,R(o.cm_max_kT)),norm=-std::expm1(-2*kappa);
      if(upper<2*kappa)omitted_probability+=probability*(std::exp(-upper)-std::exp(-2*kappa))/norm;
      std::vector<R> cuts{0};for(R z:{.25L,.5L,1.L,2.L,4.L,8.L,16.L,32.L,64.L})if(z<upper)cuts.push_back(z);cuts.push_back(upper);
      for(std::size_t j=1;j<cuts.size();++j){R center=(cuts[j]+cuts[j-1])/2,width=(cuts[j]-cuts[j-1])/2;
       for(auto a:polar){R z=center+width*a.x;
        st=angle_at(1-z/kappa,width*a.w*std::exp(-z)/norm);if(st)return st;}
      }
     }
    }
   }
  }
  auto&r=result.spectrum;
  r.reactivity_m3_s=double(rate);r.reference_reactivity_m3_s=reference.total.resolved_reactivity_m3_s;
  r.reactant_energy_moment_J_m3_s[0]=double(ea);r.reactant_energy_moment_J_m3_s[1]=double(eb);
  r.reference_reactant_energy_moment_J_m3_s[slot]=reference.total.projectile_energy_reactivity_J_m3_s;
  r.reference_reactant_energy_moment_J_m3_s[1-slot]=reference.total.target_energy_reactivity_J_m3_s;
  r.max_cm_energy_shift_fraction=double(max_shift);r.max_shell_remap_fraction=double(shell_shift);
  result.relative_retained_probability=double(relative_probability);
  result.retained_pair_probability=double(retained_probability);
  result.angular_omitted_pair_probability=double(omitted_probability);
  for(double x:{result.relative_retained_probability,result.retained_pair_probability,result.angular_omitted_pair_probability})
   if(!finite_value(x)||x<0)return PB11_STATUS_NUMERICAL_FAILURE;
  // Finish into a temporary grid: even a late failure cannot publish a source.
  std::vector<double> mapped(7*n);fusion_thermal_birth_v1 finalized{};
  st=finish_source(map,n,reaction,rate,ea,eb,r,mapped.data(),&finalized);if(st)return st;
  result.spectrum=finalized;std::copy(mapped.begin(),mapped.end(),birth);*out=result;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
