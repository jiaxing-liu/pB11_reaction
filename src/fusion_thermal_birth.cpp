#include "fusion_thermal_birth.h"
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
 void isotropic(int id,R mass,R K,R beta,R w){
  R root=std::sqrt(1-beta*beta),gamma=1/root,d=beta*beta/(root*(1+root)),rest=mass*c2;
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
int parent(int ch,R E,R C,double ma,double mb,fusion_reaction_parent_v1&out){
 R M=R(ma)+mb,mu=R(ma)*mb/M,v=std::sqrt(2*C/M),p=std::sqrt(2*mu*E);
 double pa[3]={double(p),0,double(ma*v)},pb[3]={-double(p),0,double(mb*v)};
 return fusion_c_reaction_parent(ch,FUSION_REACTANT_CLASSICAL_BUDGET,pa,pb,&out);
}
}
extern "C" int fusion_c_thermal_birth_grid(int ch,double T,const fusion_thermal_birth_options_v1*op,
 int n,const double*edges,double*birth,fusion_thermal_birth_v1*out){
 if(out)*out={};
 if(n>0&&n<=100000&&birth)std::fill(birth,birth+7*n,0.);
 if(!out||!birth)return PB11_STATUS_NULL_OUTPUT;
 if(!op||!edges||n<1||n>100000||!finite_value(T))return PB11_STATUS_INVALID_ARGUMENT;
 const auto&o=*op;
 for(double x:{o.relative_max_J,o.cm_max_kT,o.ground_state_q_J,o.cutoff_J,o.l1_fraction,o.relative_phase,o.narrow_peak_fraction,o.continuum_peak_scale})
  if(!finite_value(x))return PB11_STATUS_INVALID_ARGUMENT;
 if(T<=0||o.relative_max_J<=0||o.cm_max_kT<8||o.cm_max_kT>80||o.ground_state_q_J<0||
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
  fusion_rate_model_v1 reference{};
  st=fusion_c_thermal_pair_maxwellian_model(ch,o.continuation,o.pb_low,ma.mass_kg,mb.mass_kg,T,T,&reference);if(st)return st;
  fusion_reaction_parent_v1 maximum{};
  st=parent(ch,o.relative_max_J,R(T)*o.cm_max_kT,ma.mass_kg,mb.mass_kg,maximum);if(st)return st;
  if(ch==0&&(maximum.available_cm_energy_J>12*mev||o.ground_state_q_J>reaction.q_J))return PB11_STATUS_OUT_OF_RANGE;
  auto cm=cm_nodes(T,o);auto nodes=fusion_detail::gauss_legendre(o.relative_order);
  auto angular=fusion_detail::gauss_legendre(o.ncos);auto knots=relative_knots(T,o.relative_max_J,ch);
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
    R f0=0,flow=0,fbroad=0;
    fusion_detail::AlphaEvents low{},broad{};
    const double A0=double(R(reaction.q_J)+E);
    if(ch==0){fusion_pb_population_v1 pop{};
     st=fusion_c_pb_population(o.continuation,o.pb_low,double(E),o.narrow_peak_fraction,o.continuum_peak_scale,&pop);if(st)return st;
     f0=pop.alpha0_peak_fraction;
     if(o.remainder_policy==0){flow=pop.narrow_remainder_fraction;fbroad=pop.other_remainder_fraction;}
     else if(o.remainder_policy==1)flow=R(pop.narrow_remainder_fraction)+pop.other_remainder_fraction;
     else fbroad=R(pop.narrow_remainder_fraction)+pop.other_remainder_fraction;
     if(flow>0){st=fusion_detail::alpha_events(2,o.fsci_policy,A0,o.cutoff_J,o.l1_fraction,o.relative_phase,o.nq,o.ncos,low);if(st)return st;}
     if(fbroad>0){st=fusion_detail::alpha_events(o.broad_mode,o.fsci_policy,A0,o.cutoff_J,o.l1_fraction,o.relative_phase,o.nq,o.ncos,broad);if(st)return st;}
    }
    for(auto cq:cm){R w=wr*cq.weight;if(w==0)continue;
     fusion_reaction_parent_v1 par{};st=parent(ch,E,cq.kinetic,ma.mass_kg,mb.mass_kg,par);if(st)return st;
     R A=par.available_cm_energy_J,beta=std::abs(R(par.boost_velocity_m_s[2]))/c;
     rate+=w;ea+=w*(R(ma.mass_kg)/M*cq.kinetic+R(mb.mass_kg)/M*E);
     eb+=w*(R(mb.mass_kg)/M*cq.kinetic+R(ma.mass_kg)/M*E);
     max_shift=std::max(max_shift,std::abs(A-(R(reaction.q_J)+E))/(R(reaction.q_J)+E));
     if(ch!=0){fusion_particle_four_vector_v1 pair[2]{};double direction[3]={0,0,1};
      st=fusion_c_two_body_cm(products[0].mass_kg,products[1].mass_kg,double(A),direction,pair);if(st)return st;
      for(int j=0;j<2;++j)map.isotropic(reaction.product_ids[j],products[j].mass_kg,pair[j].kinetic_energy_J,beta,w);
     }else{
      if(f0>0)for(auto angle:angular){fusion_three_body_cm_v1 event{};
       st=fusion_c_three_equal_sequential_cm(products[0].mass_kg,double(A),o.ground_state_q_J,angle.x,&event);if(st)return st;
       for(double K:event.kinetic_energy_J)map.isotropic(4,products[0].mass_kg,K,beta,w*f0*angle.w/2);
      }
      auto alpha1=[&](const fusion_detail::AlphaEvents&events,R fraction){
       if(fraction==0)return;
       for(const auto&event:events.events){if(event.weight==0)continue;
        auto K=on_shell(event,A,R(products[0].mass_kg)*c2,shell_shift);
        for(R k:K)map.isotropic(4,products[0].mass_kg,k,beta,w*fraction*event.weight);
       }
      };
      alpha1(low,flow);alpha1(broad,fbroad);
     }
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
  result.cm_tail_energy_moment_J=double(R(T)*(1.5L*tail+2/std::sqrt(pi)*X*std::sqrt(X)*std::exp(-X)));
  result.max_cm_energy_shift_fraction=double(max_shift);result.max_shell_remap_fraction=double(shell_shift);
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
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
