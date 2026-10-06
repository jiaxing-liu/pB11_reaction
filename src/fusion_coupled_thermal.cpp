#include "fusion_source_rounding_internal.h"
#include "fusion_coupled_sources.h"
#include "fusion_beam_birth_table_internal.h"
#include "fusion_beam_birth.h"
#include "fusion_target_network.h"
#include "fusion_birth_table_internal.h"
#include "fusion_thermal_burn.h"
#include "fusion_nuclear_data.h"
#include "fusion_two_component.h"
#include "fusion_handoff.h"
#include "fusion_rate_model.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <vector>
namespace {
using R=long double;
constexpr int BAD=PB11_STATUS_INVALID_ARGUMENT,NUM=PB11_STATUS_NUMERICAL_FAILURE;
constexpr double MeV=1.602176634e-13;
bool finite_value(double x){return std::isfinite(x);}
bool nonnegative(double x){return finite_value(x)&&x>=0;}
bool trial_timing_enabled(){
 static int enabled=-1;
 if(enabled<0){const char*v=std::getenv("BALDUR_FUSION_TRIAL_TIMING");enabled=(v&&*v)?1:0;}
 return enabled==1;
}
void trial_stage(const char*name){
 if(trial_timing_enabled()){std::fprintf(stdout,"FUSION_CXX_STAGE %s\n",name);std::fflush(stdout);}
}
bool put(R x,double& d){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;d=double(x);return true;}
using fusion_detail::source_rounding::floor_source_rate;
using fusion_detail::source_rounding::source_roundoff_accumulate;
bool close(std::initializer_list<R> terms,R roundoff=0){R value=0,scale=0;for(R x:terms){if(!std::isfinite(x))return false;value+=x;scale+=std::abs(x);}return std::abs(value)<=1e-10L*scale+roundoff;}
template<class F>void ledger_fields(fusion_source_ledger_v1&l,F f){
 for(double&x:l.events_m3)f(x);
 for(double&x:l.nuclear_born_number_m3)f(x);
 for(double&x:l.nuclear_born_energy_J_m3)f(x);
 for(double&x:l.external_born_number_m3)f(x);
 for(double&x:l.external_born_energy_J_m3)f(x);
 for(double&x:l.thermal_consumed_number_m3)f(x);
 for(double&x:l.thermal_consumed_energy_J_m3)f(x);
 for(double&x:l.fast_consumed_number_m3)f(x);
 for(double&x:l.fast_consumed_energy_J_m3)f(x);
 for(double&x:l.escaped_number_m3)f(x);
 for(double&x:l.escaped_energy_J_m3)f(x);
 for(double&x:l.handed_off_number_m3)f(x);
 for(double&x:l.handed_off_energy_J_m3)f(x);
 for(double&x:l.heat_to_bath_J_m3)f(x);
 f(l.neutron_number_m3);f(l.neutron_energy_J_m3);
}
fusion_beam_birth_options_v1 beam_options_for(const fusion_thermal_birth_options_v1&b,
 const fusion_fast_target_options_v1&fast,int channel){
 fusion_beam_birth_options_v1 result{};
 result.relative_max_J=b.relative_max_J;result.angular_max_exponent=fast.angular_max_exponent;
 result.ground_state_q_J=b.ground_state_q_J;result.cutoff_J=b.cutoff_J;
 result.l1_fraction=b.l1_fraction;result.relative_phase=b.relative_phase;
 result.narrow_peak_fraction=b.narrow_peak_fraction;result.continuum_peak_scale=b.continuum_peak_scale;
 result.continuation=b.continuation;result.pb_low=channel==0?b.pb_low:FUSION_PB_LOW_TB;
 result.remainder_policy=b.remainder_policy;result.broad_mode=b.broad_mode;result.fsci_policy=b.fsci_policy;
 result.relative_order=b.relative_order;result.angular_order=fast.angular_order;result.nq=b.nq;result.ncos=b.ncos;
 return result;
}
int validate_options(const fusion_coupled_thermal_options_v1&o){
 const auto&b=o.birth;
 for(double x:{b.relative_max_J,b.cm_max_kT,b.ground_state_q_J,b.cutoff_J,b.l1_fraction,b.relative_phase,b.narrow_peak_fraction,b.continuum_peak_scale,o.max_source_rate_error,o.max_source_debit_error,o.handoff_max_L1,o.handoff_max_mean_error})if(!finite_value(x))return BAD;
 if(b.relative_max_J<=0||b.cm_max_kT<8||b.cm_max_kT>80||b.ground_state_q_J<0||b.cutoff_J<.001*MeV||b.cutoff_J>.01*MeV||b.l1_fraction<0||b.l1_fraction>1||b.narrow_peak_fraction<0||b.narrow_peak_fraction>1||b.continuum_peak_scale<0||b.continuum_peak_scale>2||o.max_source_rate_error<0||o.max_source_rate_error>1||o.max_source_debit_error<0||o.max_source_debit_error>1||o.handoff_max_L1<0||o.handoff_max_L1>2||o.handoff_max_mean_error<0||o.handoff_max_mean_error>1)return PB11_STATUS_OUT_OF_RANGE;
 if(b.relative_order<4||b.relative_order>64||b.cm_order<4||b.cm_order>32||b.nq<4||b.nq>1024||b.ncos<4||b.ncos>1024||b.remainder_policy<0||b.remainder_policy>2||(b.broad_mode!=1&&b.broad_mode!=3&&b.broad_mode!=13)||b.fsci_policy<0||b.fsci_policy>1||o.handoff_enabled<0||o.handoff_enabled>1)return BAD;
 for(int x:o.channels)if(x!=0&&x!=1)return BAD;
 double sigma=0;int st=fusion_c_cross_section_model(0,b.continuation,b.pb_low,MeV,&sigma);if(st)return st;
 fusion_nuclear_channel_v1 pb{};st=fusion_c_nuclear_channel(0,&pb);if(st)return st;
 return b.ground_state_q_J<=pb.q_J?PB11_STATUS_OK:PB11_STATUS_OUT_OF_RANGE;
}
}
namespace {
int coupled_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_birth_table_v1*const*tables,bool table_mode,int effective_charge,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out,
 fusion_handoff_diagnostics_v1*diagnostics=nullptr,bool require_diagnostics=false,const fusion_fast_target_options_v1*fastop=nullptr,bool require_fast=false,
 int beam_count=0,const fusion_beam_table_entry_v1*beam_entries=nullptr,
 fusion_beam_table_usage_v1*usage=nullptr,bool require_usage=false,
 const fusion_coupled_floor_limits_v1*floor_limits=nullptr,
 fusion_birth_floor_ledger_v1*floor_ledger=nullptr,bool require_floor=false,
 int table_policy=0,uint64_t*outside_direct=nullptr,bool require_policy=false,
 double*packet_output=nullptr,fusion_birth_packets_v1*packet_meta=nullptr,bool require_packets=false){
 trial_stage("begin");
 if(packet_meta)*packet_meta={};
 if(packet_output&&n>0&&n<=100000)std::fill(packet_output,packet_output+70*n,0.);
 if(outside_direct)*outside_direct=0;
 uint64_t outside_count=0;
 if(floor_ledger)*floor_ledger={};
 fusion_birth_floor_ledger_v1 floor_result{};
 if(usage)*usage={};
 fusion_beam_table_usage_v1 source_usage{};
 if(out)*out={};
 if(diagnostics)*diagnostics={};
 fusion_handoff_diagnostics_v1 observation{};
 if(new_thermal)std::fill(new_thermal,new_thermal+6,0.);
 if(n>0&&n<=100000){if(new_s)std::fill(new_s,new_s+6*n,0.);if(new_t)std::fill(new_t,new_t+6*n,0.);}
 if((require_packets&&(!packet_output||!packet_meta))||(require_policy&&!outside_direct)||(require_floor&&!floor_ledger)||!out||!new_thermal||!new_s||!new_t||(require_diagnostics&&!diagnostics)||(require_usage&&!usage))return PB11_STATUS_NULL_OUTPUT;
 if(table_policy!=FUSION_BEAM_TABLE_STRICT&&table_policy!=FUSION_BEAM_TABLE_DIRECT_OUTSIDE)return BAD;
 if(!op||!edges||!thermal||!charge2||!logs||!old_s||!old_t||!external||!escape||n<1||n>100000||ninert<0||ninert>32||(ninert&&!inert)||(table_mode&&!tables))return BAD;
 if((effective_charge!=0&&effective_charge!=1)||beam_count<0||beam_count>10*n||
    (beam_count&&(!beam_entries||!fastop)))return BAD;
 if(!finite_value(dt)||!finite_value(Ue)||!finite_value(Ui)||!finite_value(ne))return BAD;
 if(dt<=0||Ue<=0||Ui<=0||ne<=0)return PB11_STATUS_OUT_OF_RANGE;
 int st=validate_options(*op);if(st)return st;
 if(require_fast&&!fastop)return BAD;
 if(require_floor&&!floor_limits)return BAD;
 if(floor_limits&&(!finite_value(floor_limits->max_center_over_ion_kT)||
    !finite_value(floor_limits->max_ion_energy_fraction)))return BAD;
 if(floor_limits&&(floor_limits->max_center_over_ion_kT<0||floor_limits->max_center_over_ion_kT>1||
    floor_limits->max_ion_energy_fraction<0||floor_limits->max_ion_energy_fraction>1))return PB11_STATUS_OUT_OF_RANGE;
 if(fastop){if(!finite_value(fastop->angular_max_exponent)||fastop->angular_max_exponent<8||fastop->angular_max_exponent>80||fastop->angular_order<4||fastop->angular_order>32)return BAD;for(int ch=0;ch<5;++ch)if(fastop->channels[ch]!=0&&fastop->channels[ch]!=1)return BAD;}
 const int nb=7+ninert;
 for(int j=0;j<=n;++j)if(!nonnegative(edges[j])||(j&&edges[j]<=edges[j-1]))return BAD;
 for(int i=0;i<6*n;++i)if(!nonnegative(old_s[i])||!nonnegative(old_t[i])||!nonnegative(external[i])||!nonnegative(escape[i]))return BAD;
 for(int i=0;i<6*nb;++i)if(!finite_value(logs[i])||logs[i]<=0)return PB11_STATUS_OUT_OF_RANGE;
 try{
  std::vector<R> packet_values;
  std::array<R,140> packet_below{},packet_above{}; // flattened [N/E][source][channel][species]
  if(require_packets)packet_values.assign(70*n,0);
  std::array<fusion_nuclear_mass_v1,8> mass{};for(int i=0;i<8;++i){st=fusion_c_nuclear_mass(i,&mass[i]);if(st)return st;}
  R Npool=0,Ninert=0;for(int i=0;i<6;++i){if(!nonnegative(thermal[i])||!nonnegative(charge2[i])||(!effective_charge&&charge2[i]>mass[i].nuclear_charge*mass[i].nuclear_charge))return BAD;Npool+=thermal[i];}
  for(int j=0;j<ninert;++j){if(!nonnegative(inert[j].density_m3)||!finite_value(inert[j].mass_kg)||inert[j].mass_kg<=0||!nonnegative(inert[j].mean_charge_squared))return BAD;Ninert+=inert[j].density_m3;}
  Npool+=Ninert;if(Npool<=0)return PB11_STATUS_OUT_OF_RANGE;
  double Ti=0,Te=0;if(!put(R(Ui)/(1.5L*Npool),Ti)||!put(R(Ue)/(1.5L*ne),Te)||Ti<=0||Te<=0)return NUM;
  std::vector<R> centers(n);for(int j=0;j<n;++j)centers[j]=(R(edges[j])+edges[j+1])/2;
  std::array<R,6> oldN{},oldE{};for(int i=0;i<6;++i)for(int j=0;j<n;++j){R z=R(old_s[i*n+j])+old_t[i*n+j];oldN[i]+=z;oldE[i]+=z*centers[j];}
  std::array<fusion_nuclear_channel_v1,5> reactions{};
  std::array<fusion_birth_coefficients_v1,5> source{};
  std::array<std::vector<double>,5> grids;
  double rates[5]{},mean_a[5]{},mean_b[5]{},oldUi[6]{},trialNi[6]{},trialUi[6]{};
  fusion_coupled_thermal_v1 result{};auto&l=result.ledger;
  std::array<R,6> belowN{},belowE{};
  for(int i=0;i<6;++i)if(!put(R(1.5L)*Ti*thermal[i],oldUi[i]))return NUM;
  for(int ch=0;ch<5;++ch){st=fusion_c_nuclear_channel(ch,&reactions[ch]);if(st)return st;
   if(!op->channels[ch]||thermal[reactions[ch].reactant_ids[0]]==0||thermal[reactions[ch].reactant_ids[1]]==0)continue;
   auto options=op->birth;if(ch!=0)options.pb_low=FUSION_PB_LOW_TB;
   grids[ch].resize(7*n);
   double rate_error=0,debit_error=0;
   if(table_mode){
    if(!fusion_detail::birth_table_matches(tables[ch],ch,options,n,edges))return BAD;
    fusion_birth_table_info_v1 info{};st=fusion_c_birth_table_info(tables[ch],&info);if(st)return st;
    rate_error=info.max_sampled_direct_rate_discrepancy;debit_error=info.max_sampled_direct_debit_discrepancy;
    st=fusion_c_birth_table_evaluate(tables[ch],Ti,n,grids[ch].data(),&source[ch]);if(st)return st;
   }else{
    fusion_thermal_birth_v1 direct{};
    st=fusion_c_thermal_birth_grid(ch,Ti,&options,n,edges,grids[ch].data(),&direct);if(st)return st;
    rate_error=std::abs(direct.relative_rate_discrepancy);debit_error=direct.relative_reactant_energy_discrepancy;
    auto&c=source[ch];c.reactivity_m3_s=direct.reactivity_m3_s;
    for(int i=0;i<2;++i)c.reactant_energy_moment_J_m3_s[i]=direct.reactant_energy_moment_J_m3_s[i];
    for(int i=0;i<7;++i){c.below_number_m3_s[i]=direct.below_number_m3_s[i];c.below_energy_J_m3_s[i]=direct.below_energy_J_m3_s[i];c.above_number_m3_s[i]=direct.above_number_m3_s[i];c.above_energy_J_m3_s[i]=direct.above_energy_J_m3_s[i];}
   }
   const auto&r=source[ch];result.max_source_rate_discrepancy=std::max(result.max_source_rate_discrepancy,rate_error);
   result.max_source_debit_discrepancy=std::max(result.max_source_debit_discrepancy,debit_error);
   if(rate_error>op->max_source_rate_error||debit_error>op->max_source_debit_error)return NUM;
   for(int i=0;i<6;++i)if((!floor_limits&&r.below_number_m3_s[i]>0)||r.above_number_m3_s[i]>0||(floor_limits&&r.above_energy_J_m3_s[i]>0))return PB11_STATUS_OUT_OF_RANGE;
   rates[ch]=r.reactivity_m3_s;
   if(rates[ch]>0){mean_a[ch]=r.reactant_energy_moment_J_m3_s[0]/rates[ch];mean_b[ch]=r.reactant_energy_moment_J_m3_s[1]/rates[ch];}
  }
  trial_stage("thermal_sources_done");
  fusion_thermal_burn_v1 burn{};st=fusion_c_thermal_burn_trial(dt,thermal,oldUi,rates,mean_a,mean_b,trialNi,trialUi,&burn);if(st)return st;
  R ionU=Ui,electronU=Ue;Npool=Ninert;for(int i=0;i<6;++i){l.thermal_consumed_number_m3[i]=burn.reactant_removed_m3[i];l.thermal_consumed_energy_J_m3[i]=burn.reactant_removed_energy_J_m3[i];ionU-=burn.reactant_removed_energy_J_m3[i];Npool+=trialNi[i];}
  if(ionU<=0||Npool<=0||!put(ionU/(1.5L*Npool),Ti)||Ti<=0)return NUM;
  // Old kinetic components are consumed before nuclear/external births enter FP.
  const double*fp_s=old_s;const double*fp_t=old_t;
  std::vector<double> surviving_s,surviving_t;
  std::vector<R> fast_birth;R fast_neutronN=0,fast_neutronE=0;
  std::array<R,5> fast_events{};
  bool fast_enabled=false;if(fastop)for(int ch=0;ch<5;++ch)fast_enabled=fast_enabled||fastop->channels[ch];
  std::vector<const fusion_beam_birth_table_v1*> beam_lookup;
  if(beam_count){
   beam_lookup.assign(10*n,nullptr);
   for(int e=0;e<beam_count;++e){
    const auto&entry=beam_entries[e];
    if(entry.channel<0||entry.channel>=5||entry.projectile_slot<0||entry.projectile_slot>1||
       entry.energy_cell<0||entry.energy_cell>=n||!fastop->channels[entry.channel])return BAD;
    const auto&reaction=reactions[entry.channel];
    if(entry.projectile_slot==1&&reaction.reactant_ids[0]==reaction.reactant_ids[1])return BAD;
    auto expected=beam_options_for(op->birth,*fastop,entry.channel);
    st=fusion_c_beam_birth_resolve_support(double(centers[entry.energy_cell]),&expected,&expected);
    if(st)return st;
    if(!fusion_detail::beam_birth_table_matches(entry.table,entry.channel,entry.projectile_slot,
       double(centers[entry.energy_cell]),expected,n,edges))return BAD;
    auto&destination=beam_lookup[(entry.channel*2+entry.projectile_slot)*n+entry.energy_cell];
    if(destination)return BAD;destination=entry.table;
   }
  }
  if(fast_enabled){
   trial_stage("fast_begin");
   struct EdgeMeta{int channel,slot,index,cached;};
   struct CachedBeam{fusion_beam_birth_v1 value;std::vector<double> spectrum;};
   std::vector<CachedBeam> cache;
   // Exact intra-trial reuse only: neither backgrounds nor accepted steps share
   // cached values. Bound spectrum payload, and recompute after the cap or an
   // allocation failure. No interpolation, truncation or moment-only surrogate.
   constexpr std::size_t cache_payload_limit=32*1024*1024;
   const std::size_t spectrum_bytes=std::size_t(7)*n*sizeof(double);
   bool cache_allocation_ok=true;
   std::vector<fusion_target_network_edge_v1> links;std::vector<EdgeMeta> meta;
   std::vector<double> fast_initial(6*n),fast_energy(6*n),spectrum(7*n);
   for(int i=0;i<6;++i)for(int j=0;j<n;++j){int k=i*n+j;if(!put(R(old_s[k])+old_t[k],fast_initial[k]))return NUM;fast_energy[k]=double(centers[j]);}
   const double fastTi=Ti;
   std::array<R,6> spill_number_bound{},spill_energy_bound{};
   auto sample=[&](int ch,int slot,int k,fusion_beam_birth_v1&value,bool bound_spill)->int{
    auto base=beam_options_for(op->birth,*fastop,ch);
    fusion_beam_birth_options_v1 beam{};
    int support_status=fusion_c_beam_birth_resolve_support(fast_energy[k],&base,&beam);
    if(support_status)return support_status;
    const auto*table=beam_lookup.empty()?nullptr:beam_lookup[(ch*2+slot)*n+k%n];
    fusion_beam_birth_table_info_v1 info{};
    if(table){
     int code=fusion_c_beam_birth_table_info(table,&info);if(code)return code;
     if(table_policy==FUSION_BEAM_TABLE_DIRECT_OUTSIDE&&
        (fastTi<info.lower_kT_J||fastTi>info.upper_kT_J)){
      table=nullptr;++outside_count;
     }
    }
    if(table){
     fusion_birth_coefficients_v1 c{};
     int code=fusion_c_beam_birth_table_evaluate(table,fastTi,n,spectrum.data(),&c);if(code)return code;
     ++source_usage.table_evaluations;
     source_usage.max_validated_rate_error=std::max(source_usage.max_validated_rate_error,info.max_validated_rate_error);
     source_usage.max_validated_debit_error=std::max(source_usage.max_validated_debit_error,info.max_validated_debit_error);
     source_usage.max_validated_number_L1=std::max(source_usage.max_validated_number_L1,info.max_validated_number_L1);
     source_usage.max_validated_energy_L1=std::max(source_usage.max_validated_energy_L1,info.max_validated_energy_L1);
     auto&r=value.spectrum;r.reactivity_m3_s=c.reactivity_m3_s;
     for(int i=0;i<2;++i)r.reactant_energy_moment_J_m3_s[i]=c.reactant_energy_moment_J_m3_s[i];
     for(int i=0;i<7;++i){r.below_number_m3_s[i]=c.below_number_m3_s[i];r.below_energy_J_m3_s[i]=c.below_energy_J_m3_s[i];r.above_number_m3_s[i]=c.above_number_m3_s[i];r.above_energy_J_m3_s[i]=c.above_energy_J_m3_s[i];}
     r.relative_rate_discrepancy=info.max_sampled_direct_rate_discrepancy;
     r.relative_reactant_energy_discrepancy=info.max_sampled_direct_debit_discrepancy;
    }else{
     // The same per-projectile support is used by table matching and direct
     // fallback; it cannot accumulate across unrelated samples or channels.
     int code=fusion_c_beam_birth_grid(ch,slot,fast_energy[k],fastTi,&beam,n,edges,spectrum.data(),&value);if(code)return code;
     ++source_usage.direct_evaluations;
    }
    auto&r=value.spectrum;double rate_error=std::abs(r.relative_rate_discrepancy),debit_error=r.relative_reactant_energy_discrepancy;
    result.max_source_rate_discrepancy=std::max(result.max_source_rate_discrepancy,rate_error);result.max_source_debit_discrepancy=std::max(result.max_source_debit_discrepancy,debit_error);
    if(rate_error>op->max_source_rate_error||debit_error>op->max_source_debit_error)return NUM;
    if(bound_spill)for(int i=0;i<6;++i){
     if(floor_limits&&r.above_number_m3_s[i]==0&&r.above_energy_J_m3_s[i]==0)continue;
     if(r.below_number_m3_s[i]==0&&r.above_number_m3_s[i]==0&&
        r.below_energy_J_m3_s[i]==0&&r.above_energy_J_m3_s[i]==0)continue;
     // No edge can consume more than its entire accepted fast inventory.
     // Upper-neighbor coefficients also bound a subnormal moment rounded to
     // zero. Summing over competing channels overestimates the possible spill.
     double lower_rate=std::nextafter(r.reactivity_m3_s,0.);
     if(lower_rate<=0)return PB11_STATUS_OUT_OF_RANGE;
     auto upper=[](double x)->R{return std::nextafter(x,std::numeric_limits<double>::infinity());};
     R scale=R(fast_initial[k])/lower_rate;
     spill_number_bound[i]+=scale*((floor_limits?R(0):upper(r.below_number_m3_s[i]))+upper(r.above_number_m3_s[i]));
     spill_energy_bound[i]+=scale*((floor_limits?R(0):upper(r.below_energy_J_m3_s[i]))+upper(r.above_energy_J_m3_s[i]));
    }
    return PB11_STATUS_OK;
   };
   for(int ch=0;ch<5;++ch)if(fastop->channels[ch])for(int slot=0;slot<2;++slot){
    const auto&r=reactions[ch];if(slot==1&&r.reactant_ids[0]==r.reactant_ids[1])continue;
    int projectile=r.reactant_ids[slot],target=r.reactant_ids[1-slot];if(trialNi[target]==0)continue;
    for(int j=0;j<n;++j){int k=projectile*n+j;if(fast_initial[k]==0)continue;fusion_beam_birth_v1 value{};st=sample(ch,slot,k,value,true);if(st)return st;
     int cached=-1;
     if(cache_allocation_ok&&value.spectrum.reactivity_m3_s>0&&cache.size()<cache_payload_limit/spectrum_bytes){
      try{cache.push_back({value,spectrum});cached=int(cache.size()-1);}
      catch(const std::bad_alloc&){cache_allocation_ok=false;}
     }
     links.push_back({k,target,value.spectrum.reactivity_m3_s,value.spectrum.reactant_energy_moment_J_m3_s[1-slot]});meta.push_back({ch,slot,k,cached});
    }
   }
   trial_stage("fast_links_done");
   // Permit only spill whose conservative full-consumption bounds, including
   // step-equivalent source rates, cannot be represented by the public double
   // state. Factor two keeps margin for positive long-double arithmetic.
   // This is not an adjustable tail tolerance: every representable spill fails.
   R source_scale=std::max(R(1),R(1)/dt);
   for(int i=0;i<6;++i){
    R nb=2*source_scale*spill_number_bound[i],eb=2*source_scale*spill_energy_bound[i];
    if(!std::isfinite(nb)||!std::isfinite(eb)||nb>std::numeric_limits<double>::max()||eb>std::numeric_limits<double>::max()||
       double(nb)!=0||double(eb)!=0)return PB11_STATUS_OUT_OF_RANGE;
   }
   if(!links.empty()){
    double targetU[6]{},afterN[6]{},afterU[6]{};
    for(int i=0;i<6;++i)if(!put(1.5L*fastTi*trialNi[i],targetU[i]))return NUM;
    std::vector<double> afterFast(6*n),loss(links.size());fusion_target_network_v1 debit{};
    st=fusion_c_target_network_trial(6*n,6,int(links.size()),dt,fast_energy.data(),fast_initial.data(),trialNi,targetU,links.data(),afterFast.data(),afterN,afterU,loss.data(),&debit);if(st)return st;
    trial_stage("fast_network_done");
    // Use the same hazard denominator for S and T, avoiding loss of a small
    // component when a summed inventory rounds or the other component dominates.
    std::vector<R> hazard(6*n,0);for(size_t e=0;e<links.size();++e)hazard[links[e].fast_index]+=R(dt)*links[e].reactivity_m3_s*afterN[links[e].target_index];
    surviving_s.resize(6*n);surviving_t.resize(6*n);
    for(int k=0;k<6*n;++k){if(!put(R(old_s[k])/(1+hazard[k]),surviving_s[k])||!put(R(old_t[k])/(1+hazard[k]),surviving_t[k]))return NUM;
     R round=16*std::numeric_limits<double>::epsilon()*(R(old_s[k])+old_t[k]+afterFast[k]);
     if(!close({R(surviving_s[k]),R(surviving_t[k]),-R(afterFast[k])},round))return NUM;
    }
    fp_s=surviving_s.data();fp_t=surviving_t.data();fast_birth.assign(6*n,0);
    std::array<R,6> fastNremoved{},fastEremoved{},targetNremoved{},targetEremoved{};
    for(size_t e=0;e<links.size();++e){if(loss[e]==0)continue;const auto&m=meta[e];const auto&link=links[e];if(link.reactivity_m3_s<=0)return NUM;
     fusion_beam_birth_v1 value{};
     if(m.cached>=0){value=cache[m.cached].value;spectrum=cache[m.cached].spectrum;}
     else {st=sample(m.channel,m.slot,m.index,value,false);if(st)return st;}
     const auto&r=value.spectrum;if(r.reactivity_m3_s!=link.reactivity_m3_s||r.reactant_energy_moment_J_m3_s[1-m.slot]!=link.target_energy_reactivity_J_m3_s)return NUM;
     R amount=loss[e],scale=amount/link.reactivity_m3_s;int projectile=m.index/n;
     fast_events[m.channel]+=amount;fastNremoved[projectile]+=amount;fastEremoved[projectile]+=amount*centers[m.index%n];
     targetNremoved[link.target_index]+=amount;targetEremoved[link.target_index]+=scale*link.target_energy_reactivity_J_m3_s;
     for(int k=0;k<6*n;++k)fast_birth[k]+=scale*spectrum[k];
     if(require_packets){
      for(int k=0;k<7*n;++k)packet_values[(5+m.channel)*7*n+k]+=scale*spectrum[k];
      for(int i=0;i<7;++i){int k=(5+m.channel)*7+i;
       packet_below[k]+=scale*r.below_number_m3_s[i];packet_below[70+k]+=scale*r.below_energy_J_m3_s[i];
       packet_above[k]+=scale*r.above_number_m3_s[i];packet_above[70+k]+=scale*r.above_energy_J_m3_s[i];}
     }
     if(floor_limits)for(int i=0;i<6;++i){belowN[i]+=scale*r.below_number_m3_s[i];belowE[i]+=scale*r.below_energy_J_m3_s[i];}
     R neutronNumber=R(r.below_number_m3_s[6])+r.above_number_m3_s[6],neutronEnergy=R(r.below_energy_J_m3_s[6])+r.above_energy_J_m3_s[6];
     for(int j=0;j<n;++j){neutronNumber+=spectrum[6*n+j];neutronEnergy+=R(spectrum[6*n+j])*centers[j];}
     fast_neutronN+=scale*neutronNumber;fast_neutronE+=scale*neutronEnergy;
    }
    trial_stage("fast_losses_done");
    Npool=Ninert;
    for(int i=0;i<6;++i){trialNi[i]=afterN[i];Npool+=trialNi[i];ionU-=targetEremoved[i];
     if(!put(fastNremoved[i],l.fast_consumed_number_m3[i])||!put(fastEremoved[i],l.fast_consumed_energy_J_m3[i])||!put(R(l.thermal_consumed_number_m3[i])+targetNremoved[i],l.thermal_consumed_number_m3[i])||!put(R(l.thermal_consumed_energy_J_m3[i])+targetEremoved[i],l.thermal_consumed_energy_J_m3[i]))return NUM;
    }
    if(ionU<=0||Npool<=0||!put(ionU/(1.5L*Npool),Ti)||Ti<=0)return NUM;
   }
   trial_stage("fast_done");
  }
  if(floor_limits){
   for(int ch=0;ch<5;++ch)if(burn.events_m3[ch]>0){
    if(rates[ch]<=0)return NUM;
    R scale=R(burn.events_m3[ch])/rates[ch];
    for(int i=0;i<6;++i){belowN[i]+=scale*source[ch].below_number_m3_s[i];belowE[i]+=scale*source[ch].below_energy_J_m3_s[i];}
   }
   double N[6]{},E[6]{},reservoir=0;
   for(int i=0;i<6;++i){
    if(!put(belowN[i],N[i])||!put(belowE[i],E[i]))return NUM;
    // Do not lose a source rate representable in double when the corresponding
    // per-step amount rounds to zero. Such a case needs a rate-aware interface.
    if((N[i]==0&&double(belowN[i]/dt)!=0)||(E[i]==0&&double(belowE[i]/dt)!=0))return NUM;
    double rate=0;
    if(!floor_source_rate(N[i],dt,rate))return NUM;
    const R nr=std::abs(belowN[i]-N[i]),er=std::abs(belowE[i]-E[i]);
    const R eps=std::numeric_limits<double>::epsilon(),tiny=std::numeric_limits<double>::denorm_min();
    if(nr>eps*std::abs(belowN[i])+tiny||er>eps*std::abs(belowE[i])+tiny)return NUM;
   }
   if(!put(ionU,reservoir))return NUM;
   fusion_birth_floor_options_v1 options{double(centers[0]),Ti,reservoir,
     floor_limits->max_center_over_ion_kT,floor_limits->max_ion_energy_fraction};
   st=fusion_c_birth_floor_project(&options,N,E,&floor_result);if(st)return st;
   for(int i=0;i<6;++i){
    const R actual=R(N[i])*centers[0],mapped=floor_result.mapped_energy_J_m3[i];
    const R allowance=2*std::numeric_limits<double>::epsilon()*std::abs(actual)+
       R(std::numeric_limits<double>::denorm_min());
    if(!std::isfinite(actual)||std::abs(actual-mapped)>allowance)return NUM;
   }
   bool active=false;for(int i=0;i<6;++i){active=active||N[i]>0;ionU+=floor_result.ion_energy_correction_J_m3[i];}
   if(active&&(!put(ionU/(1.5L*Npool),Ti)||Ti<=0))return NUM;
   trial_stage("floor_done");
  }
  std::vector<double> birth(6*n),s(6*n),t(6*n);R unrepresentedN=0,unrepresentedE=0;R Q=0,neutronN=fast_neutronN,neutronE=fast_neutronE;
  for(int ch=0;ch<5;++ch){if(!put(R(burn.events_m3[ch])+fast_events[ch],l.events_m3[ch]))return NUM;Q+=R(l.events_m3[ch])*reactions[ch].q_J;}
  for(int i=0;i<6;++i){R bornN=0,bornE=0,extN=0,extE=0;
   for(int j=0;j<n;++j){R amount=fast_birth.empty()?0:fast_birth[i*n+j];
    for(int ch=0;ch<5;++ch)if(burn.events_m3[ch]>0){if(rates[ch]<=0)return NUM;amount+=R(grids[ch][i*n+j])/rates[ch]*burn.events_m3[ch];}
    bornN+=amount;bornE+=amount*centers[j];
    if(floor_limits&&j==0)amount+=floor_result.mapped_number_m3[i];
    R ex=R(dt)*external[i*n+j];extN+=ex;extE+=ex*centers[j];
    if(floor_limits){
     R lost=0;if(!floor_source_rate(amount+ex,dt,birth[i*n+j],&lost))return NUM;
     if(!source_roundoff_accumulate(lost,centers[j],unrepresentedN,unrepresentedE))return NUM;
    }
    else if(!put((amount+ex)/dt,birth[i*n+j]))return NUM;
   }
   if(floor_limits){bornN+=floor_result.born_number_m3[i];bornE+=floor_result.born_energy_J_m3[i];}
   if(!put(bornN,l.nuclear_born_number_m3[i])||!put(bornE,l.nuclear_born_energy_J_m3[i])||!put(extN,l.external_born_number_m3[i])||!put(extE,l.external_born_energy_J_m3[i]))return NUM;
  }
  for(int ch=0;ch<5;++ch)if(burn.events_m3[ch]>0){R N=R(source[ch].below_number_m3_s[6])+source[ch].above_number_m3_s[6],E=R(source[ch].below_energy_J_m3_s[6])+source[ch].above_energy_J_m3_s[6];for(int j=0;j<n;++j){N+=grids[ch][6*n+j];E+=R(grids[ch][6*n+j])*centers[j];}neutronN+=N/rates[ch]*burn.events_m3[ch];neutronE+=E/rates[ch]*burn.events_m3[ch];}
  if(!put(neutronN,l.neutron_number_m3)||!put(neutronE,l.neutron_energy_J_m3))return NUM;
  std::vector<fusion_maxwellian_bath_v1> baths(nb);
  baths[0]={ne,mass[7].mass_kg,1,Te,1};
  for(int j=0;j<6;++j)baths[j+1]={trialNi[j],mass[j].mass_kg,charge2[j],Ti,1};
  for(int j=0;j<ninert;++j)baths[j+7]={inert[j].density_m3,inert[j].mass_kg,inert[j].mean_charge_squared,Ti,1};
  std::vector<double> temperatures(nb),diffusion(nb*(n-1)),transfer(n),zero(n),heat(nb),after(n);
  for(int b=0;b<nb;++b)temperatures[b]=baths[b].kT_J;
  for(int i=0;i<6;++i){bool active=oldN[i]>0||l.nuclear_born_number_m3[i]>0||l.external_born_number_m3[i]>0;if(!active)continue;
   std::fill(transfer.begin(),transfer.end(),0.);
   for(int b=0;b<nb;++b){baths[b].coulomb_log=logs[i*nb+b];
    for(int j=0;j<n-1;++j){fusion_coulomb_energy_v1 coefficient{};st=fusion_c_coulomb_energy(edges[j+1],mass[i].mass_kg,mass[i].nuclear_charge,&baths[b],&coefficient);if(st)return st;diffusion[b*(n-1)+j]=coefficient.diffusion_J2_s;}
    if(b>0)for(int j=0;j<n;++j){double lambda=0;st=fusion_c_coulomb_transfer_rate(double(centers[j]),mass[i].mass_kg,mass[i].nuclear_charge,&baths[b],&lambda);if(st)return st;if(!put(R(transfer[j])+lambda,transfer[j]))return NUM;}
   }
   fusion_two_component_ledger_v1 fp{};
   st=fusion_c_two_component_trial(n,nb,dt,edges,fp_s+i*n,fp_t+i*n,temperatures.data(),diffusion.data(),birth.data()+i*n,zero.data(),escape+i*n,transfer.data(),s.data()+i*n,t.data()+i*n,heat.data(),&fp);if(st)return st;
   if(diagnostics){
    observation.transferred_number_m3[i]=fp.transferred_number_m3;
    observation.transferred_energy_J_m3[i]=fp.transferred_energy_J_m3;
    R candidateN=0,candidateE=0;
    for(int j=0;j<n;++j){candidateN+=t[i*n+j];candidateE+=R(t[i*n+j])*centers[j];}
    if(!put(candidateN,observation.candidate_number_m3[i])||
       !put(candidateE,observation.candidate_energy_J_m3[i]))return NUM;
   }
   l.escaped_number_m3[i]=fp.total.escaped_number_m3;l.escaped_energy_J_m3[i]=fp.total.escaped_energy_J_m3;
   electronU+=heat[0];R extra=0;for(int b=0;b<nb;++b){if(b<7)l.heat_to_bath_J_m3[i*7+b]=heat[b];else extra+=heat[b];if(b>0)ionU+=heat[b];}
   if(!put(extra,result.inert_ion_heat_J_m3[i]))return NUM;
  }
  trial_stage("two_component_done");
  if(electronU<=0||ionU<=0)return NUM;
  if(op->handoff_enabled)for(int i=0;i<6;++i){R N=0,E=0;for(int j=0;j<n;++j){N+=t[i*n+j];E+=R(t[i*n+j])*centers[j];}if(N==0)continue;
   double mixedTi=0;if(!put((ionU+E)/(1.5L*(Npool+N)),mixedTi)||mixedTi<=0)return NUM;
   if(diagnostics){observation.target_kT_J[i]=mixedTi;observation.tested[i]=1;}
   fusion_handoff_ledger_v1 handoff{};int projected=0;
   st=fusion_c_maxwellian_handoff_trial(n,mixedTi,op->handoff_max_L1,op->handoff_max_mean_error,edges,t.data()+i*n,after.data(),&projected,&handoff);if(st)return st;
   result.handoff_L1[i]=handoff.distribution_L1;result.handoff_mean_error[i]=handoff.relative_mean_energy_error;result.handoff_projected[i]=projected;
   if(projected){l.handed_off_number_m3[i]=handoff.fluid_number_m3;l.handed_off_energy_J_m3[i]=handoff.fluid_energy_J_m3;
    if(!put(R(l.heat_to_bath_J_m3[i*7+i+1])+handoff.bath_energy_correction_J_m3,l.heat_to_bath_J_m3[i*7+i+1])||!put(R(trialNi[i])+handoff.fluid_number_m3,trialNi[i]))return NUM;
    ionU+=R(handoff.fluid_energy_J_m3)+handoff.bath_energy_correction_J_m3;Npool+=handoff.fluid_number_m3;
    std::copy(after.begin(),after.end(),t.begin()+i*n);
   }
  }
  trial_stage("handoff_done");
  if(!put(electronU,result.electron_energy_J_m3)||!put(ionU,result.ion_energy_J_m3))return NUM;
  // Independent accounting of reaction, kinetic, and complete local inventories.
  R nuclearE=neutronE,removedE=0,totalOld=R(Ue)+Ui,totalNew=R(result.electron_energy_J_m3)+result.ion_energy_J_m3;
  R escapedE=0,externalE=0;
  for(int i=0;i<6;++i){R N=0,E=0,H=result.inert_ion_heat_J_m3[i];for(int b=0;b<7;++b)H+=l.heat_to_bath_J_m3[i*7+b];
   for(int j=0;j<n;++j){R z=R(s[i*n+j])+t[i*n+j];N+=z;E+=z*centers[j];}
   if(!close({N,-oldN[i],-R(l.nuclear_born_number_m3[i]),-R(l.external_born_number_m3[i]),R(l.escaped_number_m3[i]),R(l.handed_off_number_m3[i]),R(l.fast_consumed_number_m3[i])})||
      !close({E,-oldE[i],-R(l.nuclear_born_energy_J_m3[i]),-R(l.external_born_energy_J_m3[i]),R(l.escaped_energy_J_m3[i]),R(l.handed_off_energy_J_m3[i]),H,R(l.fast_consumed_energy_J_m3[i]),R(floor_result.ion_energy_correction_J_m3[i])}))return NUM;
   R dN=R(trialNi[i])-thermal[i]+N-oldN[i]-l.nuclear_born_number_m3[i]+l.thermal_consumed_number_m3[i]-l.external_born_number_m3[i]+l.escaped_number_m3[i]+l.fast_consumed_number_m3[i];
   // Difference-of-inventory residuals inherit rounding from BOTH pools,
   // including a seeded kinetic species with no thermal counterpart.
   R round=16*std::numeric_limits<double>::epsilon()*(R(thermal[i])+trialNi[i]+oldN[i]+N);
   if(!put(dN,result.particle_residual_m3[i])||!close({R(trialNi[i])-thermal[i],N-oldN[i],-R(l.nuclear_born_number_m3[i]),R(l.thermal_consumed_number_m3[i]),-R(l.external_born_number_m3[i]),R(l.escaped_number_m3[i]),R(l.fast_consumed_number_m3[i])},round))return NUM;
   nuclearE+=l.nuclear_born_energy_J_m3[i];removedE+=R(l.thermal_consumed_energy_J_m3[i])+l.fast_consumed_energy_J_m3[i];totalOld+=oldE[i];totalNew+=E;escapedE+=l.escaped_energy_J_m3[i];externalE+=l.external_born_energy_J_m3[i];
  }
  if(!close({nuclearE,-removedE,-Q}))return NUM;
  R residual=totalNew-totalOld+neutronE+escapedE-externalE-Q;
  R round=16*std::numeric_limits<double>::epsilon()*(std::abs(totalOld)+std::abs(totalNew));
  if(!put(residual,result.energy_residual_J_m3)||!close({totalNew-totalOld,neutronE,escapedE,-externalE,-Q},round))return NUM;
  bool valid=true;ledger_fields(l,[&](double x){if(!finite_value(x))valid=false;});if(!valid)return NUM;
  std::vector<double> packet_rounded;fusion_birth_packets_v1 packet_result{};
  if(require_packets){
   auto packet_put=[](R value,double&out){return value>=0&&put(value,out)&&(value==0||out!=0);};
   for(int ch=0;ch<5;++ch){
    if(!packet_put(burn.events_m3[ch],packet_result.events_m3[0][ch])||!packet_put(fast_events[ch],packet_result.events_m3[1][ch]))return NUM;
    if(burn.events_m3[ch]>0){
     if(rates[ch]<=0)return NUM;
     for(int k=0;k<7*n;++k)packet_values[ch*7*n+k]=R(grids[ch][k])/rates[ch]*burn.events_m3[ch];
     R scale=R(burn.events_m3[ch])/rates[ch];
     for(int i=0;i<7;++i){int k=ch*7+i;
      packet_below[k]=scale*source[ch].below_number_m3_s[i];packet_below[70+k]=scale*source[ch].below_energy_J_m3_s[i];
      packet_above[k]=scale*source[ch].above_number_m3_s[i];packet_above[70+k]=scale*source[ch].above_energy_J_m3_s[i];}
    }
   }
   packet_rounded.resize(70*n);
   for(int k=0;k<70*n;++k)if(!packet_put(packet_values[k],packet_rounded[k]))return NUM;
   for(int src=0;src<2;++src)for(int ch=0;ch<5;++ch)for(int i=0;i<7;++i){int k=(src*5+ch)*7+i;
    if(!packet_put(packet_below[k],packet_result.below_number_m3[src][ch][i])||!packet_put(packet_below[70+k],packet_result.below_energy_J_m3[src][ch][i])||
       !packet_put(packet_above[k],packet_result.above_number_m3[src][ch][i])||!packet_put(packet_above[70+k],packet_result.above_energy_J_m3[src][ch][i]))return NUM;
   }
   // Independently compare source packets against the production physical ledger.
   for(int i=0;i<7;++i){R N=0,E=0;
    for(int src=0;src<2;++src)for(int ch=0;ch<5;++ch){int k=(src*5+ch)*7+i;
     N+=packet_below[k]+packet_above[k];E+=packet_below[70+k]+packet_above[70+k];
     for(int j=0;j<n;++j){N+=packet_values[k*n+j];E+=packet_values[k*n+j]*centers[j];}
    }
    if(!close({N,-R(i<6?l.nuclear_born_number_m3[i]:l.neutron_number_m3)})||!close({E,-R(i<6?l.nuclear_born_energy_J_m3[i]:l.neutron_energy_J_m3)}))return NUM;
   }
   std::copy(packet_rounded.begin(),packet_rounded.end(),packet_output);*packet_meta=packet_result;
  }
  trial_stage("complete");
  std::copy(trialNi,trialNi+6,new_thermal);std::copy(s.begin(),s.end(),new_s);std::copy(t.begin(),t.end(),new_t);*out=result;if(diagnostics)*diagnostics=observation;if(usage)*usage=source_usage;if(floor_ledger)*floor_ledger=floor_result;if(outside_direct)*outside_direct=outside_count;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}

}
extern "C" int fusion_c_coupled_thermal_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out){
 return coupled_trial(dt,op,nullptr,false,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out);
}
extern "C" int fusion_c_coupled_thermal_table_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out){
 return coupled_trial(dt,op,tables,true,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out);
}

extern "C" int fusion_c_coupled_thermal_table_trial_effective_charge(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out){
 return coupled_trial(dt,op,tables,true,true,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out);
}

extern "C" int fusion_c_coupled_thermal_trial_diagnosed(double dt,const fusion_coupled_thermal_options_v1*op,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics){
 return coupled_trial(dt,op,nullptr,false,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true);
}
extern "C" int fusion_c_coupled_thermal_table_trial_diagnosed(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics){
 return coupled_trial(dt,op,tables,true,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true);
}

extern "C" int fusion_c_coupled_thermal_table_trial_effective_charge_diagnosed(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics){
 return coupled_trial(dt,op,tables,true,true,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true);
}

static int thermal_increment(const fusion_source_ledger_v1*ledger,
 const double*inert,fusion_thermal_increment_v1*out,const double*numerical=nullptr,bool require_numerical=false){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(!ledger||!inert)return BAD;
 if(require_numerical&&!numerical)return BAD;
 if(numerical)for(int i=0;i<6;++i)if(!finite_value(numerical[i]))return BAD;
 auto copy=*ledger;
 // Check all ledger fields, treating only the documented heat fields as signed.
 for(double&h:copy.heat_to_bath_J_m3){if(!finite_value(h))return BAD;h=0;}
 bool valid=true;ledger_fields(copy,[&](double v){if(!nonnegative(v))valid=false;});
 if(!valid)return BAD;
 for(int i=0;i<6;++i)if(!finite_value(inert[i]))return BAD;
 fusion_thermal_increment_v1 result{};R electron=0,ion=0;
 for(int i=0;i<6;++i){
  if(!put(R(ledger->handed_off_number_m3[i])-ledger->thermal_consumed_number_m3[i],result.thermal_number_m3[i]))return NUM;
  electron+=ledger->heat_to_bath_J_m3[i*7];
  ion+=R(ledger->handed_off_energy_J_m3[i])-ledger->thermal_consumed_energy_J_m3[i];
  ion+=inert[i];
  for(int j=1;j<7;++j)ion+=ledger->heat_to_bath_J_m3[i*7+j];
 }
 if(numerical)for(int i=0;i<6;++i)ion+=numerical[i];
 if(!put(electron,result.electron_energy_J_m3)||!put(ion,result.ion_energy_J_m3))return NUM;
 *out=result;return PB11_STATUS_OK;
}

extern "C" int fusion_c_coupled_fast_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out){
 return coupled_trial(dt,op,nullptr,false,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,nullptr,false,fastop,true);
}

extern "C" int fusion_c_coupled_fast_table_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out){
 return coupled_trial(dt,op,tables,true,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,nullptr,false,fastop,true);
}

extern "C" int fusion_c_coupled_fast_table_trial_effective_charge(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out){
 return coupled_trial(dt,op,tables,true,true,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,nullptr,false,fastop,true);
}

extern "C" int fusion_c_coupled_fast_trial_diagnosed(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics){
 return coupled_trial(dt,op,nullptr,false,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true,fastop,true);
}

extern "C" int fusion_c_coupled_fast_table_trial_diagnosed(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics){
 return coupled_trial(dt,op,tables,true,false,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true,fastop,true);
}

extern "C" int fusion_c_coupled_fast_table_trial_effective_charge_diagnosed(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,
 const fusion_birth_table_v1*const*tables,int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics){
 return coupled_trial(dt,op,tables,true,true,n,edges,thermal,Ue,Ui,ne,charge2,ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true,fastop,true);
}

extern "C" int fusion_c_coupled_sources_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,const fusion_birth_table_v1*const*tables,
 int beam_count,const fusion_beam_table_entry_v1*beam_entries,int effective_charge,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,
 fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics,fusion_beam_table_usage_v1*usage){
 return coupled_trial(dt,op,tables,tables!=nullptr,effective_charge,n,edges,thermal,Ue,Ui,ne,charge2,
  ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true,fastop,true,
  beam_count,beam_entries,usage,true);
}

extern "C" int fusion_c_coupled_sources_floor_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,const fusion_birth_table_v1*const*tables,
 int beam_count,const fusion_beam_table_entry_v1*beam_entries,int effective_charge,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,
 fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics,fusion_beam_table_usage_v1*usage,
 const fusion_coupled_floor_limits_v1*floor_limits,fusion_birth_floor_ledger_v1*floor_ledger){
 return coupled_trial(dt,op,tables,tables!=nullptr,effective_charge,n,edges,thermal,Ue,Ui,ne,charge2,
  ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true,fastop,true,
  beam_count,beam_entries,usage,true,floor_limits,floor_ledger,true);
}

extern "C" int fusion_c_coupled_sources_covered_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,const fusion_birth_table_v1*const*tables,
 int beam_count,const fusion_beam_table_entry_v1*beam_entries,int effective_charge,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,
 fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics,fusion_beam_table_usage_v1*usage,
 const fusion_coupled_floor_limits_v1*floor_limits,fusion_birth_floor_ledger_v1*floor_ledger,int table_policy,uint64_t*outside_direct){
 return coupled_trial(dt,op,tables,tables!=nullptr,effective_charge,n,edges,thermal,Ue,Ui,ne,charge2,
  ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true,fastop,true,
  beam_count,beam_entries,usage,true,floor_limits,floor_ledger,
  floor_limits!=nullptr||floor_ledger!=nullptr,table_policy,outside_direct,true);
}

extern "C" int fusion_c_coupled_sources_packets_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 const fusion_fast_target_options_v1*fastop,const fusion_birth_table_v1*const*tables,
 int beam_count,const fusion_beam_table_entry_v1*beam_entries,int effective_charge,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,
 fusion_coupled_thermal_v1*out,fusion_handoff_diagnostics_v1*diagnostics,fusion_beam_table_usage_v1*usage,
 const fusion_coupled_floor_limits_v1*floor_limits,fusion_birth_floor_ledger_v1*floor_ledger,int table_policy,uint64_t*outside_direct,double*packet_output,fusion_birth_packets_v1*packet_meta){
 return coupled_trial(dt,op,tables,tables!=nullptr,effective_charge,n,edges,thermal,Ue,Ui,ne,charge2,
  ninert,inert,logs,old_s,old_t,external,escape,new_thermal,new_s,new_t,out,diagnostics,true,fastop,true,
  beam_count,beam_entries,usage,true,floor_limits,floor_ledger,
  floor_limits!=nullptr||floor_ledger!=nullptr,table_policy,outside_direct,true,packet_output,packet_meta,true);
}

extern "C" int fusion_c_coupled_thermal_increment(const fusion_source_ledger_v1*ledger,
 const double*inert,fusion_thermal_increment_v1*out){return thermal_increment(ledger,inert,out);}
extern "C" int fusion_c_coupled_numerical_increment(const fusion_source_ledger_v1*ledger,
 const double*inert,const double*numerical,fusion_thermal_increment_v1*out){
 return thermal_increment(ledger,inert,out,numerical,true);
}
