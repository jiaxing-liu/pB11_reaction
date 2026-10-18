// Covered source accuracy contracts; finite-source fixture only.
// Small finite-source API fixture; no physics/convergence claim.
#include "fusion_coupled_sources_budget.h"
#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
namespace {
constexpr double kev=1.602176634e-16, mev=1.602176634e-13;
constexpr int n=48;
void need(bool b,const char*m){if(!b)throw std::runtime_error(m);}
template<class T>bool same(const T&a,const T&b){return std::memcmp(&a,&b,sizeof(T))==0;}
template<class T>bool zero(const T&a){const auto*p=reinterpret_cast<const unsigned char*>(&a);return std::all_of(p,p+sizeof(T),[](unsigned char x){return x==0;});}
struct Input {
 std::array<double,n+1> edges{};
 std::array<double,6> thermal{{5e19,0,0,0,0,5e19}},charge{{1,1,1,4,4,25}};
 std::array<double,42> logs{};
 std::array<double,6*n> s{},t{},external{},escape{};
 fusion_coupled_thermal_options_v1 options{};
 fusion_fast_target_options_v1 fast{};
 double temperature=.5*kev,ne=3e20;
 Input(){
  edges[0]=0;for(int j=1;j<=n;++j)edges[j]=2e-3*kev*std::pow(25*mev/(2e-3*kev),double(j-1)/(n-1));edges[n]=25*mev;
  logs.fill(15);auto&o=options.birth;o.relative_max_J=2*mev;o.cm_max_kT=40;o.ground_state_q_J=91.84*kev;o.cutoff_J=.001*mev;o.l1_fraction=.76;o.narrow_peak_fraction=.051;o.continuum_peak_scale=1;o.continuation=FUSION_ENDPOINT_S;o.pb_low=FUSION_PB_LOW_TB;o.broad_mode=13;o.relative_order=o.cm_order=4;o.nq=o.ncos=8;
  options.channels[0]=1;options.max_source_rate_error=options.max_source_debit_error=.05;options.handoff_max_L1=options.handoff_max_mean_error=.001;
  fast.angular_order=8;fast.angular_max_exponent=40;
 }
 double ui()const{long double sum=0;for(double x:thermal)sum+=x;return double(1.5L*sum*temperature);}
};
struct Output {
 std::array<double,6> thermal{};
 std::array<double,6*n> s{},t{};
 fusion_coupled_thermal_v1 result{};
 fusion_handoff_diagnostics_v1 diagnostics{};
 fusion_beam_table_usage_v1 usage{};
 fusion_birth_floor_ledger_v1 floor{};
 uint64_t outside=0;
 fusion_coupled_source_accuracy_v1 accuracy{};
};
int call(const Input&i,Output&o,const fusion_thermal_birth_budget_policy_v1*p,
 const fusion_birth_table_v1*const*tables=nullptr,bool legacy=false,bool nullaccuracy=false,int coverage=FUSION_BEAM_TABLE_STRICT,double floor_center=.5){
 fusion_coupled_floor_limits_v1 limits{floor_center,.5};
 // Exercise all covered observations, including the optional floor ledger.
 if(legacy)return fusion_c_coupled_sources_covered_trial(1e-8,&i.options,&i.fast,tables,0,nullptr,0,n,i.edges.data(),i.thermal.data(),1.5*i.ne*i.temperature,i.ui(),i.ne,i.charge.data(),0,nullptr,i.logs.data(),i.s.data(),i.t.data(),i.external.data(),i.escape.data(),o.thermal.data(),o.s.data(),o.t.data(),&o.result,&o.diagnostics,&o.usage,&limits,&o.floor,coverage,&o.outside);
 return fusion_c_coupled_sources_covered_budget_trial(1e-8,&i.options,&i.fast,tables,0,nullptr,0,n,i.edges.data(),i.thermal.data(),1.5*i.ne*i.temperature,i.ui(),i.ne,i.charge.data(),0,nullptr,i.logs.data(),i.s.data(),i.t.data(),i.external.data(),i.escape.data(),o.thermal.data(),o.s.data(),o.t.data(),&o.result,&o.diagnostics,&o.usage,&limits,&o.floor,coverage,&o.outside,p,nullaccuracy?nullptr:&o.accuracy);
}
void physical_same(const Output&a,const Output&b){
 need(a.thermal==b.thermal&&a.s==b.s&&a.t==b.t,"physical arrays parity");
 need(same(a.result,b.result)&&same(a.diagnostics,b.diagnostics)&&same(a.usage,b.usage)&&same(a.floor,b.floor)&&a.outside==b.outside,"all physical observations parity");
}
void rejected(const Input&i,const fusion_thermal_birth_budget_policy_v1*p,const fusion_birth_table_v1*const*tables=nullptr,bool legacy=false,bool nullaccuracy=false,int coverage=FUSION_BEAM_TABLE_STRICT,double floor_center=.5){
 Output out;std::memset(&out,0x5a,sizeof(out));
 need(call(i,out,p,tables,legacy,nullaccuracy,coverage,floor_center)!=0,"invalid trial accepted");
 if(legacy||nullaccuracy)out.accuracy={};
 need(zero(out.thermal)&&zero(out.s)&&zero(out.t)&&zero(out.result)&&zero(out.diagnostics)&&zero(out.usage)&&zero(out.floor)&&zero(out.accuracy)&&out.outside==0,"failure did not clear every output");
}
void direct_contract(){
 Input i;auto original=i;fusion_thermal_birth_budget_policy_v1 p{1,4,4,0,0,0};Output full,exact,mixed;
 need(call(i,full,&p,nullptr,true)==0,"full physical fixture");need(call(i,exact,&p)==0,"zero-budget physical fixture");physical_same(full,exact);
 need(exact.accuracy.version==1&&exact.accuracy.mode==1&&zero(exact.accuracy.table),"direct mode metadata");
 need(exact.accuracy.direct.low_pairs==0&&exact.accuracy.direct.high_pairs>0,"zero budget full parents");
 p.number_L1_tolerance=p.energy_L1_tolerance=2;
 need(call(i,mixed,&p)==0,"nonzero budget physical fixture");auto&q=mixed.accuracy.direct;
 need(q.version==1&&q.low_pairs>0&&q.high_pairs==0&&q.mixed==1&&q.joint_nodes==q.low_pairs,"actual low parent use");
 need(mixed.accuracy.number_L1_estimate>=q.number_L1_bound+q.number_roundoff_allowance&&mixed.accuracy.energy_L1_estimate>=q.energy_L1_bound+q.energy_roundoff_allowance,"actual report plus allowance");
 need(full.result.ledger.events_m3[0]==mixed.result.ledger.events_m3[0],"unchanged burn rate");
 need(i.s==original.s&&i.t==original.t&&i.thermal==original.thermal,"inputs preserved");
 rejected(i,nullptr);rejected(i,&p,nullptr,false,true);rejected(i,&p,nullptr,false,false,77);rejected(i,&p,nullptr,false,false,FUSION_BEAM_TABLE_STRICT,1.01);
 Output coverage;need(call(i,coverage,&p,nullptr,false,false,FUSION_BEAM_TABLE_DIRECT_OUTSIDE)==0,"empty beam direct-outside policy");physical_same(mixed,coverage);need(coverage.outside==0,"no invented outside-domain calls");
 // Successful source above establishes this same fixture evaluates pB. A zero
 // direct discrepancy gate deliberately fails after that source has returned.
 auto late=i;late.options.max_source_rate_error=late.options.max_source_debit_error=0;rejected(late,&p);
 for(int absent=0;absent<2;++absent){auto inactive=i;if(absent)inactive.thermal[5]=0;else inactive.options.channels[0]=0;Output x;need(call(inactive,x,&p)==0,"inactive valid policy");need(x.accuracy.version==1&&x.accuracy.mode==0&&zero(x.accuracy.direct)&&zero(x.accuracy.table)&&x.accuracy.number_L1_estimate==0&&x.accuracy.energy_L1_estimate==0,"inactive no fictitious report");auto bad=p;bad.version=2;rejected(inactive,&bad);}
 // Other channel retains full direct arithmetic even with a generous pB policy.
 i.options.channels[0]=0;i.options.channels[3]=1;i.thermal={{0,5e19,5e19,0,0,0}};i.temperature=10*kev;i.fast.channels[3]=1;i.s[1*n+n/2]=1e12;i.t[1*n+n/2]=1e12;Output dtfull,dtbudget;need(call(i,dtfull,&p,nullptr,true)==0&&call(i,dtbudget,&p)==0,"DT fixture");physical_same(dtfull,dtbudget);need(dtbudget.accuracy.mode==0&&dtbudget.result.ledger.events_m3[3]>0,"DT fast source remains full without activating pB report");
}
void table_contract(){
 Input i;fusion_thermal_birth_budget_policy_v1 p{1,4,4,0,1e-4,1e-4};
 fusion_birth_table_control_v1 control{.01,.01,.01,.01,.05,.05,8,64,3};
 double lo=i.temperature,hi=lo*1.000001;fusion_birth_table_v1*raw=nullptr;
 need(fusion_c_pb_birth_table_create_budget(lo,hi,&i.options.birth,&p,&control,n,i.edges.data(),&raw)==0,"tiny budget table create");
 using Owner=std::unique_ptr<fusion_birth_table_v1,decltype(&fusion_c_birth_table_destroy)>;Owner table(raw,fusion_c_birth_table_destroy);
 std::array<const fusion_birth_table_v1*,5> tables{{table.get(),nullptr,nullptr,nullptr,nullptr}};
 fusion_birth_table_budget_info_v1 info{};need(fusion_c_pb_birth_table_budget_info(table.get(),&info)==0,"budget info");int matches=0;need(fusion_c_pb_birth_table_budget_matches_request(table.get(),lo,hi,&i.options.birth,&p,&control,n,i.edges.data(),&matches)==0&&matches==1,"exact table request");
 i.temperature=std::sqrt(lo*hi);Output out;need(call(i,out,&p,tables.data())==0,"budget table evaluation");need(out.accuracy.version==1&&out.accuracy.mode==2&&zero(out.accuracy.direct)&&same(out.accuracy.table,info),"table observation equals immutable sampled info");need(out.accuracy.number_L1_estimate==info.base.max_validated_number_L1&&out.accuracy.energy_L1_estimate==info.base.max_validated_energy_L1,"composite sampled estimates");
 rejected(i,&p,tables.data(),true);auto mismatch=p;mismatch.energy_L1_tolerance=std::nextafter(p.energy_L1_tolerance,1.);rejected(i,&mismatch,tables.data());
 auto source=i;source.options.birth.cutoff_J*=1.001;rejected(source,&p,tables.data());auto outside=i;outside.temperature=hi*1.001;rejected(outside,&p,tables.data());
 raw=nullptr;need(fusion_c_birth_table_create(0,lo,hi,&i.options.birth,&control,n,i.edges.data(),&raw)==0,"tiny full table create");Owner full(raw,fusion_c_birth_table_destroy);tables[0]=full.get();rejected(i,&p,tables.data());
}
}
#define main original_fast_tests_main
#include "test_fusion_coupled_fast.cpp"
#undef main
int main(){try{
 direct_contract();table_contract();
 auto in=make_subnormal_spill_inputs(2.9601203997565862e-291);
 auto old=make_trial(in),budget=make_trial(in),packet=make_trial(in);
 fusion_handoff_diagnostics_v1 d{};fusion_beam_table_usage_v1 u{};
 uint64_t outside=0;fusion_thermal_birth_budget_policy_v1 policy{1,4,4,0,0,0};
 fusion_coupled_source_accuracy_v1 accuracy{};
 auto evaluate=[&](Trial& o,int mode){
 const int cells=in.grid.cells();std::vector<double> mapped(70*cells);fusion_birth_packets_v1 meta{};
#define COVERED_ARGS 1e-4,&in.options,&in.fast_options,nullptr,0,nullptr,0,cells,in.grid.edges.data(),in.thermal_number.data(),in.electron_energy_J_m3,in.ion_energy_J_m3,in.electron_density_m3,in.thermal_charge_squared.data(),0,nullptr,in.coulomb_logs.data(),in.old_s.data(),in.old_t.data(),in.external_birth.data(),in.escape.data(),o.thermal_number.data(),o.s.data(),o.t.data(),&o.result,&d,&u,nullptr,nullptr,FUSION_BEAM_TABLE_STRICT,&outside
 if(mode==0)return fusion_c_coupled_sources_covered_trial(COVERED_ARGS);
 if(mode==1)return fusion_c_coupled_sources_covered_budget_trial(COVERED_ARGS,&policy,&accuracy);
 return fusion_c_coupled_sources_packets_trial(COVERED_ARGS,mapped.data(),&meta);
#undef COVERED_ARGS
 };
 require(evaluate(old,0)==0,"weak source old covered succeeds");
 require(evaluate(packet,2)!=0,"weak positive unrepresentable packets reject");
 require(evaluate(budget,1)==0,"weak source covered budget succeeds");
 require(old.s==budget.s&&old.t==budget.t&&old.thermal_number==budget.thermal_number&&same_result(old.result,budget.result),"weak covered physical parity");
 require(old.result.ledger.events_m3[3]>0&&accuracy.mode==0,"real weak DT source, inactive pB");
 std::cout<<"PASS covered budget finite-source contracts\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
