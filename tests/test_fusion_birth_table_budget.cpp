// API/serialization engineering fixture, not a source-convergence certificate.
// Target-local link requirement: -Wl,--wrap=fusion_c_pb_thermal_birth_budget_grid
#include "fusion_birth_table_budget.h"
#include "fusion_birth_table_internal.h"
#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include <utility>
namespace {unsigned direct_calls=0;}
extern "C" int __real_fusion_c_pb_thermal_birth_budget_grid(double,const fusion_thermal_birth_options_v1*,const fusion_thermal_birth_budget_policy_v1*,int,const double*,double*,fusion_thermal_birth_v1*,fusion_thermal_birth_budget_report_v1*);
extern "C" int __wrap_fusion_c_pb_thermal_birth_budget_grid(double t,const fusion_thermal_birth_options_v1*o,const fusion_thermal_birth_budget_policy_v1*p,int n,const double*e,double*b,fusion_thermal_birth_v1*r,fusion_thermal_birth_budget_report_v1*q){++direct_calls;return __real_fusion_c_pb_thermal_birth_budget_grid(t,o,p,n,e,b,r,q);}
namespace {
constexpr double kev=1.602176634e-16,mev=1.602176634e-13;
constexpr int cells=8;constexpr double lower=.5*kev,upper=lower*1.000001;
using Grid=std::array<double,7*cells>;using Bytes=std::vector<unsigned char>;
using Owner=std::unique_ptr<fusion_birth_table_v1,decltype(&fusion_c_birth_table_destroy)>;
const std::array<double,cells+1> edges{{0,.1*mev,.2*mev,.4*mev,.8*mev,1.2*mev,1.6*mev,2*mev,2.4*mev}};
void need(bool b,const char*m){if(!b)throw std::runtime_error(m);}
void ok(int s,const char*m){need(s==PB11_STATUS_OK,m);}
template<class T>bool cleared(const T&v){const auto*p=reinterpret_cast<const unsigned char*>(&v);return std::all_of(p,p+sizeof(v),[](unsigned char c){return c==0;});}
fusion_thermal_birth_options_v1 source(){fusion_thermal_birth_options_v1 o{};o.relative_max_J=2*mev;o.cm_max_kT=40;o.ground_state_q_J=91.84*kev;o.cutoff_J=.001*mev;o.l1_fraction=.76;o.narrow_peak_fraction=.051;o.continuum_peak_scale=1;o.continuation=FUSION_ENDPOINT_S;o.pb_low=FUSION_PB_LOW_TB;o.broad_mode=13;o.relative_order=4;o.cm_order=4;o.nq=o.ncos=8;return o;}
fusion_birth_table_control_v1 control(){
 // Declared before execution: this tiny interval probes API/interpolation,
 // not the accuracy of relative/CM4 integration of the natural rate model.
 return {.01,.01,.01,.01,.05,.05,8,64,3};
}
fusion_thermal_birth_budget_policy_v1 policy(){return {1,4,4,0,1e-4,1e-4};}
Bytes pack(const fusion_birth_table_v1*t){size_t n=0,w=99;ok(fusion_c_birth_table_pack_size(t,&n),"pack size");Bytes b(n);ok(fusion_c_birth_table_pack(t,b.data(),b.size(),&w),"pack");need(w==n,"written bytes");return b;}
uint64_t word(const Bytes&b,size_t off){need(off+8<=b.size(),"offset");uint64_t x=0;for(int i=0;i<8;++i)x|=uint64_t(b[off+i])<<(8*i);return x;}
void put(Bytes&b,size_t off,uint64_t x){need(off+8<=b.size(),"offset write");for(int i=0;i<8;++i)b[off+i]=static_cast<unsigned char>(x>>(8*i));}
void realword(Bytes&b,size_t off,double x){uint64_t bits;std::memcpy(&bits,&x,8);put(b,off,bits);}
void checksum(Bytes&b){uint64_t h=14695981039346656037ULL;for(size_t i=0;i+8<b.size();++i){h^=b[i];h*=1099511628211ULL;}put(b,b.size()-8,h);}
void reject(Bytes b,bool repair=true){if(repair)checksum(b);fusion_birth_table_v1*t=reinterpret_cast<fusion_birth_table_v1*>(uintptr_t(1));unsigned before=direct_calls;need(fusion_c_birth_table_unpack(b.data(),b.size(),&t)!=0,"malformed table imported");need(t==nullptr,"failed import atomic null");need(direct_calls==before,"import repeated source quadrature");}
Owner create(const fusion_thermal_birth_budget_policy_v1&p){auto o=source();auto c=control();fusion_birth_table_v1*t=nullptr;ok(fusion_c_pb_birth_table_create_budget(lower,upper,&o,&p,&c,cells,edges.data(),&t),"budget create");need(t!=nullptr,"published handle");return Owner(t,fusion_c_birth_table_destroy);}
void match(const fusion_birth_table_v1*t){auto o=source();auto c=control();auto p=policy();int matches=-1;
 auto check=[&](double lo,double hi,const fusion_thermal_birth_options_v1&s,const fusion_thermal_birth_budget_policy_v1&q,const fusion_birth_table_control_v1&ctrl,const double*e,bool expected){matches=99;ok(fusion_c_pb_birth_table_budget_matches_request(t,lo,hi,&s,&q,&ctrl,cells,e,&matches),"explicit request compare");need(matches==int(expected),"exact request matching");};
 check(lower,upper,o,p,c,edges.data(),true);
 for(int field=0;field<6;++field){auto q=p;switch(field){case 0:q.version=2;break;case 1:q.coarse_nq=5;break;case 2:q.coarse_ncos=5;break;case 3:q.reserved=1;break;case 4:q.number_L1_tolerance*=2;break;case 5:q.energy_L1_tolerance*=2;break;}matches=99;int rc=fusion_c_pb_birth_table_budget_matches_request(t,lower,upper,&o,&q,&c,cells,edges.data(),&matches);need((rc==0&&matches==0)||(rc!=0&&matches==0),"every policy field distinguished");}
 auto changed=o;changed.narrow_peak_fraction+=.001;check(lower,upper,changed,p,c,edges.data(),false);auto cc=c;cc.max_energy_L1*=.5;check(lower,upper,o,p,cc,edges.data(),false);auto ee=edges;ee[1]*=.9;check(lower,upper,o,p,c,ee.data(),false);check(lower*.99,upper,o,p,c,edges.data(),false);check(lower,upper*1.01,o,p,c,edges.data(),false);
 matches=99;ok(fusion_c_birth_table_matches_request(t,0,lower,upper,&o,&c,cells,edges.data(),&matches),"old matcher");need(matches==0,"old public matcher rejects budget");need(!fusion_detail::birth_table_matches(t,0,o,cells,edges.data()),"old internal matcher rejects budget");
}
void check_reports_and_interpolation(const fusion_birth_table_v1*t){auto o=source();auto p=policy();auto c=control();fusion_birth_table_budget_info_v1 info{};ok(fusion_c_pb_birth_table_budget_info(t,&info),"budget info");fusion_birth_table_info_v1 generic{};ok(fusion_c_birth_table_info(t,&generic),"generic info");need(generic.knots==info.base.knots&&generic.max_validated_number_L1==info.base.max_validated_number_L1&&generic.max_validated_energy_L1==info.base.max_validated_energy_L1,"generic combined info");
 need(info.base.max_validated_number_L1<=c.max_number_L1&&info.base.max_validated_energy_L1<=c.max_energy_L1,"combined gates");need(info.max_sampled_source_number_bound<=c.max_number_L1&&info.max_sampled_source_energy_bound<=c.max_energy_L1,"direct source gates");need(info.base.max_validated_number_L1>=info.max_validated_interpolation_number_L1&&info.base.max_validated_energy_L1>=info.max_validated_interpolation_energy_L1,"separate interpolation accounting");
 for(int k=0;k<info.base.knots;++k){double T=0;fusion_thermal_birth_budget_report_v1 report{};ok(fusion_c_pb_birth_table_budget_knot_report(t,k,&T,&report),"knot report");need(T>=lower&&T<=upper&&report.version==1&&report.high_pairs+report.low_pairs==report.joint_nodes,"knot metadata");need(report.full_nq==8&&report.full_ncos==8&&report.coarse_nq==4&&report.coarse_ncos==4,"knot orders");need(report.number_L1_bound+report.number_roundoff_allowance<=info.max_sampled_source_number_bound+1e-18&&report.energy_L1_bound+report.energy_roundoff_allowance<=info.max_sampled_source_energy_bound+1e-18,"effective source accounting includes allowance");Grid direct{},table{};fusion_thermal_birth_v1 out{};fusion_thermal_birth_budget_report_v1 r{};fusion_birth_coefficients_v1 coeff{};ok(fusion_c_pb_thermal_birth_budget_grid(T,&o,&p,cells,edges.data(),direct.data(),&out,&r),"direct stored knot");ok(fusion_c_birth_table_evaluate(t,T,cells,table.data(),&coeff),"stored knot evaluate");need(std::memcmp(direct.data(),table.data(),sizeof(Grid))==0,"stored knot exact source");need(r.number_L1_bound==report.number_L1_bound&&r.energy_L1_bound==report.energy_L1_bound&&r.high_pairs==report.high_pairs,"stored report reproduction");}
 const double T=std::sqrt(lower*upper);Grid direct{},interpolated{};fusion_thermal_birth_v1 d{};fusion_thermal_birth_budget_report_v1 r{};fusion_birth_coefficients_v1 q{};ok(fusion_c_pb_thermal_birth_budget_grid(T,&o,&p,cells,edges.data(),direct.data(),&d,&r),"midpoint direct");ok(fusion_c_birth_table_evaluate(t,T,cells,interpolated.data(),&q),"midpoint interpolate");long double dn=0,du=0;for(int s=0;s<7;++s){for(int i=0;i<cells;++i){auto n=s*cells+i;need(std::isfinite(interpolated[n])&&interpolated[n]>=0,"positive interpolation");long double delta=std::abs(static_cast<long double>(direct[n])-interpolated[n]);dn+=delta;du+=delta*(static_cast<long double>(edges[i])+edges[i+1])/2;}for(int side=0;side<2;++side){dn+=std::abs(static_cast<long double>(side?d.above_number_m3_s[s]:d.below_number_m3_s[s])-(side?q.above_number_m3_s[s]:q.below_number_m3_s[s]));du+=std::abs(static_cast<long double>(side?d.above_energy_J_m3_s[s]:d.below_energy_J_m3_s[s])-(side?q.above_energy_J_m3_s[s]:q.below_energy_J_m3_s[s]));}}
 fusion_nuclear_channel_v1 reaction{};ok(fusion_c_nuclear_channel(0,&reaction),"Q");long double N=3*static_cast<long double>(d.reactivity_m3_s),U=static_cast<long double>(reaction.q_J)*d.reactivity_m3_s+d.reactant_energy_moment_J_m3_s[0]+d.reactant_energy_moment_J_m3_s[1];need(N>0&&U>0,"expected ledger nonzero");need(dn/N+r.number_L1_bound+1e-12L<=c.max_number_L1&&du/U+r.energy_L1_bound+1e-12L<=c.max_energy_L1,"raw interpolation plus source budget plus engineering allowance");
}
void bytes(const fusion_birth_table_v1*t){auto b=pack(t);fusion_birth_table_budget_info_v1 info{};ok(fusion_c_pb_birth_table_budget_info(t,&info),"info");need(word(b,8)==2&&word(b,16)==b.size(),"budget format2 header");const size_t header=88,basewords=38,policyOffset=header+basewords*8,edgeOffset=policyOffset+12*8,nodeOffset=edgeOffset+(cells+1)*8,nodeBytes=(51+7*cells)*8;
 need(b.size()==nodeOffset+info.base.knots*nodeBytes+8,"format2 exact schema size");
 fusion_birth_table_v1*raw=nullptr;unsigned before=direct_calls;ok(fusion_c_birth_table_unpack(b.data(),b.size(),&raw),"generic budget unpack");Owner copy(raw,fusion_c_birth_table_destroy);need(direct_calls==before,"unpack performs no new source integration");need(pack(copy.get())==b,"all serialized bytes roundtrip");match(copy.get());
 std::array<unsigned char,32> da{},db{};size_t na=0,nb=0;ok(fusion_c_birth_table_content_digest(t,da.data(),&na),"original digest");ok(fusion_c_birth_table_content_digest(copy.get(),db.data(),&nb),"import digest");need(da==db&&na==nb&&na==b.size(),"digest/length stable");
 for(double T:{lower,std::sqrt(lower*upper),upper}){Grid a{},c{};fusion_birth_coefficients_v1 oa{},oc{};ok(fusion_c_birth_table_evaluate(t,T,cells,a.data(),&oa),"original evaluate");ok(fusion_c_birth_table_evaluate(copy.get(),T,cells,c.data(),&oc),"copy evaluate");need(std::memcmp(a.data(),c.data(),sizeof(Grid))==0&&std::memcmp(&oa,&oc,sizeof(oa))==0,"generic coefficients exact after import");}
 for(int k=0;k<info.base.knots;++k){double x=0,y=0;fusion_thermal_birth_budget_report_v1 a{},c{};ok(fusion_c_pb_birth_table_budget_knot_report(t,k,&x,&a),"original knot");ok(fusion_c_pb_birth_table_budget_knot_report(copy.get(),k,&y,&c),"copy knot");need(x==y&&a.high_pairs==c.high_pairs&&a.number_L1_bound==c.number_L1_bound&&a.energy_L1_bound==c.energy_L1_bound,"import reports exact");}
 Bytes small(b.size()-1,0x5a);size_t written=99;need(fusion_c_birth_table_pack(t,small.data(),small.size(),&written)!=0&&written==0&&std::all_of(small.begin(),small.end(),[](unsigned char x){return x==0x5a;}),"short pack untouched");
 auto bad=b;bad.back()^=1;reject(bad,false);bad=b;bad.pop_back();reject(bad,false);bad=b;bad.push_back(0);reject(bad,false);bad=b;put(bad,8,3);reject(bad);bad=b;bad[24]^=1;reject(bad);bad=b;put(bad,16,b.size()-8);reject(bad);
 for(auto change: {std::pair<size_t,uint64_t>{policyOffset,2},{policyOffset+24,1},{nodeOffset+(1+31)*8,2},{nodeOffset+(1+31+6+2)*8,UINT64_MAX}}){bad=b;put(bad,change.first,change.second);reject(bad);}
 bad=b;realword(bad,policyOffset+32,3);reject(bad);bad=b;realword(bad,nodeOffset+(1+31+6+7)*8,-1);reject(bad);
 // Coordinated, checksum-valid forgery: a zero policy may not retain low
 // parents even if both stored bounds are changed to zero consistently.
 bad=b;realword(bad,policyOffset+32,0);realword(bad,policyOffset+40,0);
 realword(bad,policyOffset+48,1e-12);realword(bad,policyOffset+56,1e-12);
 bool retained_low=false;
 for(int k=0;k<info.base.knots;++k){const size_t node=nodeOffset+size_t(k)*nodeBytes;
  retained_low|=word(bad,node+(1+31+6+4)*8)>0;
  realword(bad,node+(1+31+6+7)*8,0);realword(bad,node+(1+31+6+7+1)*8,0);
 }
 need(retained_low,"zero-policy regression fixture must retain nonzero low parents");
 reject(bad); // reject() repairs FNV; v1 deliberately fails this assertion.

}
void errors(const fusion_birth_table_v1*t){auto o=source();auto p=policy();auto c=control();fusion_birth_table_v1*out=reinterpret_cast<fusion_birth_table_v1*>(uintptr_t(1));p.version=2;need(fusion_c_pb_birth_table_create_budget(lower,upper,&o,&p,&c,cells,edges.data(),&out)!=0&&out==nullptr,"failed create nulls handle");need(fusion_c_pb_birth_table_create_budget(lower,upper,&o,&p,&c,cells,edges.data(),nullptr)==PB11_STATUS_NULL_OUTPUT,"null create output");fusion_birth_table_budget_info_v1 info;std::memset(&info,0x5a,sizeof(info));need(fusion_c_pb_birth_table_budget_info(nullptr,&info)!=0&&cleared(info),"null info clears");
 fusion_thermal_birth_budget_report_v1 r;double T=99;std::memset(&r,0x5a,sizeof(r));need(fusion_c_pb_birth_table_budget_knot_report(t,-1,&T,&r)!=0&&T==0&&cleared(r),"bad index clears both");T=99;need(fusion_c_pb_birth_table_budget_knot_report(t,0,&T,nullptr)==PB11_STATUS_NULL_OUTPUT&&T==0,"null knot report clears T");std::memset(&r,0x5a,sizeof(r));need(fusion_c_pb_birth_table_budget_knot_report(t,0,nullptr,&r)==PB11_STATUS_NULL_OUTPUT&&cleared(r),"null knot T clears report");int m=99;need(fusion_c_pb_birth_table_budget_matches_request(nullptr,lower,upper,&o,&p,&c,cells,edges.data(),&m)!=0&&m==0,"failed match clears");
}
void source_gate(){auto o=source();auto c=control();auto p=policy();p.number_L1_tolerance=p.energy_L1_tolerance=2;fusion_birth_table_v1*t=reinterpret_cast<fusion_birth_table_v1*>(uintptr_t(1));need(fusion_c_pb_birth_table_create_budget(lower,upper,&o,&p,&c,cells,edges.data(),&t)!=0&&t==nullptr,"endpoint source bound alone must fit controls");c.max_number_L1=1.01;p=policy();t=reinterpret_cast<fusion_birth_table_v1*>(uintptr_t(1));need(fusion_c_pb_birth_table_create_budget(lower,upper,&o,&p,&c,cells,edges.data(),&t)!=0&&t==nullptr,"controls retain range0..1");}
void full_and_zero(){auto o=source();auto c=control();fusion_birth_table_v1*raw=nullptr;ok(fusion_c_birth_table_create(0,lower,upper,&o,&c,cells,edges.data(),&raw),"legacy full constructor");Owner full(raw,fusion_c_birth_table_destroy);need(word(pack(full.get()),8)==1,"full remains format1");fusion_birth_table_budget_info_v1 inf;std::memset(&inf,0x5a,sizeof(inf));need(fusion_c_pb_birth_table_budget_info(full.get(),&inf)!=0&&cleared(inf),"full does not acquire implicit policy");auto p=policy();p.number_L1_tolerance=p.energy_L1_tolerance=0;auto zero=create(p);need(word(pack(zero.get()),8)==2,"zero budget still format2");int m=99;ok(fusion_c_birth_table_matches_request(zero.get(),0,lower,upper,&o,&c,cells,edges.data(),&m),"zero old matcher");need(m==0&&!fusion_detail::birth_table_matches(zero.get(),0,o,cells,edges.data()),"zero budget remains explicit");for(double T:{lower,upper}){Grid a{},b{};fusion_birth_coefficients_v1 x{},y{};ok(fusion_c_birth_table_evaluate(full.get(),T,cells,a.data(),&x),"full endpoint");ok(fusion_c_birth_table_evaluate(zero.get(),T,cells,b.data(),&y),"zero endpoint");need(std::memcmp(a.data(),b.data(),sizeof(Grid))==0&&std::memcmp(&x,&y,sizeof(x))==0,"zero full source identical");}}
}
int main(){try{auto table=create(policy());need(direct_calls>0,"test wrapper sees actual creation");match(table.get());check_reports_and_interpolation(table.get());bytes(table.get());errors(table.get());source_gate();full_and_zero();std::cout<<"PASS explicit budget birth-table API, interpolation/error accounting, versioned import and no-import-quadrature contracts; engineering fixture, not physical convergence\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
