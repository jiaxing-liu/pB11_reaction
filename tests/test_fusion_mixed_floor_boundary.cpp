#include "fusion_coupled_sources.h"
#include "fusion_thermal_burn.h"
#include "fusion_nuclear_data.h"
#include "fusion_rate_model.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
static int provider_calls=0;
extern "C" int __wrap_fusion_c_thermal_birth_grid(int ch,double,const fusion_thermal_birth_options_v1*,int n,const double*,double*grid,fusion_thermal_birth_v1*out){
 ++provider_calls; if(ch!=3) return PB11_STATUS_INVALID_ARGUMENT;
 std::fill(grid,grid+7*n,0.); *out={}; out->reactivity_m3_s=1e-10;
 const double q=std::numeric_limits<double>::denorm_min();
 out->below_number_m3_s[0]=1990*q; out->below_number_m3_s[4]=2*q;
 return 0;
}
void check(bool x,const char*m){if(!x)throw std::runtime_error(m);}
template<class T>void poison(T&v){std::memset(&v,0xa5,sizeof(v));}
template<class T>void zero(T&v,const char*m){const unsigned char*p=reinterpret_cast<const unsigned char*>(&v); for(size_t i=0;i<sizeof(v);++i)check(p[i]==0,m);}
int main(int argc,char**argv){try{
 check(argc==2,"mode required"); const bool high=std::string(argv[1])=="high";
 const double tiny=std::numeric_limits<double>::denorm_min(),dt=.1,Ui=.1991,center=.001;
 double edges[2]={0,2*center},thermal[6]={0,1,1,0,0,0},charge[6]={1,1,1,4,4,25};
 double old[6]={},logs[42];std::fill(logs,logs+42,15.);
 fusion_coupled_thermal_options_v1 op{};
 op.birth.relative_max_J=5*1.602176634e-13;op.birth.cm_max_kT=40;op.birth.ground_state_q_J=91.84*1.602176634e-16;op.birth.cutoff_J=.002*1.602176634e-13;
 op.birth.l1_fraction=.76;op.birth.narrow_peak_fraction=.051;op.birth.continuum_peak_scale=1;
 op.birth.continuation=FUSION_ENDPOINT_S;op.birth.pb_low=FUSION_PB_LOW_TB;op.birth.broad_mode=13;
 op.birth.relative_order=8;op.birth.cm_order=8;op.birth.nq=8;op.birth.ncos=8;
 op.max_source_rate_error=.002;op.max_source_debit_error=.002;op.handoff_max_L1=.001;op.handoff_max_mean_error=.001;op.channels[3]=1;
 fusion_fast_target_options_v1 fast{};fast.angular_order=8;fast.angular_max_exponent=40;
 // Real burn independently determines the actual source scale used by coupled.
 double energy[6]={0,Ui/2,Ui/2,0,0,0},rates[5]={0,0,0,1e-10,0},means[5]={},newN[6],newU[6]; fusion_thermal_burn_v1 burn{};
 check(fusion_c_thermal_burn_trial(dt,thermal,energy,rates,means,means,newN,newU,&burn)==0,"real burn");
 const long double scale=(long double)burn.events_m3[3]/rates[3];
 const long double normalN=scale*(1990*tiny),preciseN=scale*(2*tiny);
 const long double normalBorrow=(long double)(double)normalN*center,preciseBorrow=preciseN*center,limit=(long double)tiny*Ui;
 check((double)normalN>0 && (double)normalBorrow==0,"normal mapped energy underflows");
 check(preciseN>0 && (double)preciseN==0 && (double)(preciseN/dt)>0,"precise rate-only number");
 check(normalBorrow<limit && preciseBorrow<limit && normalBorrow+preciseBorrow>limit,"mixed strict boundary separates");
 std::cout<<std::setprecision(24)<<"fixture scale="<<scale<<" normalN="<<normalN<<" preciseN="<<preciseN<<" normalBorrow="<<normalBorrow<<" preciseBorrow="<<preciseBorrow<<" limit="<<limit<<" total="<<normalBorrow+preciseBorrow<<"\n";
 double outputN[6],s[6],t[6],packets[70];fusion_coupled_thermal_v1 result;fusion_handoff_diagnostics_v1 diagnostics;fusion_beam_table_usage_v1 usage;fusion_birth_floor_ledger_v1 ledger;fusion_birth_packets_v1 meta;uint64_t outside;
 poison(outputN);poison(s);poison(t);poison(packets);poison(result);poison(diagnostics);poison(usage);poison(ledger);poison(meta);poison(outside);
 fusion_coupled_floor_limits_v1 limits{.1,high?2*tiny:tiny};
 int status=fusion_c_coupled_sources_packets_trial(dt,&op,&fast,nullptr,0,nullptr,0,1,edges,thermal,.2,Ui,2,charge,0,nullptr,logs,old,old,old,old,outputN,s,t,&result,&diagnostics,&usage,&limits,&ledger,FUSION_BEAM_TABLE_STRICT,&outside,packets,&meta);
 check(provider_calls==1,"actual direct thermal provider call");
 check(status==(high?PB11_STATUS_NUMERICAL_FAILURE:PB11_STATUS_OUT_OF_RANGE),"synthetic mixed borrowing boundary status");
 zero(outputN,"clear thermal");zero(s,"clear s");zero(t,"clear t");zero(packets,"clear packet array");zero(result,"clear result");zero(diagnostics,"clear diagnostics");zero(usage,"clear usage");zero(ledger,"clear floor");zero(meta,"clear metadata");zero(outside,"clear outside counter");
 std::cout<<"RESULT mode="<<argv[1]<<" status="<<status<<" all_outputs_zero=1 physics_accepted=0\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
