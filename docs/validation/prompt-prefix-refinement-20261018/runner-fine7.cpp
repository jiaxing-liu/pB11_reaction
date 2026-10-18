#include "fusion_prompt_orbit.h"
#include "fusion_boundary_segment.h"
#include "fusion_flux_inverse.h"
#include "fusion_flux_boundary.h"
#include "fusion_magnetic_push.h"
#include "fusion_nuclear_data.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <stdexcept>
#ifdef FUSION_USE_PREPARED
struct PreparedCache {
 const int sizes[5]={512,2048,8192,32768,65536};
 fusion_flux_boundary_prepared* contexts[5]{};
 ~PreparedCache(){for(auto*p:contexts)fusion_c_flux_boundary_destroy(p);}
 const fusion_flux_boundary_prepared* at(int n)const{for(int i=0;i<5;++i)if(sizes[i]==n)return contexts[i];return nullptr;}
} prepared;
#endif
constexpr double light=299792458.,e=1.602176634e-19;
int kcore,nodes,harmonics; std::vector<double> xi_data,coef_data,slope_data,F_data,C_data;
const double *xi,*coef,*slope,*F,*C;
long long field_calls=0,outside_field_requests=0;
#ifdef FUSION_TRACE_SEGMENTS
std::ofstream segment_trace;
int traced_marker=-1;
long long traced_call=0;
#endif
double length(const double*x){return std::hypot(std::hypot(x[0],x[1]),x[2]);}
double separation(const double*a,const double*b){double d[3];for(int j=0;j<3;++j)d[j]=a[j]-b[j];return length(d);}
int query_boundary(int N,double R,double Z,fusion_flux_boundary_value_v1*out){
#ifdef FUSION_USE_PREPARED
 return fusion_c_flux_boundary_prepared_point(prepared.at(N),R,Z,out);
#else
 return fusion_c_flux_boundary(harmonics,coef+(nodes-1)*(1+2*harmonics),N,R,Z,1e-12,1e-12,out);
#endif
}
int query_segment(int N,const double*a,const double*b,double tol,fusion_boundary_segment_value_v1*out){
#ifdef FUSION_USE_PREPARED
 return fusion_c_boundary_prepared_segment(prepared.at(N),a,b,tol,4096,out);
#else
 return fusion_c_boundary_segment(harmonics,coef+(nodes-1)*(1+2*harmonics),N,a,b,1e-12,1e-12,tol,4096,out);
#endif
}
int boundary(const double*x,double&clear){
 for(int N:{512,2048,8192,32768,65536}){
  fusion_flux_boundary_value_v1 b{};
  int rc=query_boundary(N,std::hypot(x[0],x[1]),x[2],&b);
  if(rc)return rc;
  clear=b.polygon_distance_m-b.uncertainty_m;
  if(b.classification==FUSION_BOUNDARY_INSIDE)return 0;
  if(b.classification==FUSION_BOUNDARY_OUTSIDE)return 2001;
 }
 return 2001;
}
int field(void*,const double*x,double*B){
 double clear;int rc=boundary(x,clear);if(rc){++outside_field_requests;return rc;}
 double R=std::hypot(x[0],x[1]);fusion_flux_inverse_value_v1 inv{};
 rc=fusion_c_flux_inverse(nodes,harmonics,xi,coef,slope,F,C,kcore,R,x[2],2e-14,60,0,0,0,&inv);if(rc)return rc;
 fusion_flux_field_value_v1 b{};rc=fusion_c_flux_core_eval(nodes,harmonics,xi,coef,slope,F,C,kcore,inv.xi,inv.theta_rad,&b);if(rc)return rc;
 B[0]=b.BR_T*x[0]/R-b.Bphi_T*x[1]/R;
 B[1]=b.BR_T*x[1]/R+b.Bphi_T*x[0]/R;B[2]=b.BZ_T;++field_calls;return 0;
}


int point_callback(void*,const double*x,int*kind){
 for(int N:{512,2048,8192,32768,65536}){
  fusion_flux_boundary_value_v1 b{};int rc=query_boundary(N,std::hypot(x[0],x[1]),x[2],&b);if(rc)return rc;
  *kind=b.classification;if(*kind!=FUSION_BOUNDARY_AMBIGUOUS)return 0;
 }return 0;
}
int segment_callback(void*,const double*a,const double*b,double tol,fusion_boundary_segment_value_v1*out){
#ifdef FUSION_TRACE_SEGMENTS
 ++traced_call;
 auto trace=[&](int N){
  segment_trace<<traced_marker<<','<<traced_call<<','<<N<<','<<tol<<','<<out->result<<','<<out->termination_reason<<','<<out->tolerance_met<<','<<out->intervals_visited<<','<<out->fraction_lo<<','<<out->fraction_hi<<','<<out->boundary_band_m;
  for(int j=0;j<3;++j)segment_trace<<','<<a[j];
  for(int j=0;j<3;++j)segment_trace<<','<<b[j];
  segment_trace<<'\n';segment_trace.flush();
 };
#endif
 for(int N:{512,2048,8192,32768,65536}){
  int rc=query_segment(N,a,b,tol,out);if(rc)return rc;
  if(out->result==FUSION_SEGMENT_CLEAR)return 0;
  if(out->result==FUSION_SEGMENT_GEOMETRIC_CANDIDATE&&out->fraction_hi-out->fraction_lo<=4*tol){
#ifdef FUSION_TRACE_SEGMENTS
   trace(N);
#endif
   return 0;
  }
 }
#ifdef FUSION_TRACE_SEGMENTS
 trace(65536);
#endif
 return 0;
}
int main(int argc,char**argv){try{
 if(argc!=2&&argc!=4&&argc!=5&&argc!=6)throw std::runtime_error("usage: runner case.dat [horizon_s initial_dt_s [max_levels [max_steps_per_level]]]");
 // Short, static E=0 chain-integration gate. These are not physical loss settings.
 fusion_prompt_options_v1 o{1e-8,1e-9,1e-14,1e-12,1e-5,1e-5,1e-4,7,7,1000000};
 if(argc>=4){
  size_t used=0;o.horizon_s=std::stod(argv[2],&used);if(used!=std::string(argv[2]).size())throw std::runtime_error("horizon text");
  o.initial_dt_s=std::stod(argv[3],&used);if(used!=std::string(argv[3]).size())throw std::runtime_error("dt text");
 }
 if(argc>=5){
  size_t used=0;int levels=std::stoi(argv[4],&used);
  if(used!=std::string(argv[4]).size()||levels<o.min_levels||levels>20)throw std::runtime_error("max levels");
  o.max_levels=levels;
 }
 if(argc==6){
  size_t used=0;int steps=std::stoi(argv[5],&used);
  if(used!=std::string(argv[5]).size()||steps<1||steps>1000000)throw std::runtime_error("max steps per level");
  o.max_steps_per_level=steps;
 }
 if(!std::isfinite(o.horizon_s)||o.horizon_s<=0||!std::isfinite(o.initial_dt_s)||o.initial_dt_s<=0)throw std::runtime_error("time controls");
 std::ifstream input(argv[1]);input>>nodes>>harmonics>>kcore;
 if(nodes<3||nodes>300||harmonics<1||harmonics>7||kcore<1||kcore>=nodes-1)throw std::runtime_error("shape");
 auto read=[&](std::vector<double>&v,int n){v.resize(n);for(double&x:v){input>>x;if(!input||!std::isfinite(x))throw std::runtime_error("input");}};
 read(xi_data,nodes);read(coef_data,nodes*(1+2*harmonics));read(slope_data,nodes*(1+2*harmonics));read(F_data,nodes);read(C_data,nodes);
 xi=xi_data.data();coef=coef_data.data();slope=slope_data.data();F=F_data.data();C=C_data.data();
#ifdef FUSION_USE_PREPARED
 for(int i=0;i<5;++i)if(fusion_c_flux_boundary_prepare(harmonics,coef+(nodes-1)*(1+2*harmonics),prepared.sizes[i],1e-12,1e-12,&prepared.contexts[i]))throw std::runtime_error("prepare boundary");
#endif
 int count;input>>count;if(count<1||count>1000000)throw std::runtime_error("count");
#ifdef FUSION_TRACE_SEGMENTS
 segment_trace.open("segment-trace.csv");
 if(!segment_trace)throw std::runtime_error("segment trace output");
 segment_trace<<std::setprecision(17)<<"marker,callback,segments,fraction_tolerance,result,termination_reason,tolerance_met,intervals_visited,fraction_lo,fraction_hi,boundary_band_m,ax,ay,az,bx,by,bz\n";
#endif
 std::ofstream csv("driver.csv");csv<<std::setprecision(17)<<"marker,number,kinetic_J,api_status,outcome,reason,levels,steps,time_lo_s,time_hi_s,time_spread_s,position_spread_m,field_calls,outside_field_requests\n";
 std::ofstream settings("options.json");settings<<std::setprecision(17)<<"{\"horizon_s\":"<<o.horizon_s<<",\"initial_dt_s\":"<<o.initial_dt_s
 <<",\"geometry_time_budget_s\":"<<o.geometry_time_budget_s<<",\"event_time_tolerance_s\":"<<o.event_time_tolerance_s
 <<",\"event_position_tolerance_m\":"<<o.event_position_tolerance_m<<",\"final_position_tolerance_m\":"<<o.final_position_tolerance_m
 <<",\"final_u_relative_tolerance\":"<<o.final_u_relative_tolerance<<",\"min_levels\":"<<o.min_levels<<",\"max_levels\":"<<o.max_levels
 <<",\"max_steps_per_level\":"<<o.max_steps_per_level<<"}\n";settings.close();
 std::ofstream states("final-state.csv");states<<std::setprecision(17)<<"marker,outcome,state_valid,final_time_s,x,y,z,ux,uy,uz,endpoint_clearance_m,u_relative_drift,estimated_u_relative_spread\n";
 for(int j=0;j<count;++j){
#ifdef FUSION_TRACE_SEGMENTS
  traced_marker=j;
#endif
  int id;double mass,charge,N,K;fusion_orbit_state_v1 initial{};
  input>>id>>mass>>charge>>N>>K;for(double&x:initial.x_m)input>>x;for(double&u:initial.u_m_s)input>>u;
  if(!input||id!=j)throw std::runtime_error("marker input");
  auto before=field_calls;auto outside=outside_field_requests;
  fusion_prompt_result_v1 r{};int rc=fusion_c_prompt_orbit(mass,charge,&initial,&o,field,point_callback,segment_callback,nullptr,&r);
  csv<<id<<','<<N<<','<<K<<','<<rc<<','<<r.outcome<<','<<r.reason<<','<<r.levels_used<<','<<r.steps_last_level<<','<<r.event_time_lo_s<<','<<r.event_time_hi_s<<','<<r.estimated_time_spread_s<<','<<r.estimated_position_spread_m<<','<<field_calls-before<<','<<outside_field_requests-outside<<'\n';csv.flush();
  if(rc)throw std::runtime_error("orbit API error; do not count as unresolved");
  bool valid=r.outcome==FUSION_PROMPT_RETAINED;
  states<<id<<','<<r.outcome<<','<<valid;
  if(valid){
   double clearance=0;if(boundary(r.final_state.x_m,clearance))throw std::runtime_error("retained final boundary");
   states<<','<<r.final_time_s;for(double v:r.final_state.x_m)states<<','<<v;for(double v:r.final_state.u_m_s)states<<','<<v;
   double initial_u=length(initial.u_m_s),final_u=length(r.final_state.u_m_s);
   if(initial_u==0&&final_u!=0)throw std::runtime_error("stationary drift");
   double drift=initial_u>0?std::abs(final_u/initial_u-1):0.;
   states<<','<<clearance<<','<<drift<<','<<r.estimated_u_relative_spread;
  }else states<<",,,,,,,,,,"; // Empty unsupported fields, not physical zero values.
  states<<'\n';states.flush();
  if(j%50==0)std::cout<<j<<'/'<<count<<" outcome="<<r.outcome<<std::endl;
 }
 std::string extra;if(input>>extra)throw std::runtime_error("trailing input");
 if(outside_field_requests)throw std::runtime_error("outside field request");
 std::cout<<"PASS selected subset only; no total loss fraction"<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
