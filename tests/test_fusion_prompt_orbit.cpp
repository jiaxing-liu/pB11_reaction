#include "fusion_prompt_orbit.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <algorithm>
void ck(bool ok,const char*msg){if(!ok)throw std::runtime_error(msg);}
struct Context {double radius=2.4,B=0;int fields=0,points=0,segments=0;int field_error=0,segment_error=0;bool ambiguous=false,false_clear=false,malformed=false;};
double norm(const double*x){return std::hypot(std::hypot(x[0],x[1]),x[2]);}
int point(void*p,const double*x,int*k){auto&c=*static_cast<Context*>(p);++c.points;double r=norm(x);*k=c.malformed?99:(c.ambiguous||r==c.radius?FUSION_BOUNDARY_AMBIGUOUS:r<c.radius?FUSION_BOUNDARY_INSIDE:FUSION_BOUNDARY_OUTSIDE);return 0;}
int field(void*p,const double*x,double*B){auto&c=*static_cast<Context*>(p);++c.fields;if(norm(x)>=c.radius)return 999;B[0]=B[1]=0;B[2]=c.B;return c.field_error;}
int segment(void*p,const double*a,const double*b,double tol,fusion_boundary_segment_value_v1*out){
 auto&c=*static_cast<Context*>(p);++c.segments;*out={};if(c.segment_error)return c.segment_error;
 if(c.ambiguous)return 0;
 if(c.false_clear||norm(b)<c.radius){out->result=FUSION_SEGMENT_CLEAR;return 0;}
 // Independent analytic convex sphere: bracket first quadratic root by bisection.
 if(norm(a)>=c.radius)return PB11_STATUS_OUT_OF_RANGE;
 double lo=0,hi=1;
 for(int i=0;i<60&&hi-lo>tol*.5;++i){double m=(lo+hi)/2,x[3];if(m==lo||m==hi)break;for(int j=0;j<3;++j)x[j]=(1-m)*a[j]+m*b[j];if(norm(x)<c.radius)lo=m;else hi=m;}
 double x[3];for(int j=0;j<3;++j)x[j]=(1-hi)*a[j]+hi*b[j];
 if(norm(x)==c.radius){hi=std::min(1.,hi+tol*.125);for(int j=0;j<3;++j)x[j]=(1-hi)*a[j]+hi*b[j];}
 if(norm(x)<=c.radius)return 0;
 out->result=FUSION_SEGMENT_GEOMETRIC_CANDIDATE;out->fraction_lo=lo;out->fraction_hi=hi;out->tolerance_met=hi-lo<=tol;return 0;
}
fusion_prompt_options_v1 options(){return {1.2,.1,1e-10,1e-6,1e-6,1e-7,1e-8,3,14,100000};}
int call(const fusion_orbit_state_v1&s,const fusion_prompt_options_v1&o,Context&c,fusion_prompt_result_v1&r){auto copy=s;int rc=fusion_c_prompt_orbit(1,1,&s,&o,field,point,segment,&c,&r);ck(std::memcmp(&s,&copy,sizeof(s))==0,"immutable input");return rc;}
int main(){try{
 fusion_orbit_state_v1 s{{2,0,0},{.5,0,0}};fusion_prompt_result_v1 r{};auto o=options();Context c;
 ck(call(s,o,c,r)==0&&r.outcome==FUSION_PROMPT_EVENT,"straight event");
 double exact=.4/(.5/std::hypot(1.,.5/299792458.));
 ck(r.event_time_lo_s<=exact+1e-14&&r.event_time_hi_s>=exact-1e-14,"straight analytic time");ck(r.levels_used>=3,"minimum levels");
 auto saved=r;Context repeat;ck(call(s,o,repeat,r)==0&&std::memcmp(&saved,&r,sizeof(r))==0,"deterministic replay");
 o.horizon_s=.2;c={};ck(call(s,o,c,r)==0&&r.outcome==FUSION_PROMPT_RETAINED,"straight retained");ck(r.final_time_s==.2&&std::abs(r.final_state.x_m[0]-2.1)<2e-14,"retained exact");
 // Same B=0 path would agree immediately; min_levels prevents early return.
 o.min_levels=5;c={};ck(call(s,o,c,r)==0&&r.levels_used==5,"explicit minimum level");
 o=options();c={};c.B=1;ck(call(s,o,c,r)==0&&r.outcome==FUSION_PROMPT_EVENT,"uniform B event");
 double a=0,b=1.2;for(int i=0;i<60;++i){double t=(a+b)/2,x=2+.5*std::sin(t),y=.5*(std::cos(t)-1);if(std::hypot(x,y)<2.4)a=t;else b=t;}
 double t=(a+b)/2,tc=(r.event_time_lo_s+r.event_time_hi_s)/2;
 ck(std::abs(tc-t)<o.event_time_tolerance_s,"independent uniform B time");
 double xc=(r.event_x_lo_m[0]+r.event_x_hi_m[0])/2,yc=(r.event_x_lo_m[1]+r.event_x_hi_m[1])/2;
 ck(std::hypot(xc-(2+.5*std::sin(t)),yc-.5*(std::cos(t)-1))<o.event_position_tolerance_m,"independent uniform B position");
 std::cout<<"uniform B levels="<<r.levels_used<<" time_error="<<std::abs(tc-t)<<" spread="<<r.estimated_time_spread_s<<'\n';
 s={{2.399,0,0},{1,0,0}};o=options();o.initial_dt_s=1;o.horizon_s=1;c={};ck(call(s,o,c,r)==0&&r.outcome==FUSION_PROMPT_EVENT&&r.phase==1&&c.fields==0,"phase1 no field query");
 s={{2,0,0},{.5,0,0}};o=options();c={};c.ambiguous=true;ck(call(s,o,c,r)==0&&r.outcome==0&&c.fields==0&&r.reason==FUSION_PROMPT_BOUNDARY_UNRESOLVED,"ambiguity not exit");
 c={};c.field_error=1234;r=saved;ck(call(s,o,c,r)==1234,"field failure propagated");fusion_prompt_result_v1 zero{};ck(std::memcmp(&r,&zero,sizeof(r))==0,"error clears result");
 c={};c.segment_error=1235;ck(call(s,o,c,r)==1235&&std::memcmp(&r,&zero,sizeof(r))==0,"segment failure clears result");
 c={};c.false_clear=true;o.initial_dt_s=1;o.horizon_s=1;ck(call(s,o,c,r)==FUSION_PROMPT_CALLBACK_CONTRACT,"false clear rejected");
 c={};c.malformed=true;ck(call(s,o,c,r)==FUSION_PROMPT_CALLBACK_CONTRACT,"malformed point rejected");
 c={};o=options();o.max_steps_per_level=1;ck(call(s,o,c,r)==0&&r.reason==FUSION_PROMPT_STEP_BUDGET&&r.outcome==0,"resource budget");
 o=options();o.max_levels=o.min_levels=2;o.event_time_tolerance_s=o.event_position_tolerance_m=1e-14;c={};c.B=1;ck(call(s,o,c,r)==0&&r.outcome==0&&r.reason==FUSION_PROMPT_REFINEMENT_LIMIT,"strict tolerance remains unresolved");
 o=options();o.horizon_s=std::numeric_limits<double>::max();o.initial_dt_s=std::numeric_limits<double>::denorm_min();c={};ck(call(s,o,c,r)==0&&r.reason==FUSION_PROMPT_STEP_BUDGET&&c.fields==0,"extreme ratio no integer overflow");
 o=options();o.horizon_s=o.initial_dt_s=std::numeric_limits<double>::denorm_min();c={};ck(call(s,o,c,r)==0&&r.reason==FUSION_PROMPT_TIME_RESOLUTION&&c.fields==0,"subnormal time no zero drift");
 o=options();o.horizon_s=std::numeric_limits<double>::quiet_NaN();ck(call(s,o,c,r)==PB11_STATUS_OUT_OF_RANGE,"invalid options");
 std::cout<<"prompt orbit analytic/domain/failure/budget tests pass\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
