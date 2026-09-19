#include "fusion_prompt_orbit.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
namespace {
using L=long double;
L distance(const double*a,const double*b){return std::hypot(std::hypot(L(a[0])-b[0],L(a[1])-b[1]),L(a[2])-b[2]);}
L norm(const double*a){const double z[3]={};return distance(a,z);}
bool finite(const double*x){for(int j=0;j<3;++j)if(!std::isfinite(x[j]))return false;return true;}
bool positive(double x){return std::isfinite(x)&&x>0;}
int classify(fusion_orbit_point_callback_v1 f,void*c,const double*x,int&kind){
 kind=-1;int rc=f(c,x,&kind);if(rc)return rc;
 return kind>=FUSION_BOUNDARY_AMBIGUOUS&&kind<=FUSION_BOUNDARY_OUTSIDE?0:FUSION_PROMPT_CALLBACK_CONTRACT;
}
struct Level {fusion_prompt_result_v1 r{};};
int drift(const double*a,const double*b,L t0,double half,int phase,
 const fusion_prompt_options_v1&o,fusion_orbit_point_callback_v1 point,
 fusion_orbit_segment_callback_v1 segment,void*c,Level&lv,bool&clear){
 clear=false;fusion_boundary_segment_value_v1 ev{};
 double tol=double(std::min(1.L,L(o.geometry_time_budget_s)/(4*L(half))));
 if(tol==0){lv.r.reason=FUSION_PROMPT_TIME_RESOLUTION;return 0;}
 int rc=segment(c,a,b,tol,&ev);if(rc)return rc;
 if(ev.result<FUSION_SEGMENT_UNRESOLVED||ev.result>FUSION_SEGMENT_GEOMETRIC_CANDIDATE)return FUSION_PROMPT_CALLBACK_CONTRACT;
 if(ev.result==FUSION_SEGMENT_UNRESOLVED){lv.r.reason=FUSION_PROMPT_BOUNDARY_UNRESOLVED;return 0;}
 int ka,kb;
 if(ev.result==FUSION_SEGMENT_CLEAR){
  if((rc=classify(point,c,a,ka))||(rc=classify(point,c,b,kb)))return rc;
  if(ka!=FUSION_BOUNDARY_INSIDE||kb!=FUSION_BOUNDARY_INSIDE)return FUSION_PROMPT_CALLBACK_CONTRACT;
  clear=true;return 0;
 }
 if(!std::isfinite(ev.fraction_lo)||!std::isfinite(ev.fraction_hi)||ev.fraction_lo<0||ev.fraction_hi>1||ev.fraction_lo>=ev.fraction_hi)return FUSION_PROMPT_CALLBACK_CONTRACT;
 double lo[3],hi[3];for(int j=0;j<3;++j){lo[j]=double((1-L(ev.fraction_lo))*a[j]+L(ev.fraction_lo)*b[j]);hi[j]=double((1-L(ev.fraction_hi))*a[j]+L(ev.fraction_hi)*b[j]);}
 if((rc=classify(point,c,lo,ka))||(rc=classify(point,c,hi,kb)))return rc;
 if(ka!=FUSION_BOUNDARY_INSIDE||kb!=FUSION_BOUNDARY_OUTSIDE)return FUSION_PROMPT_CALLBACK_CONTRACT;
 if((L(ev.fraction_hi)-ev.fraction_lo)*half>o.geometry_time_budget_s){lv.r.reason=FUSION_PROMPT_BOUNDARY_UNRESOLVED;return 0;}
 lv.r.event_time_lo_s=double(t0+L(ev.fraction_lo)*half);lv.r.event_time_hi_s=double(t0+L(ev.fraction_hi)*half);
 if(lv.r.event_time_hi_s<=lv.r.event_time_lo_s){lv.r={};lv.r.reason=FUSION_PROMPT_TIME_RESOLUTION;return 0;}
 if(L(lv.r.event_time_hi_s)-lv.r.event_time_lo_s>o.geometry_time_budget_s){lv.r={};lv.r.reason=FUSION_PROMPT_BOUNDARY_UNRESOLVED;return 0;}
 for(int j=0;j<3;++j){lv.r.event_x_lo_m[j]=lo[j];lv.r.event_x_hi_m[j]=hi[j];}
 lv.r.outcome=FUSION_PROMPT_EVENT;lv.r.phase=phase;return 0;
}
int run(double mass,double charge,const fusion_orbit_state_v1&initial,
 const fusion_prompt_options_v1&o,std::uint64_t n,fusion_magnetic_field_callback_v1 field,
 fusion_orbit_point_callback_v1 point,fusion_orbit_segment_callback_v1 segment,void*c,Level&lv){
 if(n>std::uint64_t(o.max_steps_per_level)){lv.r.reason=FUSION_PROMPT_STEP_BUDGET;return 0;}
 fusion_orbit_state_v1 state=initial,next{};
 for(std::uint64_t i=0;i<n;++i){
  // Same exact requested horizon at every level; adjacent rounded times define dt.
  double t0=double(L(o.horizon_s)*i/n),t1=double(L(o.horizon_s)*(i+1)/n),dt=t1-t0;
  if(dt<=0||dt*.5==0){lv.r.reason=FUSION_PROMPT_TIME_RESOLUTION;return 0;}
  double mid[3];int rc=fusion_c_magnetic_midpoint(&state,dt,mid);if(rc)return rc;
  bool clear=false;rc=drift(state.x_m,mid,t0,dt*.5,1,o,point,segment,c,lv,clear);
  lv.r.steps_last_level=int(i);if(rc||!clear)return rc;
  rc=fusion_c_magnetic_push(mass,charge,dt,&state,field,c,&next,mid);if(rc)return rc;
  rc=drift(mid,next.x_m,L(t0)+L(dt)*.5,dt*.5,2,o,point,segment,c,lv,clear);
  if(rc||!clear)return rc;
  state=next;
 }
 lv.r.outcome=FUSION_PROMPT_RETAINED;lv.r.final_state=state;lv.r.final_time_s=o.horizon_s;lv.r.steps_last_level=int(n);return 0;
}
}
extern "C" int fusion_c_prompt_orbit(double mass,double charge,
 const fusion_orbit_state_v1*initial,const fusion_prompt_options_v1*o,
 fusion_magnetic_field_callback_v1 field,fusion_orbit_point_callback_v1 point,
 fusion_orbit_segment_callback_v1 segment,void*c,fusion_prompt_result_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(!initial||!o||!field||!point||!segment||!std::isfinite(mass)||!std::isfinite(charge)||!finite(initial->x_m)||!finite(initial->u_m_s))return PB11_STATUS_INVALID_ARGUMENT;
 if(mass<=0||!positive(o->horizon_s)||!positive(o->initial_dt_s)||!positive(o->geometry_time_budget_s)||!positive(o->event_time_tolerance_s)||!positive(o->event_position_tolerance_m)||!positive(o->final_position_tolerance_m)||!positive(o->final_u_relative_tolerance)||o->min_levels<2||o->max_levels<o->min_levels||o->max_levels>20||o->max_steps_per_level<1||o->max_steps_per_level>1000000)return PB11_STATUS_OUT_OF_RANGE;
 int kind;int rc=classify(point,c,initial->x_m,kind);if(rc)return rc;
 if(kind==FUSION_BOUNDARY_OUTSIDE)return PB11_STATUS_OUT_OF_RANGE;
 if(kind==FUSION_BOUNDARY_AMBIGUOUS){out->reason=FUSION_PROMPT_BOUNDARY_UNRESOLVED;return 0;}
 L base=std::ceil(L(o->horizon_s)/o->initial_dt_s);
 if(base>o->max_steps_per_level){out->reason=FUSION_PROMPT_STEP_BUDGET;return 0;}
 // base is checked <=1e6 before conversion; doubling <=20 levels stays
 // below 2^40, so uint64_t cannot wrap before the per-level budget test.
 std::uint64_t n=std::uint64_t(std::max(1.L,base));Level previous,current;
 for(int level=1;level<=o->max_levels;++level,n*=2){
  current={};rc=run(mass,charge,*initial,*o,n,field,point,segment,c,current);if(rc)return rc;
  current.r.levels_used=level;
  if(level>=o->min_levels&&current.r.outcome&&current.r.outcome==previous.r.outcome){
   L dx=0,dt=0,du=0;
   if(current.r.outcome==FUSION_PROMPT_EVENT){
    // max cross-endpoint time difference equals center difference + half widths.
    dt=std::max(std::abs(L(current.r.event_time_hi_s)-previous.r.event_time_lo_s),std::abs(L(previous.r.event_time_hi_s)-current.r.event_time_lo_s));
    for(const double*a:{current.r.event_x_lo_m,current.r.event_x_hi_m})for(const double*b:{previous.r.event_x_lo_m,previous.r.event_x_hi_m})dx=std::max(dx,distance(a,b));
   }else{
    dx=distance(current.r.final_state.x_m,previous.r.final_state.x_m);
    L den=std::max(norm(current.r.final_state.u_m_s),norm(previous.r.final_state.u_m_s));
    du=den==0?0:distance(current.r.final_state.u_m_s,previous.r.final_state.u_m_s)/den;
   }
   bool passed=current.r.outcome==FUSION_PROMPT_EVENT?(dt<=o->event_time_tolerance_s&&dx<=o->event_position_tolerance_m):(dx<=o->final_position_tolerance_m&&du<=o->final_u_relative_tolerance);
   if(passed){current.r.estimated_time_spread_s=double(dt);current.r.estimated_position_spread_m=double(dx);current.r.estimated_u_relative_spread=double(du);*out=current.r;return 0;}
  }
  previous=current;
  if(current.r.reason==FUSION_PROMPT_STEP_BUDGET||current.r.reason==FUSION_PROMPT_TIME_RESOLUTION)break;
 }
 out->levels_used=current.r.levels_used;out->steps_last_level=current.r.steps_last_level;
 out->reason=current.r.reason?current.r.reason:FUSION_PROMPT_REFINEMENT_LIMIT;return 0;
}
