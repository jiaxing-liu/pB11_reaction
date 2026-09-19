#include "fusion_thermal_markers.h"
#include "fusion_prompt_reduce.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
using L=long double;
constexpr double pi=3.1415926535897932384626433832795;
void ck(bool b,const char*s){if(!b)throw std::runtime_error(s);}
double norm(const double*x){return std::hypot(x[0],x[1],x[2]);}
int point(void*,const double*x,int*k){double r=norm(x);*k=r<1?FUSION_BOUNDARY_INSIDE:r>1?FUSION_BOUNDARY_OUTSIDE:FUSION_BOUNDARY_AMBIGUOUS;return 0;}
int field(void*ctx,const double*x,double*b){if(norm(x)>=1)return 999;b[0]=b[1]=0;b[2]=*static_cast<double*>(ctx);return 0;}
int segment(void*,const double*a,const double*b,double tol,fusion_boundary_segment_value_v1*r){
 *r={};if(norm(a)>=1)return 998;
 // A sphere is convex, so inside endpoints prove the entire chord inside.
 if(norm(b)<1){r->result=FUSION_SEGMENT_CLEAR;return 0;}
 double lo=0,hi=1;
 for(int i=0;i<64&&hi-lo>tol*.25;++i){double m=(lo+hi)/2,x[3];if(m==lo||m==hi)break;for(int k=0;k<3;++k)x[k]=a[k]+m*(b[k]-a[k]);if(norm(x)<1)lo=m;else hi=m;}
 hi=std::min(1.,hi+tol*.125);double x[3];for(int k=0;k<3;++k)x[k]=a[k]+hi*(b[k]-a[k]);
 if(norm(x)<=1)return 0;
 r->result=FUSION_SEGMENT_GEOMETRIC_CANDIDATE;r->fraction_lo=lo;r->fraction_hi=hi;r->tolerance_met=hi-lo<=tol;return 0;
}
struct Exact {bool event;L theta;};
// Planar uniform-B helix: g(theta)=A+B sin(theta)+C cos(theta).
// All derivative roots are atan2(B,C)+k*pi. Splitting at ALL of them prevents
// missing an exit followed by reentry, even when the full-period endpoints agree.
Exact first_exit(L radius,L rho,L phi,L horizon){
 L A=radius*radius+2*radius*rho*std::sin(phi)+2*rho*rho-1;
 L B=2*rho*radius*std::cos(phi),C=-2*rho*radius*std::sin(phi)-2*rho*rho;
 auto g=[&](L t){return A+B*std::sin(t)+C*std::cos(t);};
 std::vector<L> knots{0,horizon};L theta=std::atan2(B,C),p=std::acos(-1.L);
 for(int k=-4;k<=4;++k){L t=theta+k*p;if(t>0&&t<horizon)knots.push_back(t);}
 std::sort(knots.begin(),knots.end());
 for(size_t i=1;i<knots.size();++i){L a=knots[i-1],b=knots[i];if(g(a)<0&&g(b)>0){
   for(int it=0;it<100;++it){L m=(a+b)/2;if(m==a||m==b)break;if(g(m)<0)a=m;else b=m;}
   return {true,(a+b)/2};
 }}return {false,horizon};
}
void analytic_position(L radius,L rho,L phi,L theta,double*x){
 x[0]=double(radius+rho*(std::sin(phi)-std::sin(phi-theta)));
 x[1]=double(rho*(std::cos(phi-theta)-std::cos(phi)));x[2]=0;
}
int main(){try{
 constexpr double mass=1.67262192369e-27,q=1.602176634e-19,K=1.602176634e-15,c=299792458.,rho=.75;
 double u=std::sqrt(K/mass)*std::sqrt(2+K/(mass*c*c)),gamma=std::hypot(1.,u/c),v=u/gamma,omega=v/rho,B=gamma*mass*omega/q;
 auto known=first_exit(.5,rho,0,1);ck(known.event&&std::abs(known.theta-.701758821744515L)<1e-13,"canonical first root");
 auto return_case=first_exit(.5,rho,0,2*pi);ck(return_case.event&&std::abs(return_case.theta-known.theta)<1e-15,"first exit despite full-period reentry");
 std::cout<<std::setprecision(16)<<"horizon_phase,initial_dt_phase,event_count,retained_count,unresolved_count,max_event_phase_error,max_event_position_error,max_retained_position_error\n";
 for(double horizon:{1.,2*pi}){
  double prior=1;
  for(double phase:{.2,.1,.05,.025,.0125}){
   fusion_prompt_options_v1 o{horizon/omega,phase/omega,1e-11/omega,1./omega,1.,1.,1.,2,2,4096};
   double max_t=0,max_x=0,max_final=0;int events=0,retained=0,unresolved=0;
   std::vector<fusion_prompt_weight_v1>w;double exact_events=0;
   for(double radius:{.2,.5,.7})for(int k=0;k<4;++k){
    double phi=k*pi/2;fusion_orbit_state_v1 initial{{radius,0,0},{u*std::cos(phi),u*std::sin(phi),0}};
    Exact ref=first_exit(radius,rho,phi,horizon);if(ref.event)exact_events+=1./12;
    fusion_prompt_result_v1 r{};int status=fusion_c_prompt_orbit(mass,q,&initial,&o,field,point,segment,&B,&r);
    ck(status==0,"uniform B driver error");
    w.push_back({1./12,K,status,r.outcome});
    if(r.outcome==FUSION_PROMPT_UNRESOLVED){++unresolved;continue;}
    ck((r.outcome==FUSION_PROMPT_EVENT)==ref.event,"uniform B exact classification");
    double x[3];analytic_position(radius,rho,phi,ref.theta,x);
    if(ref.event){++events;double t=(r.event_time_lo_s+r.event_time_hi_s)*omega/2;max_t=std::max(max_t,std::abs(t-double(ref.theta)));
      double dx[3];for(int j=0;j<3;++j)dx[j]=(r.event_x_lo_m[j]+r.event_x_hi_m[j])/2-x[j];max_x=std::max(max_x,norm(dx));
    }else{++retained;double dx[3];for(int j=0;j<3;++j)dx[j]=r.final_state.x_m[j]-x[j];max_final=std::max(max_final,norm(dx));}
   }
   fusion_prompt_totals_v1 totals{};ck(fusion_c_prompt_reduce(int(w.size()),w.data(),&totals)==0,"uniform ensemble reduce");
   ck(exact_events>=totals.number_fraction_low-1e-14&&exact_events<=totals.number_fraction_high+1e-14,"uniform reference weight bounds");
   ck(std::abs(totals.total_number-1)<1e-14&&std::abs(totals.total_energy_J/K-1)<1e-14,"uniform ensemble closure");
   ck(unresolved==0,"non-grazing fixed set resolved");
   double error=std::max({max_t,max_x,max_final});ck(error<prior,"uniform refinement error decreases");prior=error;
   std::cout<<horizon<<','<<phase<<','<<events<<','<<retained<<','<<unresolved<<','<<max_t<<','<<max_x<<','<<max_final<<'\n';
  }
  ck(prior<2e-5,"finest absolute orbit accuracy");
 }
 // Independently known non-crossing orbit with maximum radius 1-2e-8.
 // Near contact may remain unresolved; it must not be asserted to escape.
 {
  double radius=.25,small_rho=.375-1e-8,phi=pi/2,om=v/small_rho,field_B=gamma*mass*om/q;
  ck(!first_exit(radius,small_rho,phi,2*pi).event,"near tangent reference retained");
  fusion_orbit_state_v1 initial{{radius,0,0},{u*std::cos(phi),u*std::sin(phi),0}};
  fusion_prompt_options_v1 o{2*pi/om,.1/om,1e-12/om,1e-10/om,1e-10,1e-10,1e-11,3,12,20000};
  fusion_prompt_result_v1 result{};
  ck(fusion_c_prompt_orbit(mass,q,&initial,&o,field,point,segment,&field_B,&result)==0,"near tangent API");
  ck(result.outcome!=FUSION_PROMPT_EVENT,"non-crossing near tangent falsely escaped");
  std::cout<<"# near_tangent outcome="<<result.outcome<<" reason="<<result.reason<<" levels="<<result.levels_used<<'\n';
 }
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
