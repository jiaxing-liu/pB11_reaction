#include "fusion_pb_population.h"
#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>
namespace {
constexpr double u=1.602176634e-16;
void need(bool b,const char*msg){if(!b)throw std::runtime_error(msg);}
bool near(double a,double b,double tol=2e-8){return std::abs(a-b)<=tol*std::max(std::abs(a),std::abs(b));}
void budget(const fusion_pb_population_rates_v1&r){
 double fusion_beam_window_v1::*f[]={&fusion_beam_window_v1::resolved_reactivity_m3_s,&fusion_beam_window_v1::projectile_energy_reactivity_J_m3_s,&fusion_beam_window_v1::target_energy_reactivity_J_m3_s,&fusion_beam_window_v1::relative_energy_reactivity_J_m3_s,&fusion_beam_window_v1::cm_energy_reactivity_J_m3_s,&fusion_beam_window_v1::resolved_pair_probability};
 for(auto field:f){double sum=r.alpha0_peak.*field+r.narrow_remainder.*field+r.other_remainder.*field;need(near(sum,r.total.*field),"population sum");need(r.continuum_extrapolated.*field>=0 && r.continuum_extrapolated.*field<=r.total.*field*(1+1e-8),"diagnostic subset");}
 for(auto v:{r.total,r.alpha0_peak,r.narrow_remainder,r.other_remainder,r.continuum_extrapolated})need(near(v.projectile_energy_reactivity_J_m3_s+v.target_energy_reactivity_J_m3_s,v.relative_energy_reactivity_J_m3_s+v.cm_energy_reactivity_J_m3_s),"selected energy identity");
}
}
int main(){try{
 fusion_nuclear_mass_v1 p{},b{};need(!fusion_c_nuclear_mass(0,&p)&&!fusion_c_nuclear_mass(5,&b),"mass");
 for(int low:{0,1,2,3})for(double E:{0.,22.,148.,399.,401.,640.,2000.,10000.}){
  fusion_pb_population_v1 r{};need(!fusion_c_pb_population(1,low,E*u,.051,1,&r),"point status");
  double total=0;need(!fusion_c_cross_section_model(0,1,low,E*u,&total)&&total==r.total_cross_section_m2,"total unchanged");
  need(near(r.alpha0_peak_fraction+r.narrow_remainder_fraction+r.other_remainder_fraction,1.,2e-15),"point sum");
  need(r.narrow_fit_cross_section_m2>=0&&r.narrow_fit_cross_section_m2<=total,"fit subset");
  if(E>400)need(r.narrow_fit_fraction==0,"fit-piece boundary");
  if(E==148){double c0=197+(low==2?-12:low==3?12:0),c1=low==1?.240:.269,c2=low==1?2.31e-4:2.54e-4;double n=1.82e4/(2.35*2.35);need(near(r.narrow_fit_fraction,n/(c0+c1*148+c2*148*148+n),1e-14),"narrow peak reference");}
 }
 {fusion_pb_population_v1 a{},d{};double E=565*u/(1+p.mass_kg/b.mass_kg);need(!fusion_c_pb_population(1,0,E,.051,1,&a),"table node");need(near(a.continuum_peak_fraction,2.53/(2.53+598),2e-14),"table reference");need(!fusion_c_pb_population(1,0,E,.051,2,&d)&&near(d.continuum_peak_fraction,2*a.continuum_peak_fraction,2e-15),"scale");}
 for(bool thermal:{false,true})for(bool reverse:{false,true})for(double T:{0.,10.,150.}){
  double ma=reverse?b.mass_kg:p.mass_kg,mb=reverse?p.mass_kg:b.mass_kg;double ta=T*u,tb=T==0?0:3*u;
  fusion_pb_population_rates_v1 r{};int status=thermal?fusion_c_pb_thermal_population_rates(1,0,ma,mb,ta,tb,.051,1,&r):fusion_c_pb_beam_population_rates(1,0,ma,mb,ta,tb,.051,1,&r);need(!status,"rate status");budget(r);
  fusion_rate_model_v1 old{};status=thermal?fusion_c_thermal_pair_maxwellian_model(0,1,0,ma,mb,ta,tb,&old):fusion_c_beam_maxwellian_model(0,1,0,ma,mb,ta,tb,&old);need(!status&&r.total.resolved_reactivity_m3_s==old.total.resolved_reactivity_m3_s,"old rate parity");
 }
 {double beam=170*u;fusion_pb_population_rates_v1 r{};fusion_pb_population_v1 f{};need(!fusion_c_pb_beam_population_rates(1,0,p.mass_kg,b.mass_kg,beam,0,.051,1,&r),"cold beam");double E=beam/(1+p.mass_kg/b.mass_kg);need(!fusion_c_pb_population(1,0,E,.051,1,&f),"cold point");need(near(r.alpha0_peak.resolved_reactivity_m3_s/r.total.resolved_reactivity_m3_s,f.alpha0_peak_fraction,1e-13),"cold conditional fraction");budget(r);}
 {fusion_pb_population_v1 r{};r.total_cross_section_m2=99;need(fusion_c_pb_population(1,0,-1,.051,1,&r)==PB11_STATUS_OUT_OF_RANGE&&r.total_cross_section_m2==0,"point clear");need(fusion_c_pb_population(1,0,0,1.1,1,&r)==PB11_STATUS_OUT_OF_RANGE,"peak guard");need(fusion_c_pb_population(1,0,0,.051,2.1,&r)==PB11_STATUS_OUT_OF_RANGE,"scale guard");need(fusion_c_pb_population(0,0,0,.051,1,&r)==PB11_STATUS_INVALID_ARGUMENT,"policy guard");fusion_pb_population_rates_v1 q{};q.total.resolved_reactivity_m3_s=99;need(fusion_c_pb_thermal_population_rates(1,0,2*p.mass_kg,b.mass_kg,u,u,.051,1,&q)==PB11_STATUS_OUT_OF_RANGE&&q.total.resolved_reactivity_m3_s==0,"mass and clear");}
 std::cout<<"Population weighting tests passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
