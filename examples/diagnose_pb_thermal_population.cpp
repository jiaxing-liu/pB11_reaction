// Independent public-API diagnostic. Input: label kTa_J kTb_J per line.
// Source-conditional fractions are undefined when the total rate is zero.
#include "fusion_pb_population.h"
#include "fusion_nuclear_data.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
int main(int argc,char**argv){try{
 if(argc!=5)throw std::runtime_error("usage: diagnostic continuation pb_low narrow_peak continuum_scale < states.txt");
 auto number=[](const char*s){std::size_t n;double v=std::stod(s,&n);if(n!=std::string(s).size()||!std::isfinite(v))throw std::runtime_error("invalid option");return v;};
 double cp=number(argv[1]),lp=number(argv[2]);
 if(cp!=std::floor(cp)||cp<1||cp>2||lp!=std::floor(lp)||lp<0||lp>3)throw std::runtime_error("invalid policies");
 int continuation=int(cp),low=int(lp);double peak=number(argv[3]),scale=number(argv[4]);
 fusion_nuclear_mass_v1 p{},b{};
 if(fusion_c_nuclear_mass(0,&p)||fusion_c_nuclear_mass(5,&b))throw std::runtime_error("canonical masses");
 std::cout<<std::setprecision(17)<<"label,kTa_J,kTb_J,total_reactivity_m3_s,alpha0_peak_reactivity_m3_s,narrow_remainder_reactivity_m3_s,other_remainder_reactivity_m3_s,continuum_extrapolated_reactivity_m3_s,conditional_defined,mean_relative_energy_J,alpha0_peak_event_fraction,continuum_extrapolated_event_fraction\n";
 std::string line;int count=0;
 while(std::getline(std::cin,line)){
  if(line.empty())continue;
  std::istringstream in(line);std::string label,extra;double a,t;
  if(!(in>>label>>a>>t)||(in>>extra)||label.find_first_of(",\"\r\n")!=std::string::npos)throw std::runtime_error("malformed state");
  fusion_pb_population_rates_v1 r{};
  int rc=fusion_c_pb_thermal_population_rates(continuation,low,p.mass_kg,b.mass_kg,a,t,peak,scale,&r);
  if(rc)throw std::runtime_error("population API failed for "+label+": "+std::to_string(rc));
  double rate=r.total.resolved_reactivity_m3_s;
  std::cout<<label<<','<<a<<','<<t<<','<<rate<<','<<r.alpha0_peak.resolved_reactivity_m3_s<<','<<r.narrow_remainder.resolved_reactivity_m3_s<<','<<r.other_remainder.resolved_reactivity_m3_s<<','<<r.continuum_extrapolated.resolved_reactivity_m3_s<<','<<(rate>0);
  if(rate>0)std::cout<<','<<r.total.relative_energy_reactivity_J_m3_s/rate<<','<<r.alpha0_peak.resolved_reactivity_m3_s/rate<<','<<r.continuum_extrapolated.resolved_reactivity_m3_s/rate;
  else std::cout<<",,,";
  std::cout<<'\n';++count;
 }
 if(!std::cin.eof()||count==0)throw std::runtime_error("empty or unreadable input");
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
