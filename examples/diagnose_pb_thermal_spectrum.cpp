// Independent public-API, whole scalar marginal sensitivity diagnostic.
// No density factor, sigma modification, event correlations, or host mutation.
// stdin: cells; edges[cells+1]; 8 double + 9 integer birth options in struct
// order; then label kT_J pairs. Baseline uses input nq/ncos; refinement doubles
// each (capped at1024). Five explicit remainder/continuum alternatives.
#include "fusion_thermal_birth.h"
#include "fusion_nuclear_data.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
struct State {std::string label;double T;};
std::string token(){std::string s;if(!(std::cin>>s))throw std::runtime_error("truncated input");return s;}
double number(const std::string&s){std::size_t n=0;double x=std::stod(s,&n);if(n!=s.size()||!std::isfinite(x))throw std::runtime_error("invalid finite number: "+s);return x;}
int integer(){double x=number(token());if(x!=std::floor(x)||x<std::numeric_limits<int>::min()||x>std::numeric_limits<int>::max())throw std::runtime_error("invalid integer");return int(x);}
void options(std::ostream&f,const fusion_thermal_birth_options_v1&o){
 f<<','<<o.relative_max_J<<','<<o.cm_max_kT<<','<<o.ground_state_q_J<<','<<o.cutoff_J<<','<<o.l1_fraction<<','<<o.relative_phase<<','<<o.narrow_peak_fraction<<','<<o.continuum_peak_scale<<','<<o.continuation<<','<<o.pb_low<<','<<o.remainder_policy<<','<<o.broad_mode<<','<<o.fsci_policy<<','<<o.relative_order<<','<<o.cm_order<<','<<o.nq<<','<<o.ncos;
}
}
int main(int argc,char**argv){try{
 if(argc!=3&&argc!=4)throw std::runtime_error("usage: diagnostic spectrum.csv metadata.csv [base-only] < input.txt");
 const bool base_only=argc==4;
 if(base_only&&std::string(argv[3])!="base-only")throw std::runtime_error("unknown variant selection");
 if(std::string(argv[1])==argv[2]||std::ifstream(argv[1]).good()||std::ifstream(argv[2]).good())throw std::runtime_error("output path exists or paths coincide");
 const int cells=integer();if(cells<1||cells>100000)throw std::runtime_error("cells outside1..100000");
 std::vector<double> edges(cells+1);for(double&x:edges)x=number(token());
 if(edges.front()<0)throw std::runtime_error("negative edge");
 for(int i=1;i<=cells;++i)if(!(edges[i]>edges[i-1]))throw std::runtime_error("nonincreasing edges");
 fusion_thermal_birth_options_v1 base{};
 double*ds[]={&base.relative_max_J,&base.cm_max_kT,&base.ground_state_q_J,&base.cutoff_J,&base.l1_fraction,&base.relative_phase,&base.narrow_peak_fraction,&base.continuum_peak_scale};for(auto p:ds)*p=number(token());
 int*is[]={&base.continuation,&base.pb_low,&base.remainder_policy,&base.broad_mode,&base.fsci_policy,&base.relative_order,&base.cm_order,&base.nq,&base.ncos};for(auto p:is)*p=integer();
 // Variants intentionally replace remainder and continuum scale only. Refuse
 // inconsistent base settings rather than quietly calling a different base.
 if(base.remainder_policy!=0||base.continuum_peak_scale!=1)throw std::runtime_error("base requires remainder0 scale1");
 if(base.nq<4||base.nq>1024||base.ncos<4||base.ncos>1024)throw std::runtime_error("invalid angular orders");
 std::vector<State> states;std::set<std::string> labels;std::string label;
 while(std::cin>>label){if(label.find_first_of(",\"\r\n")!=std::string::npos||!labels.insert(label).second)throw std::runtime_error("unsafe or duplicate label");double T=number(token());if(T<=0)throw std::runtime_error("nonpositive temperature");states.push_back({label,T});}
 if(!std::cin.eof()||states.empty())throw std::runtime_error("empty or unreadable states");
 fusion_nuclear_channel_v1 reaction{};int rc=fusion_c_nuclear_channel(0,&reaction);if(rc)throw std::runtime_error("canonical channel lookup");
 std::ofstream spectrum(argv[1]),meta(argv[2]);if(!spectrum||!meta)throw std::runtime_error("cannot open outputs");spectrum<<std::setprecision(17);meta<<std::setprecision(17);
 spectrum<<"label,kT_J,variant,level,nq,ncos,species,cell,lower_J,upper_J,center_J,birth_number_m3_s\n";
 meta<<"label,kT_J,variant,level,status,input_nq,input_ncos,relative_max_J,cm_max_kT,ground_state_q_J,cutoff_J,l1_fraction,relative_phase,narrow_peak_fraction,continuum_peak_scale,continuation,pb_low,remainder_policy,broad_mode,fsci_policy,relative_order,cm_order,nq,ncos,rate,reference_rate,debit0,debit1,reference_debit0,reference_debit1,product_energy,number_residual,energy_residual,relative_rate_discrepancy,relative_debit_discrepancy,cm_retained_probability,cm_tail_probability,cm_tail_energy,max_cm_shift,max_shell_remap,mapped_charged_N,total_charged_N,total_product_E,expected_charged_N,expected_product_E,N_closure,E_closure,product_E_closure";
 for(int sp=0;sp<7;++sp){meta<<",below_N_"<<sp<<",below_E_"<<sp<<",above_N_"<<sp<<",above_E_"<<sp;}meta<<'\n';
 const char*names[]={"base","all_low","all_broad","continuum_zero","continuum_double"};
 for(const auto&s:states)for(int variant=0;variant<(base_only?1:5);++variant)for(int level=0;level<2;++level){
  auto o=base;if(variant==1)o.remainder_policy=1;if(variant==2)o.remainder_policy=2;if(variant==3)o.continuum_peak_scale=0;if(variant==4)o.continuum_peak_scale=2;
  if(level){o.nq=std::min(1024,base.nq*2);o.ncos=std::min(1024,base.ncos*2);}
  std::vector<double> birth(7*cells);fusion_thermal_birth_v1 r{};rc=fusion_c_thermal_birth_grid(0,s.T,&o,cells,edges.data(),birth.data(),&r);
  meta<<s.label<<','<<s.T<<','<<names[variant]<<','<<level<<','<<rc<<','<<base.nq<<','<<base.ncos;options(meta,o);
  if(rc){for(int missing=0;missing<52;++missing)meta<<',';meta<<'\n';meta.flush();throw std::runtime_error("thermal birth API failed label="+s.label+" variant="+names[variant]+" level="+std::to_string(level)+" status="+std::to_string(rc));}
  long double mappedN=0,totalN=0,totalE=0;
  for(int sp=0;sp<7;++sp){for(int cell=0;cell<cells;++cell){const double b=birth[sp*cells+cell];long double center=(static_cast<long double>(edges[cell])+edges[cell+1])/2;
    spectrum<<s.label<<','<<s.T<<','<<names[variant]<<','<<level<<','<<o.nq<<','<<o.ncos<<','<<sp<<','<<cell<<','<<edges[cell]<<','<<edges[cell+1]<<','<<center<<','<<b<<'\n';
    if(sp<6){mappedN+=b;totalN+=b;}totalE+=b*center;
   }if(sp<6)totalN+=static_cast<long double>(r.below_number_m3_s[sp])+r.above_number_m3_s[sp];totalE+=static_cast<long double>(r.below_energy_J_m3_s[sp])+r.above_energy_J_m3_s[sp];}
  const long double expectedN=3*static_cast<long double>(r.reactivity_m3_s),expectedE=static_cast<long double>(reaction.q_J)*r.reactivity_m3_s+r.reactant_energy_moment_J_m3_s[0]+r.reactant_energy_moment_J_m3_s[1];
  meta<<','<<r.reactivity_m3_s<<','<<r.reference_reactivity_m3_s<<','<<r.reactant_energy_moment_J_m3_s[0]<<','<<r.reactant_energy_moment_J_m3_s[1]<<','<<r.reference_reactant_energy_moment_J_m3_s[0]<<','<<r.reference_reactant_energy_moment_J_m3_s[1]<<','<<r.product_energy_moment_J_m3_s<<','<<r.number_residual_m3_s<<','<<r.energy_residual_J_m3_s<<','<<r.relative_rate_discrepancy<<','<<r.relative_reactant_energy_discrepancy<<','<<r.cm_retained_probability<<','<<r.cm_tail_probability<<','<<r.cm_tail_energy_moment_J<<','<<r.max_cm_energy_shift_fraction<<','<<r.max_shell_remap_fraction<<','<<mappedN<<','<<totalN<<','<<totalE<<','<<expectedN<<','<<expectedE<<','<<totalN-expectedN<<','<<totalE-expectedE<<','<<totalE-r.product_energy_moment_J_m3_s;
  for(int sp=0;sp<7;++sp){meta<<','<<r.below_number_m3_s[sp]<<','<<r.below_energy_J_m3_s[sp]<<','<<r.above_number_m3_s[sp]<<','<<r.above_energy_J_m3_s[sp];}meta<<'\n';
  meta.flush();
  if(!spectrum||!meta)throw std::runtime_error("output write failed");
  std::cerr<<"completed label="<<s.label<<" variant="<<names[variant]<<" level="<<level<<" nq="<<o.nq<<" ncos="<<o.ncos<<std::endl;
 }
 spectrum.flush();meta.flush();if(!spectrum||!meta)throw std::runtime_error("output flush failed");spectrum.close();meta.close();if(!spectrum||!meta)throw std::runtime_error("output close failed");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
