#include "fusion_birth_table.h"
#include "fusion_nuclear_data.h"
#include <vector>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <chrono>
#include <memory>
void ok(int s,const char*m){if(s)throw std::runtime_error(std::string(m)+" status="+std::to_string(s));}
using Clock=std::chrono::steady_clock;using R=long double;
int main(int argc,char**argv){try{
 int ch=argc>1?std::stoi(argv[1]):3,n=argc>2?std::stoi(argv[2]):400;
 double lower=argc>3?std::stod(argv[3]):19.5,upper=argc>4?std::stod(argv[4]):23;
 double tol=argc>5?std::stod(argv[5]):.001;
 if(n<100||n>10000)throw std::runtime_error("cells");
 constexpr double keV=1.602176634e-16,MeV=1.602176634e-13;
 std::vector<double> edges(n+1),grid(7*n),reference(7*n);
 int low=n/3;for(int i=1;i<=low;++i)edges[i]=1e-18*keV*std::pow(1e20,double(i-1)/(low-1));
 for(int i=low+1;i<=n;++i)edges[i]=(100.+24900.*(i-low)/(n-low))*keV;
 fusion_thermal_birth_options_v1 source{2.5*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,8,8};
 fusion_birth_table_control_v1 control{tol,tol,tol,tol,1e-5,1e-5,128,1024,16};fusion_birth_table_v1*raw=nullptr;
 auto start=Clock::now();ok(fusion_c_birth_table_create(ch,lower*keV,upper*keV,&source,&control,n,edges.data(),&raw),"create");
 std::unique_ptr<fusion_birth_table_v1,decltype(&fusion_c_birth_table_destroy)> table(raw,fusion_c_birth_table_destroy);
 double seconds=std::chrono::duration<double>(Clock::now()-start).count();fusion_birth_table_info_v1 info{};ok(fusion_c_birth_table_info(raw,&info),"info");
 std::fprintf(stderr,"channel=%d cells=%d range_keV=%.9g:%.9g tolerance=%.9g knots=%d direct_evaluations=%d build_seconds=%.9g accepted_max_rate=%.9g debit=%.9g number_L1=%.9g energy_L1=%.9g\n",ch,n,lower,upper,tol,info.knots,info.direct_evaluations,seconds,info.max_validated_rate_error,info.max_validated_debit_error,info.max_validated_number_L1,info.max_validated_energy_L1);
 std::puts("T_keV,rate_error,debit_error,max_number_L1,max_energy_L1,direct_seconds,table_seconds");
 for(double fraction:{.137,.371,.613,.893}){
  double T=std::exp((1-fraction)*std::log(lower)+fraction*std::log(upper))*keV;
  fusion_thermal_birth_v1 ref{};start=Clock::now();ok(fusion_c_thermal_birth_grid(ch,T,&source,n,edges.data(),reference.data(),&ref),"direct");double direct_s=std::chrono::duration<double>(Clock::now()-start).count();
  fusion_birth_coefficients_v1 c{};start=Clock::now();for(int repeat=0;repeat<1000;++repeat)ok(fusion_c_birth_table_evaluate(raw,T,n,grid.data(),&c),"evaluate");double table_s=std::chrono::duration<double>(Clock::now()-start).count()/1000;
  auto rel=[](R diff,R denom){return denom>0?diff/denom:(diff==0?0:INFINITY);};
  R er=rel(std::abs(R(c.reactivity_m3_s)-ref.reactivity_m3_s),ref.reactivity_m3_s),ed=0,en=0,ee=0;
  for(int i=0;i<2;++i)ed=std::max(ed,rel(std::abs(R(c.reactant_energy_moment_J_m3_s[i])-ref.reactant_energy_moment_J_m3_s[i]),ref.reactant_energy_moment_J_m3_s[i]));
  for(int i=0;i<7;++i){R dn=std::abs(R(c.below_number_m3_s[i])-ref.below_number_m3_s[i])+std::abs(R(c.above_number_m3_s[i])-ref.above_number_m3_s[i]);R de=std::abs(R(c.below_energy_J_m3_s[i])-ref.below_energy_J_m3_s[i])+std::abs(R(c.above_energy_J_m3_s[i])-ref.above_energy_J_m3_s[i]);R N=R(ref.below_number_m3_s[i])+ref.above_number_m3_s[i],E=R(ref.below_energy_J_m3_s[i])+ref.above_energy_J_m3_s[i];
   for(int j=0;j<n;++j){int k=i*n+j;R center=(R(edges[j])+edges[j+1])/2;dn+=std::abs(R(grid[k])-reference[k]);de+=center*std::abs(R(grid[k])-reference[k]);N+=reference[k];E+=center*reference[k];}en=std::max(en,rel(dn,N));ee=std::max(ee,rel(de,E));}
  std::printf("%.12g,%.12Lg,%.12Lg,%.12Lg,%.12Lg,%.9g,%.9g\n",T/keV,er,ed,en,ee,direct_s,table_s);
  if(std::max({er,ed,en,ee})>tol)throw std::runtime_error("off-construction validation exceeds declared tolerance");
 }
 std::fprintf(stderr,"off-construction checks passed\n");return 0;
}catch(const std::exception&e){std::fprintf(stderr,"ERROR %s\n",e.what());return 1;}}
