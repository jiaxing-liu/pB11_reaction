#include "fusion_cross_sections.h"
#include "fusion_nuclear_data.h"
#include "fusion_beam.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace {
constexpr double kev=1.602176634e-16,barn=1e-28;
using GK=boost::math::quadrature::gauss_kronrod<double,61>;
void check(bool x){if(!x)throw std::runtime_error("window study check");}
struct Model{
 int ch; double lo,hi,slo,shi,bg;
 double sigma(double E,int mode)const{ // barns; explicit illustrative continuation
  if(E==0)return 0;
  if(mode>=4){
   check(ch==0 && E<=400);
   const double c1=mode==4?.240:.269,c2=mode==4?2.31e-4:2.54e-4;
   const double c0=197+(mode==5?12:(mode==6?-12:0));
   double S=c0+c1*E+c2*E*E+1.82e4/((E-148)*(E-148)+2.35*2.35);
   return 1000*S/E*std::exp(-bg/std::sqrt(E));
  }
  if(mode==0){double v=0;check(fusion_c_cross_section(ch,E*kev,&v)==0);return v/barn;}
  if(mode==3)return shi; // high constant-sigma sensitivity only
  double boundary=mode==1?lo:hi,sigma_boundary=mode==1?slo:shi;
  // Constant S at endpoint with the existing channel's fixed Gamow exponent.
  return sigma_boundary*std::exp(std::log(boundary/E)+bg*(1/std::sqrt(boundary)-1/std::sqrt(E)));
 }
 double integral(double T,double Elo,double Ehi,int mode,int moment,double tol)const{
  double a=std::max(0.,Elo/T),b=std::min(800.,Ehi/T);if(b<=a)return 0;
  std::vector<double>cuts{a,b};
  auto add=[&](double x){if(x>a&&x<b)cuts.push_back(x);};
  for(double x:{.01,.1,.25,.5,1.,2.,4.,8.,16.,32.,64.,128.,256.,512.})add(x);
  for(double E:{101.,148.,195.,400.,530.,640.9,668.,900.,1211.,2340.,3294.,3480.,5700.})add(E/T);
  std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
  long double sum=0;
  for(size_t j=1;j<cuts.size();++j){
   auto f=[&](double x){if(x==0)return 0.;double v=sigma(T*x,mode)*x*std::exp(-x);return moment?v*x:v;};
   double error=0;double value=GK::integrate(f,cuts[j-1],cuts[j],18,tol,&error);
   check(std::isfinite(value)&&value>=0&&std::isfinite(error));sum+=value;
  }
  return double(sum);
 }
};
}
int main(int argc,char**argv){
 double tolerance=argc>1?std::stod(argv[1]):1e-10;
 check(std::isfinite(tolerance)&&tolerance>0);
 std::cout<<std::setprecision(17)<<"channel,T_keV,K_fit_m3_s,K_low_Sconst_m3_s,K_high_Sconst_m3_s,K_high_flat_m3_s,relative_E_fit_keV,relative_E_low_keV,relative_E_high_keV,pb_K_below_SW140,pb_K_gap3480_5700,window_rate_relative_difference,window_moment_relative_difference,pb_K_below_Becker22,pb_K_low400_TB,pb_K_low400_NS,pb_K_low400_C0plus12,pb_K_low400_C0minus12\n";
 const double bg[5]={std::sqrt(22589.),31.3970,31.3970,34.3827,68.7508};
 for(int ch=0;ch<5;++ch){
  double low=0,high=0,slo=0,shi=0;check(fusion_c_cross_section_domain(ch,&low,&high)==0);
  if(low>0)check(fusion_c_cross_section(ch,low,&slo)==0);
  check(fusion_c_cross_section(ch,high,&shi)==0);
  Model m{ch,low/kev,high/kev,slo/barn,shi/barn,bg[ch]};
  fusion_nuclear_channel_v1 r{};fusion_nuclear_mass_v1 a{},b{};
  check(fusion_c_nuclear_channel(ch,&r)==0);
  check(fusion_c_nuclear_mass(r.reactant_ids[0],&a)==0&&fusion_c_nuclear_mass(r.reactant_ids[1],&b)==0);
  const double mu=a.mass_kg/(1+a.mass_kg/b.mass_kg);
  for(double T:{.01,.03,.05,.1,.2,.5,1.,3.,10.,30.,100.,190.,500.}){
   const double scale=std::sqrt(8*T*kev/(std::acos(-1.)*mu))*barn;
   double fit=m.integral(T,m.lo,m.hi,0,0,tolerance),fm=m.integral(T,m.lo,m.hi,0,1,tolerance);
   double lower=m.lo>0?m.integral(T,0,m.lo,1,0,tolerance):0;
   double lm=m.lo>0?m.integral(T,0,m.lo,1,1,tolerance):0;
   double upper=m.integral(T,m.hi,800*T,2,0,tolerance),um=m.integral(T,m.hi,800*T,2,1,tolerance);
   double flat=m.integral(T,m.hi,800*T,3,0,tolerance);
   double pblo=ch==0?m.integral(T,0,140,0,0,tolerance):0;
   double pbgap=ch==0?m.integral(T,3480,5700,0,0,tolerance):0;
   fusion_beam_window_v1 w{};check(fusion_c_thermal_pair_maxwellian_window(ch,a.mass_kg,b.mass_kg,T*kev,T*kev,&w)==0);
   double kd=fit>0?std::abs(w.resolved_reactivity_m3_s-scale*fit)/(scale*fit):0;
   double md=fm>0?std::abs(w.relative_energy_reactivity_J_m3_s-scale*T*kev*fm)/(scale*T*kev*fm):0;
   check(kd<2e-7&&md<2e-7);
   std::cout<<ch<<','<<T<<','<<scale*fit<<','<<scale*lower<<','<<scale*upper<<','<<scale*flat<<','
    <<(fit>0?T*fm/fit:0)<<','<<(lower>0?T*lm/lower:0)<<','<<(upper>0?T*um/upper:0)<<','
    <<scale*pblo<<','<<scale*pbgap<<','<<kd<<','<<md;
   double b22=0,tb=0,ns=0,plus=0,minus=0;
   if(ch==0){
    b22=m.integral(T,0,22,0,0,tolerance);
    tb=m.integral(T,0,400,0,0,tolerance);
    ns=m.integral(T,0,400,4,0,tolerance);
    plus=m.integral(T,0,400,5,0,tolerance);
    minus=m.integral(T,0,400,6,0,tolerance);
    check(std::abs((plus+minus)/2-tb)<=1e-10*tb);
   }
   std::cout<<','<<scale*b22<<','<<scale*tb<<','<<scale*ns<<','<<scale*plus<<','<<scale*minus<<'\n';
  }
 }
}
