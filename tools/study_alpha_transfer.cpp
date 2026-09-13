// Absorbing-boundary sensitivity against the full resolved trace-alpha FP.
// This diagnostic does NOT by itself label an energy cutoff physical ash.
#include "fusion_coulomb.h"
#include "fusion_kinetics.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
constexpr double kev=1.602176634e-16;
void run(int n,int steps,double cut,double penalty,int grid=0){
 const double pop=1e12,initial=3000*kev,T=10*kev,dt=1./steps,ma=6.644657345e-27;
 fusion_maxwellian_bath_v1 baths[3]={{1e20,9.1093837139e-31,1,T,15},{5e19,3.3435837768e-27,1,T,15},{5e19,5.0073567512e-27,1,T,15}};
 double ts[]={T,T,T};std::vector<double>e(n+1),c(n),old(n),next(n),zero(n),D(3*(n-1)),eq(n);
 const double lo=cut>0?cut:.001,hi=10000.;
 for(int i=0;i<=n;++i)e[i]=lo*kev*std::pow(hi/lo,double(i)/n);
 if(grid==1){
  const double join=std::max(20.,lo);const int split=lo<20?n/4:0;
  for(int i=0;i<=n;++i)e[i]=kev*(i<split?lo*std::pow(join/lo,double(i)/split):join+(6000.-join)*double(i-split)/(n-split));
 }
 double norm=0;for(int i=0;i<n;++i){c[i]=(e[i]+e[i+1])/2;eq[i]=(e[i+1]-e[i])*std::sqrt(c[i])*std::exp(-c[i]/T);norm+=eq[i];}
 for(double &v:eq)v/=norm;
 int r=int(std::upper_bound(c.begin(),c.end(),initial)-c.begin());
 double w=(initial-c[r-1])/(c[r]-c[r-1]);old[r]=pop*w;old[r-1]=pop*(1-w);
 for(int b=0;b<3;++b)for(int i=0;i<n-1;++i){fusion_coulomb_energy_v1 o{};
  if(fusion_c_coulomb_energy(e[i+1],ma,2,&baths[b],&o))throw std::runtime_error("coefficient");D[b*(n-1)+i]=o.diffusion_J2_s;}
 double H[3]={},ash=0,ashE=0;
 for(int s=1;s<=steps;++s){fusion_kinetic_ledger_v1 l{};double h[3];
  int status=fusion_c_energy_fp_trial(n,3,dt,e.data(),old.data(),ts,D.data(),zero.data(),zero.data(),cut>0?penalty/dt:0,next.data(),h,&l);
  if(status){std::cerr<<n<<","<<steps<<","<<cut<<","<<s<<","<<status<<'\n';throw std::runtime_error("FP rejected");}
  for(int b=0;b<3;++b)H[b]+=h[b];ash+=l.thermalized_number_m3;ashE+=l.thermalized_energy_J_m3;old.swap(next);
  if(s%std::max(1,steps/200)==0){double N=0,E=0,overlap=0,below=0;
   for(int i=0;i<n;++i){N+=old[i];E+=old[i]*c[i];overlap+=std::min(old[i]/pop,eq[i]);if(c[i]<=3*T)below+=old[i];}
   double ne=(N+ash)/pop-1,ee=(E+ashE+H[0]+H[1]+H[2])/(pop*initial)-1;
   if(std::abs(ne)>1e-9 || std::abs(ee)>1e-9)throw std::runtime_error("cumulative conservation");
   // If the absorbed particles are mixed into a Ti Maxwellian, the remaining
   // carried energy minus 3Ti/2 must also be exchanged with the ion pool.
   double ionMixed=H[1]+H[2]+ashE-1.5*T*ash;
   std::cout<<grid<<','<<n<<','<<steps<<','<<cut<<','<<c[0]/kev<<','<<penalty<<','<<s*dt<<','<<N/pop<<','<<ash/pop<<','<<E/(pop*initial)<<','<<ashE/(pop*initial)<<','<<H[0]/(pop*initial)<<','<<ionMixed/(pop*initial)<<','<<below/pop<<','<<overlap<<','<<ne<<','<<ee<<'\n';
  }
 }
}
int main(int argc,char**argv){try{
 std::cout<<std::setprecision(16)<<"grid,cells,steps,cut_keV,sink_energy_keV,penalty,time_s,kinetic_N_fraction,absorbed_N_fraction,kinetic_energy_fraction,carried_energy_fraction,electron_heat_fraction,ion_heat_after_Ti_mixing_fraction,kinetic_below_3Ti_fraction,truncated_Maxwellian_overlap,number_residual,energy_residual\n";
 if(argc==5){run(std::stoi(argv[1]),std::stoi(argv[2]),std::stod(argv[3]),1e6,std::stoi(argv[4]));return 0;}
 for(double cut:{0.,1.,2.5,5.,10.,20.,40.,80.})run(400,1000,cut,1e6);
 for(double cut:{0.,2.5,10.,40.})run(800,2000,cut,1e6);
 for(double penalty:{1e4,1e8})run(400,1000,10.,penalty);
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
