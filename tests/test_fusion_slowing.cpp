// Trace-alpha relaxation benchmark, NOT a self-heated fusion or device case.
// Fixed Maxwellian electron/D/T baths at 10 keV, supplied lnLambda=15.
// No births, escape or ash conversion: equilibrium remains on the kinetic grid.
// Initial 3 MeV energy is split between adjacent centers preserving N and E.
// Tests independent accumulated budgets, the analytic 3kT/2 mean, and
// convergence to the sampled Maxwellian; also covers tiny-tail underflow.
#include "fusion_coulomb.h"
#include "fusion_kinetics.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <array>
constexpr double kev=1.602176634e-16;
struct Result { double electron_fraction, equilibrium_error; };
Result run(int n,int steps) {
 double duration=1., ma=6.644657345e-27, initial=3000*kev, pop=1e12;
 double ts[3]={10*kev,10*kev,10*kev};
 fusion_maxwellian_bath_v1 baths[3]={{1e20,9.1093837139e-31,1,ts[0],15},{5e19,3.3435837768e-27,1,ts[1],15},{5e19,5.0073567512e-27,1,ts[2],15}};
 std::vector<double> e(n+1),c(n),old(n),next(n),zero(n),d(3*(n-1));
 for(int i=0;i<=n;i++)e[i]=.001*kev*std::pow(1e7,double(i)/n);
 for(int i=0;i<n;i++)c[i]=.5*(e[i]+e[i+1]);
 int right=std::upper_bound(c.begin(),c.end(),initial)-c.begin();
 double w=(initial-c[right-1])/(c[right]-c[right-1]);old[right]=pop*w;old[right-1]=pop*(1-w);
 for(int b=0;b<3;b++)for(int i=0;i<n-1;i++){
 fusion_coulomb_energy_v1 o{};if(fusion_c_coulomb_energy(e[i+1],ma,2,&baths[b],&o)) throw std::runtime_error("Coulomb coefficient failed"); d[b*(n-1)+i]=o.diffusion_J2_s;
 }
 double heats[3]={},maxerr=0;
 for(int s=0;s<steps;s++){
 double h[3]; fusion_kinetic_ledger_v1 l{};
 int status=fusion_c_energy_fp_trial(n,3,duration/steps,e.data(),old.data(),ts,d.data(),zero.data(),zero.data(),0,next.data(),h,&l); if(status){std::cerr<<"failed "<<n<<" "<<steps<<" step "<<s<<" status "<<status<<"\n";throw std::runtime_error("trial rejected");}
 for(int b=0;b<3;b++) { heats[b]+=h[b]; }
 maxerr=std::max(maxerr,std::abs(l.energy_balance_error_J_m3)/(pop*initial));old.swap(next);
 }
 double energy=0,number=0,eqnorm=0,err=0;
 for(int i=0;i<n;i++){number+=old[i];energy+=c[i]*old[i];eqnorm+=(e[i+1]-e[i])*sqrt(c[i])*exp(-c[i]/ts[0]);}
 for(int i=0;i<n;i++)err+=fabs(old[i]-pop*(e[i+1]-e[i])*sqrt(c[i])*exp(-c[i]/ts[0])/eqnorm)/pop;
 if (std::abs(number/pop-1)>1e-10 ||
     std::abs((energy+heats[0]+heats[1]+heats[2]-pop*initial)/(pop*initial))>1e-10)
     throw std::runtime_error("independent accumulated budget failed");
 if (std::abs(energy/pop/kev-15)>0.002 || err>2e-5)
     throw std::runtime_error("common-temperature Maxwellian not approached");
 if (heats[0]<=0 || heats[1]<=0 || heats[2]<=0)
     throw std::runtime_error("3 MeV alpha did not heat each cold bath");
 std::cout<<n<<','<<steps<<','<<energy/pop/kev<<','<<heats[0]/(pop*initial)<<','<<(heats[1]+heats[2])/(pop*initial)<<','<<number/pop-1<<','<<(energy+heats[0]+heats[1]+heats[2]-pop*initial)/(pop*initial)<<','<<err<<','<<maxerr<<'\n';
 return {heats[0]/(pop*initial),err};
}
int main() {
 try {
 std::cout<<std::setprecision(12);
 std::cout<<"cells,steps,mean_energy_keV,electron_fraction,ion_fraction,particle_residual,energy_residual,maxwellian_L1,max_step_energy_residual\n";
 std::array<double,3> fractions{};
 int k=0;
 for(int n:{200,400,800}) {
   double previous_error=1.;
   for(int steps:{250,500,1000}) {
     const auto r=run(n,steps);
     if(r.equilibrium_error>=previous_error)
       throw std::runtime_error("time refinement failed to reduce relaxation error");
     previous_error=r.equilibrium_error;
     fractions[k]=r.electron_fraction;
   }
   ++k;
 }
 const double coarse=std::abs(fractions[1]-fractions[0]);
 const double fine=std::abs(fractions[2]-fractions[1]);
 if(fine>0.4*coarse || fine>1e-4)
   throw std::runtime_error("electron/ion partition failed grid convergence");
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
 return 0;
}
