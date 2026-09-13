// Root research prototype: fixed equal-temperature baths, explicit Maxwellian
// projection tolerance. Not yet a general evolving-background ash interface.
#include "fusion_two_component.h"
#include <boost/math/special_functions/gamma.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
using R=long double;
constexpr double kev=1.602176634e-16, temp=10*kev, ma=6.644657345e-27,N0=1e12,E0=3000*kev;
void ok(bool b){if(!b)throw std::runtime_error("handoff prototype check");}
R sum(const std::vector<double>& x){R v=0;for(double q:x)v+=q;return v;}
int main(int argc,char** argv){try{
 int n=argc>1?atoi(argv[1]):2000,steps=argc>2?atoi(argv[2]):2000;double tol=argc>3?atof(argv[3]):.001;ok(n>=100&&steps>0&&tol>0);
 double dt=1./steps;
 fusion_maxwellian_bath_v1 baths[3]={{1e20,9.1093837139e-31,1,temp,15},{5e19,3.3435837768e-27,1,temp,15},{5e19,5.0073567512e-27,1,temp,15}};
 double ts[3]={temp,temp,temp};std::vector<double> edges(n+1),c(n),s(n),t(n),f(n),sn(n),tn(n),fn(n),z(n),birth(n),lam(n),D(3*(n-1)),q(n);
 int low=n/4;edges[0]=0;for(int i=1;i<=low;i++)edges[i]=.001*kev*pow(20000.,double(i-1)/(low-1));for(int i=low+1;i<=n;i++)edges[i]=(20.+5980.*(i-low)/(n-low))*kev;
 R qsum=0,qU=0;for(int i=0;i<n;i++){
  c[i]=(edges[i]+edges[i+1])/2;R a=edges[i]/temp,b=edges[i+1]/temp;
  R qi=a>=1.5L?boost::math::gamma_q(1.5L,a)-boost::math::gamma_q(1.5L,b):boost::math::gamma_p(1.5L,b)-boost::math::gamma_p(1.5L,a);ok(qi>=0);q[i]=qi;qsum+=qi;qU+=qi*c[i];
  for(int j=1;j<3;j++){double rate;ok(!fusion_c_coulomb_transfer_rate(c[i],ma,2,&baths[j],&rate));lam[i]+=rate;}
 }
 for(int j=0;j<3;j++)for(int i=0;i<n-1;i++){fusion_coulomb_energy_v1 out{};ok(!fusion_c_coulomb_energy(edges[i+1],ma,2,&baths[j],&out));D[j*(n-1)+i]=out.diffusion_J2_s;}
 auto r=std::upper_bound(c.begin(),c.end(),E0)-c.begin();double w=(E0-c[r-1])/(c[r]-c[r-1]);s[r]=w*N0;s[r-1]=(1-w)*N0;f=s;
 R fluid=0,projectionU=0,shapeBudget=0,maxL1=0,maxE=0,maxN=0;std::array<R,3> heat{},fullheat{};int events=0;double first=-1;
 printf("time_s,fluid_N_fraction,ST_N_fraction,TH_N_fraction,half_L1_vs_full,projection_shape_bound,electron_heat_fraction,ion_heat_fraction,fluid_U_fraction,kinetic_U_fraction,number_residual,energy_residual,projection_heat_correction_fraction,fluid_grid_energy_bias_fraction,events\n");
 for(int step=1;step<=steps;step++){
  double hs[3],ht[3],hf[3];fusion_kinetic_ledger_v1 ls{},lt{},lf{};
  ok(!fusion_c_energy_fp_trial(n,3,dt,edges.data(),s.data(),ts,D.data(),z.data(),lam.data(),0,sn.data(),hs,&ls));
  for(int i=0;i<n;i++)birth[i]=lam[i]*sn[i];
  ok(!fusion_c_energy_fp_trial(n,3,dt,edges.data(),t.data(),ts,D.data(),birth.data(),z.data(),0,tn.data(),ht,&lt));
  ok(!fusion_c_energy_fp_trial(n,3,dt,edges.data(),f.data(),ts,D.data(),z.data(),z.data(),0,fn.data(),hf,&lf));
  for(int j=0;j<3;j++){heat[j]+=hs[j]+ht[j];fullheat[j]+=hf[j];}
  R Nt=sum(tn),Ut=0,shape=1-qsum;for(int i=0;i<n;i++){Ut+=R(tn[i])*c[i];if(Nt>0)shape+=fabsl(R(tn[i])/Nt-q[i]);}
  // Full bin-probability L1 includes the exact Maxwellian tail outside grid.
  // A separate energy bound prevents number-norm proximity hiding a hot tail.
  R energyError=Nt>0?fabsl(Ut/(1.5L*temp*Nt)-1):0;
  if(Nt>0 && shape<=tol && energyError<=tol){
   fluid+=Nt;R dU=Ut-1.5L*temp*Nt;projectionU+=dU;heat[1]+=dU*.5L;heat[2]+=dU*.5L;
   // D/T have common Ti, so a declared equal-number heat capacity split.
   shapeBudget+=.5L*Nt*shape;events++;if(first<0)first=step*dt;std::fill(tn.begin(),tn.end(),0.);
  }
  s.swap(sn);t.swap(tn);f.swap(fn);R kinU=0,l1=0;
  for(int i=0;i<n;i++){kinU+=(R(s[i])+t[i])*c[i];l1+=fabsl(R(s[i])+t[i]+fluid*q[i]-f[i]);}
  l1=.5L*(l1+fluid*(1-qsum))/N0;maxL1=std::max(maxL1,l1);
  R nr=(sum(s)+sum(t)+fluid)/N0-1,er=(kinU+1.5L*temp*fluid+heat[0]+heat[1]+heat[2])/(N0*E0)-1;maxN=std::max(maxN,fabsl(nr));maxE=std::max(maxE,fabsl(er));
  if(step%std::max(1,steps/200)==0)printf("%.9g,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%.14Lg,%d\n",step*dt,fluid/N0,sum(s)/N0,sum(t)/N0,l1,shapeBudget/N0,heat[0]/(N0*E0),(heat[1]+heat[2])/(N0*E0),1.5L*temp*fluid/(N0*E0),kinU/(N0*E0),nr,er,projectionU/(N0*E0),fluid*(qU-1.5L*temp)/(N0*E0),events);
 }
 fprintf(stderr,"n=%d steps=%d tol=%.6g first_projection_s=%.9g events=%d final_fluid_N=%.12Lg max_half_L1=%.12Lg max_N=%.12Lg max_E=%.12Lg exact_bin_grid_mean_error=%.12Lg heat_delta_e=%.12Lg heat_delta_ions=%.12Lg\n",n,steps,tol,first,events,fluid/N0,maxL1,maxN,maxE,qU/(1.5L*temp)-1,(heat[0]-fullheat[0])/(N0*E0),(heat[1]+heat[2]-fullheat[1]-fullheat[2])/(N0*E0));
 ok(maxN<1e-10&&maxE<1e-10);
 }catch(const std::exception&e){fprintf(stderr,"ERROR %s\n",e.what());return 1;}}
