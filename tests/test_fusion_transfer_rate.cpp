#include "fusion_two_component.h"
#include <cmath>
#include <cstdio>
#include <limits>
#include <initializer_list>
#include <stdexcept>
void require(bool b){if(!b)throw std::runtime_error("transfer coefficient validation");}
int main(){
 const double pi=acos(-1.),ma=6.6446573357e-27,keV=1.602176634e-16;
 double worst=0;
 for(double mb:{3.3435837768e-27,5.0073567446e-27}) {
  fusion_maxwellian_bath_v1 b{5e19,mb,1,10*keV,15};double vt=sqrt(2*b.kT_J/mb);
  auto radial=[&](double y){double v=y*vt;fusion_coulomb_energy_v1 o{};require(!fusion_c_coulomb_energy(.5*ma*v*v,ma,2,&b,&o));return v*o.diffusion_J2_s/(ma*b.kT_J);};
  double r0;require(!fusion_c_coulomb_transfer_rate(0,ma,2,&b,&r0));require(r0>0);
  for(double y:{.001,.01,.1,.5,1.,2.,3.}) {
   double h=1e-3*y,r;require(!fusion_c_coulomb_transfer_rate(.5*ma*y*y*vt*vt,ma,2,&b,&r));
   double num=(-radial(y+2*h)+8*radial(y+h)-8*radial(y-h)+radial(y-2*h))/(12*h*vt*pow(y*vt,2));
   worst=fmax(worst,fabs(num/r-1));require(fabs(num/r-1)<1e-7);
   require(fabs(r/r0-exp(-y*y))<1e-14);
  }
  // Integrate returned rate over velocity space, compared to independently
  // measured high-speed friction flux: int lambda d3v = 4*pi*nu.
  long double integral=0;int n=10000;for(int j=0;j<=n;j++){
   long double y=10.L*j/n;double r;require(!fusion_c_coulomb_transfer_rate(.5*ma*y*y*vt*vt,ma,2,&b,&r));
   integral+=(j==0||j==n?1:j%2?4:2)*r*y*y;
  }
  integral*=4*pi*pow((long double)vt,3)*10.L/n/3;
  require(fabsl(integral/(4*pi*radial(10))-1)<1e-12);
  double q;require(!fusion_c_coulomb_transfer_rate(1e6*keV,ma,2,&b,&q));require(q==0);
  require(!fusion_c_coulomb_transfer_rate(0,ma,0,&b,&q));require(q==0);
  require(fusion_c_coulomb_transfer_rate(-1,ma,2,&b,&q)==PB11_STATUS_OUT_OF_RANGE&&q==0);
  require(fusion_c_coulomb_transfer_rate(0,ma,2,nullptr,&q)==PB11_STATUS_INVALID_ARGUMENT&&q==0);
  require(fusion_c_coulomb_transfer_rate(0,ma,2,&b,nullptr)==PB11_STATUS_NULL_OUTPUT);
  require(fusion_c_coulomb_transfer_rate(NAN,ma,2,&b,&q)==PB11_STATUS_INVALID_ARGUMENT&&q==0);
  auto bad=b;bad.kT_J=-1;require(fusion_c_coulomb_transfer_rate(0,ma,2,&bad,&q)==PB11_STATUS_OUT_OF_RANGE&&q==0);
  bad=b;bad.density_m3=std::numeric_limits<double>::max();require(fusion_c_coulomb_transfer_rate(0,ma,1e100,&bad,&q)==PB11_STATUS_NUMERICAL_FAILURE&&q==0);
 }
 printf("max friction-divergence relative error %.12g; radial integral and ABI/domain checks pass\n",worst);
}
