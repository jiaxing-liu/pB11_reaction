#include "fusion_flux_volume_diagnostic.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

extern "C" int fusion_c_flux_core_enclosed_volume_diagnostic_v1(
 int n,int m,const double*x,const double*a,const double*d,const double*F,
 const double*C,int core,double q,fusion_flux_volume_diagnostic_v1*out) {
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 try {
  if(n<3||n>10000||m<1||m>64||core<1||core>=n-1)
   return PB11_STATUS_OUT_OF_RANGE;
  fusion_flux_field_value_v1 value{};
  int rc=fusion_c_flux_core_eval(n,m,x,a,d,F,C,core,q,0,&value);
  if(rc)return rc;
  if(q==0){out->version=1;out->angular_nodes=1;return PB11_STATUS_OK;}
  int count=1;
  while(count<=3*m)count*=2;
  std::array<double,256> rtheta{},ztheta{};
  std::array<long double,64> radial_coeff{};
  const long double pi=std::acos(-1.L);
  long double mean_radius=0;
  int orientation=0;
  double min_jacobian=std::numeric_limits<double>::max();
  for(int j=0;j<count;++j){
   if(j){
    rc=fusion_c_flux_core_eval(n,m,x,a,d,F,C,core,q,
                             double(2*pi*j/count),&value);
    if(rc)return rc;
   }
   const double J=value.signed_jacobian_m3;
   const int sign=J>0?1:(J<0?-1:0);
   if(!sign||(orientation&&sign!=orientation))return PB11_STATUS_NUMERICAL_FAILURE;
   orientation=sign;
   min_jacobian=std::min(min_jacobian,std::abs(J));
   rtheta[j]=value.R_theta_m;ztheta[j]=value.Z_theta_m;
   mean_radius+=static_cast<long double>(value.R_m)/count;
  }
  // Recover the zero-mean radius from its angular derivative. Subtracting
  // nearly equal absolute radii would lose relative accuracy near the axis.
  for(int h=1;h<=m;++h){
   long double projection=0;
   for(int j=0;j<count;++j)
    projection+=rtheta[j]*std::sin(h*2*pi*j/count);
   radial_coeff[h-1]=-2*projection/(count*h);
  }
  long double sum=0;
  for(int j=0;j<count;++j){
   long double delta=0;
   for(int h=1;h<=m;++h)delta+=radial_coeff[h-1]*std::cos(h*2*pi*j/count);
   sum+=(2*mean_radius*delta+delta*delta)*ztheta[j];
  }
  const long double volume=2*pi*pi*sum/count;
  if(!std::isfinite(volume)||std::abs(volume)>std::numeric_limits<double>::max()||
     volume==0||(volume>0?1:-1)!=orientation||static_cast<double>(volume)==0)
   return PB11_STATUS_NUMERICAL_FAILURE;
  out->version=1;out->angular_nodes=count;out->boundary_orientation=orientation;
  out->signed_volume_m3=static_cast<double>(volume);
  out->absolute_volume_m3=std::abs(out->signed_volume_m3);
  out->min_abs_boundary_jacobian_m3=min_jacobian;
  return PB11_STATUS_OK;
 }catch(...){*out={};return PB11_STATUS_EXCEPTION;}
}
