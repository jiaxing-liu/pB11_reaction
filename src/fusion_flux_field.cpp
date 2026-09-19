#include "fusion_flux_field.h"
#include <algorithm>
#include <cmath>
#include <limits>
namespace {
using L=long double;
bool finite_value(double x){return std::isfinite(x);}
struct Pair{L y,d;};
Pair hermite(L y0,L y1,L d0,L d1,L h,L t){
 if(t==0)return {y0,d0};
 if(t==1)return {y1,d1};
 L t2=t*t,t3=t2*t;
 return {(2*t3-3*t2+1)*y0+(t3-2*t2+t)*h*d0+(-2*t3+3*t2)*y1+(t3-t2)*h*d1,
 ((6*t2-6*t)*y0+(-6*t2+6*t)*y1)/h+(3*t2-4*t+1)*d0+(3*t2-2*t)*d1};
}
}
extern "C" int fusion_c_flux_field_eval(int n,int m,const double*x,
 const double*a,const double*d,const double*F,const double*C,double q,double theta,
 fusion_flux_field_value_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(n<2||n>10000||m<1||m>64)return PB11_STATUS_OUT_OF_RANGE;
 if(!x||!a||!d||!F||!C||!finite_value(q)||!finite_value(theta))return PB11_STATUS_INVALID_ARGUMENT;
 const int stride=1+2*m;
 for(int i=0;i<n;++i){
  if(!finite_value(x[i])||!finite_value(F[i])||!finite_value(C[i])||(i&&x[i]<=x[i-1]))return PB11_STATUS_INVALID_ARGUMENT;
  for(int k=0;k<stride;++k)if(!finite_value(a[i*stride+k])||!finite_value(d[i*stride+k]))return PB11_STATUS_INVALID_ARGUMENT;
 }
 if(x[0]!=0)return PB11_STATUS_INVALID_ARGUMENT;
 if(q<0||q>x[n-1])return PB11_STATUS_OUT_OF_RANGE;
 if(q==0)return FUSION_FIELD_AXIS_COORDINATE_SINGULAR;
 int i=std::min(int(std::upper_bound(x,x+n,q)-x)-1,n-2);
 L h=L(x[i+1])-x[i],t=(L(q)-x[i])/h;
 auto p=[&](int k){return hermite(a[i*stride+k],a[(i+1)*stride+k],d[i*stride+k],d[(i+1)*stride+k],h,t);};
 auto base=p(0);L R=base.y,Rx=base.d,Z=0,Zx=0,Rt=0,Zt=0;
 L angle=std::remainder(L(theta),2*std::acos(-1.L));
 for(int k=1;k<=m;++k){auto r=p(k),z=p(m+k);L s=std::sin(k*angle),c=std::cos(k*angle);
  R+=r.y*c;Rx+=r.d*c;Z+=z.y*s;Zx+=z.d*s;Rt-=k*r.y*s;Zt+=k*z.y*c;
 }
 L det=Zt*Rx-Rt*Zx,scale=std::hypot(Rt,Zt)*std::hypot(Rx,Zx),J=R*det;
 if(!(R>0)||!std::isfinite(J)||!std::isfinite(scale)||scale==0||std::abs(det)<=64*std::numeric_limits<double>::epsilon()*scale)return PB11_STATUS_NUMERICAL_FAILURE;
 L f=(1-t)*F[i]+t*F[i+1],c=(1-t)*C[i]+t*C[i+1];
 if(static_cast<double>(J)==0)return PB11_STATUS_NUMERICAL_FAILURE;
 L values[]={R,Z,Rx,Zx,Rt,Zt,J,c*Rt/J,c*Zt/J,f/R};
 for(L v:values)if(!std::isfinite(v)||std::abs(v)>std::numeric_limits<double>::max())return PB11_STATUS_NUMERICAL_FAILURE;
 *out={double(R),double(Z),double(Rx),double(Zx),double(Rt),double(Zt),double(J),double(values[7]),double(values[8]),double(values[9])};
 return PB11_STATUS_OK;
}
