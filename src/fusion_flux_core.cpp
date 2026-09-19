#include "fusion_flux_core.h"
#include <cmath>
#include <limits>
namespace {
using L=long double;
struct Pair {L y,d;};
Pair even(L a,L v,L slope,L xc,L t){
 L delta=v-a,s=xc*slope,alpha=2*delta-s/2,beta=s/2-delta,t2=t*t;
 return {a+t2*(alpha+beta*t2),(2*alpha*t+4*beta*t*t2)/xc};
}
Pair mode(int m,L v,L slope,L xc,L t){
 L beta=(xc*slope-m*v)/2,alpha=v-beta;
 return {std::pow(t,m)*(alpha+beta*t*t),std::pow(t,m-1)*(m*alpha+(m+2)*beta*t*t)/xc};
}
Pair first(L axis_slope,L v,L slope,L xc,L t){
 L G=v/xc,b=(5*G-slope-4*axis_slope)/2,d=(slope-3*G+2*axis_slope)/2,t2=t*t;
 return {xc*t*(axis_slope+b*t2+d*t2*t2),axis_slope+3*b*t2+5*d*t2*t2};
}
}
extern "C" int fusion_c_flux_core_eval(int n,int m,const double*x,
 const double*a,const double*d,const double*F,const double*C,int k,double q,double theta,
 fusion_flux_field_value_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(n<3||n>10000||m<1||m>64||k<1||k>=n-1)return PB11_STATUS_OUT_OF_RANGE;
 if(!x||!a||!d||!F||!C||!std::isfinite(q)||!std::isfinite(theta))return PB11_STATUS_INVALID_ARGUMENT;
 fusion_flux_field_value_v1 boundary{};
 int status=fusion_c_flux_field_eval(n,m,x,a,d,F,C,x[k],theta,&boundary);
 if(status)return status;
 if(q<0||q>x[n-1])return PB11_STATUS_OUT_OF_RANGE;
 int stride=1+2*m;
 if(!(a[0]>0)||C[0]!=0)return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=1;j<stride;++j)if(a[j]!=0)return PB11_STATUS_INVALID_ARGUMENT;
 L ar=d[1],az=d[m+1],largest=std::fmax(std::abs(ar),std::abs(az));
 if(largest==0||std::fmin(std::abs(ar),std::abs(az))/largest<=64*std::numeric_limits<double>::epsilon())return PB11_STATUS_NUMERICAL_FAILURE;
 if(q>=x[k])return fusion_c_flux_field_eval(n,m,x,a,d,F,C,q,theta,out);
 L xc=x[k],t=L(q)/xc,angle=std::remainder(L(theta),2*std::acos(-1.L));
 auto geom=[&](int j,int harmonic){
  if(harmonic==0)return even(a[0],a[k*stride],d[k*stride],xc,t);
  if(harmonic==1)return first(d[j],a[k*stride+j],d[k*stride+j],xc,t);
  return mode(harmonic,a[k*stride+j],d[k*stride+j],xc,t);
 };
 auto base=geom(0,0);L R=base.y,Rx=base.d,Z=0,Zx=0,Rt=0,Zt=0;
 for(int j=1;j<=m;++j){auto r=geom(j,j),z=geom(m+j,j);L s=std::sin(j*angle),c=std::cos(j*angle);
  R+=r.y*c;Rx+=r.d*c;Z+=z.y*s;Zx+=z.d*s;Rt-=j*r.y*s;Zt+=j*z.y*c;
 }
 L h=L(x[k+1])-x[k];
 L f=even(F[0],F[k],(L(F[k+1])-F[k])/h,xc,t).y;
 L c=mode(1,C[k],(L(C[k+1])-C[k])/h,xc,t).y;
 L det=Zt*Rx-Rt*Zx,J=R*det,br=0,bz=0;
 if(!(R>0))return PB11_STATUS_NUMERICAL_FAILURE;
 if(q!=0){
  L scale=std::hypot(Rt,Zt)*std::hypot(Rx,Zx);
  if(!std::isfinite(J)||!std::isfinite(scale)||scale==0||std::abs(det)<=64*std::numeric_limits<double>::epsilon()*scale||((det>0)!=(ar*az>0))||static_cast<double>(J)==0)return PB11_STATUS_NUMERICAL_FAILURE;
  br=c*Rt/J;bz=c*Zt/J;
 }
 L values[]={R,Z,Rx,Zx,Rt,Zt,J,br,bz,f/R};
 for(L v:values)if(!std::isfinite(v)||std::abs(v)>std::numeric_limits<double>::max())return PB11_STATUS_NUMERICAL_FAILURE;
 *out={double(R),double(Z),double(Rx),double(Zx),double(Rt),double(Zt),double(J),double(br),double(bz),double(f/R)};
 return PB11_STATUS_OK;
}
