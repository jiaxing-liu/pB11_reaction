#include "fusion_flux_inverse.h"
#include <algorithm>
#include <cmath>
#include <limits>
extern "C" int fusion_c_flux_inverse(int n,int m,const double*x,
 const double*a,const double*d,const double*F,const double*C,int k,
 double R,double Z,double tol,int limit,int seed_mode,double seed_u,double seed_v,
 fusion_flux_inverse_value_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(!std::isfinite(R)||!std::isfinite(Z)||!std::isfinite(tol)||tol<=0||limit<1||limit>100||(seed_mode!=0&&seed_mode!=1)||!std::isfinite(seed_u)||!std::isfinite(seed_v))return PB11_STATUS_INVALID_ARGUMENT;
 fusion_flux_field_value_v1 axis{};
 int status=fusion_c_flux_core_eval(n,m,x,a,d,F,C,k,0,0,&axis);
 if(status)return status;
 if(R==axis.R_m&&Z==0){out->iterations=0;return PB11_STATUS_OK;}
 using L=long double;
 L u=seed_mode?seed_u:(L(R)-a[0])/d[1],v=seed_mode?seed_v:L(Z)/d[m+1];
 auto project=[&](L& U,L& V){L r=std::hypot(U,V);if(r>x[n-1]){U*=L(x[n-1])/r;V*=L(x[n-1])/r;}};
 project(u,v);
 auto eval=[&](L U,L V,fusion_flux_field_value_v1&value){
  double q=static_cast<double>(std::hypot(U,V));
  // Only conversion-level rounding: the long-double trial is projected first.
  q=std::min(q,x[n-1]);
  return fusion_c_flux_core_eval(n,m,x,a,d,F,C,k,q,std::atan2(V,U),&value);
 };
 for(int iter=0;iter<limit;++iter){
  fusion_flux_field_value_v1 b{};status=eval(u,v,b);if(status)return status;
  L rR=L(R)-b.R_m,rZ=L(Z)-b.Z_m,res=std::hypot(rR,rZ),q=std::hypot(u,v);
  if(res<=tol){*out={double(u),double(v),double(q),q==0?0:double(std::atan2(v,u)),double(res),iter};return PB11_STATUS_OK;}
  L Ru,Rv,Zu,Zv;
  if(q==0){Ru=d[1];Rv=0;Zu=0;Zv=d[m+1];}
  else{L ct=u/q,st=v/q;Ru=b.R_xi_m*ct-b.R_theta_m*st/q;Rv=b.R_xi_m*st+b.R_theta_m*ct/q;Zu=b.Z_xi_m*ct-b.Z_theta_m*st/q;Zv=b.Z_xi_m*st+b.Z_theta_m*ct/q;}
  L det=Ru*Zv-Rv*Zu,scale=std::hypot(Ru,Rv)*std::hypot(Zu,Zv);
  if(!std::isfinite(det)||scale==0||std::abs(det)<=64*std::numeric_limits<double>::epsilon()*scale)return FUSION_FIELD_NOT_CONVERGED;
  L du=(rR*Zv-rZ*Rv)/det,dv=(Ru*rZ-Zu*rR)/det;
  bool accepted=false;
  for(int back=0;back<24;++back){
   L alpha=std::ldexp(1.L,-back),nu=u+alpha*du,nv=v+alpha*dv;
   if(!std::isfinite(nu)||!std::isfinite(nv))continue;
   project(nu,nv);fusion_flux_field_value_v1 trial{};
   if(eval(nu,nv,trial)!=PB11_STATUS_OK)continue;
   L nr=std::hypot(L(R)-trial.R_m,L(Z)-trial.Z_m);
   if(nr<res){u=nu;v=nv;accepted=true;break;}
  }
  if(!accepted)return FUSION_FIELD_NOT_CONVERGED;
 }
 return FUSION_FIELD_NOT_CONVERGED;
}
