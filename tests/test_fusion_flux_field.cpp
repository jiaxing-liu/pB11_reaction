#include "fusion_flux_field.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
void check(bool x){if(!x)throw std::runtime_error("flux field check failed");}
void near(double a,double b,double tol=3e-13){check(std::abs(a-b)<=tol*std::max(1.,std::abs(b)));}
int main(){try{
 const double pi=std::acos(-1.);double x[]={0,.5,1},a[]={2,0,0,2,.2,.2,2,.4,.4},d[]={0,.4,.4,0,.4,.4,0,.4,.4};
 double F[]={1.6,1.6,1.6},C[]={0,.04,.08};fusion_flux_field_value_v1 v{};
 auto eval=[&](double q,double t){check(fusion_c_flux_field_eval(3,1,x,a,d,F,C,q,t,&v)==0);return v;};
 for(double q:{1e-8,.1,.25,.5,.75,1.})for(double t:{-4.,0.,.3,pi,2*pi}){
  auto b=eval(q,t);double r=.4*q,R=2+r*std::cos(t);
  near(b.R_m,R);near(b.Z_m,r*std::sin(t));near(b.R_xi_m,.4*std::cos(t));near(b.Z_xi_m,.4*std::sin(t));
  near(b.signed_jacobian_m3,R*.16*q);near(b.BR_T,-.2*q*std::sin(t)/R);near(b.BZ_T,.2*q*std::cos(t)/R);near(b.Bphi_T,1.6/R);
  near(b.BR_T*b.Z_theta_m-b.BZ_T*b.R_theta_m,0);
  auto p=eval(q,t+2*pi);near(p.BR_T,b.BR_T);near(p.BZ_T,b.BZ_T);
 }
 // Independent cylindrical finite-difference divergence using analytic inverse circle.
 auto field=[&](double R,double Z){return eval(std::hypot(R-2,Z)/.4,std::atan2(Z,R-2));};
 for(double h:{1e-3,5e-4,2.5e-4}){
  double R=2.13,Z=.07;
  double div=((R+h)*field(R+h,Z).BR_T-(R-h)*field(R-h,Z).BR_T)/(2*h*R)+(field(R,Z+h).BZ_T-field(R,Z-h).BZ_T)/(2*h);
  check(std::abs(div)<1e-11);
 }
 // Nonlinear radial radius: cubic Hermite must reproduce coefficient AND slope.
 for(int i=0;i<3;++i){double q=x[i];a[3*i+1]=a[3*i+2]=.3*q+.1*q*q*q;d[3*i+1]=d[3*i+2]=.3+.3*q*q;}
 for(double q:{.125,.375,.625,.875}){auto b=eval(q,.7);near(b.R_m,2+(.3*q+.1*q*q*q)*std::cos(.7));near(b.R_xi_m,(.3+.3*q*q)*std::cos(.7));}
 // Reflection of geometry reverses J; signed field follows the mapping.
 auto before=eval(.4,.7);for(int i=0;i<3;++i){a[3*i+2]*=-1;d[3*i+2]*=-1;}auto after=eval(.4,.7);
 near(after.signed_jacobian_m3,-before.signed_jacobian_m3);near(after.BR_T,-before.BR_T);near(after.BZ_T,before.BZ_T);
 auto bad=[&](double q,double t,int expected){v.R_m=123;check(fusion_c_flux_field_eval(3,1,x,a,d,F,C,q,t,&v)==expected);check(v.R_m==0&&v.BR_T==0&&v.Bphi_T==0);};
 bad(0,0,FUSION_FIELD_AXIS_COORDINATE_SINGULAR);bad(-.1,0,PB11_STATUS_OUT_OF_RANGE);bad(1.1,0,PB11_STATUS_OUT_OF_RANGE);bad(.2,std::numeric_limits<double>::quiet_NaN(),PB11_STATUS_INVALID_ARGUMENT);
 x[1]=0;bad(.2,0,PB11_STATUS_INVALID_ARGUMENT);x[1]=.5;
 d[0]=std::numeric_limits<double>::infinity();bad(.2,0,PB11_STATUS_INVALID_ARGUMENT);d[0]=0;
 for(int i=0;i<3;++i){a[3*i+2]=0;d[3*i+2]=0;}bad(.2,0,PB11_STATUS_NUMERICAL_FAILURE);
 check(fusion_c_flux_field_eval(3,1,x,a,d,F,C,.2,0,nullptr)==PB11_STATUS_NULL_OUTPUT);
 std::cout<<"flux geometry, signed field, cubic slopes, periodicity, divergence and failure contracts passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
