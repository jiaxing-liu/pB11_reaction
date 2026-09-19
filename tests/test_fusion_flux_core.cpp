#include "fusion_flux_core.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <algorithm>
void ck(bool x){if(!x)throw std::runtime_error("regular core check failed");}
void near(double a,double b,double tol=3e-13){ck(std::abs(a-b)<=tol*std::max(1.,std::abs(b)));}
int main(){try{
 double x[]={0,.2,.5,1},a[12],d[12],F[]={1.6,1.6,1.6,1.6},C[]={0,.016,.04,.08};
 for(int i=0;i<4;++i){a[3*i]=2;a[3*i+1]=.4*x[i];a[3*i+2]=.6*x[i];d[3*i]=0;d[3*i+1]=.4;d[3*i+2]=.6;}
 fusion_flux_field_value_v1 v{},old{};
 for(int k:{1,2})for(double q:{0.,1e-10,1e-8,.01,.15,.2,.3,.5,.8,1.})for(double t:{-4.,0.,.7,3.14,7.}){
  ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,k,q,t,&v)==0);
  double R=2+.4*q*std::cos(t);
  near(v.R_m,R);near(v.Z_m,.6*q*std::sin(t));near(v.BR_T,-.08*q*std::sin(t)/(R*.6));near(v.BZ_T,.08*q*std::cos(t)/(R*.4));near(v.Bphi_T,1.6/R);
  if(q==0){ck(v.BR_T==0&&v.BZ_T==0&&v.signed_jacobian_m3==0);near(v.Bphi_T,.8);}
  if(q>=x[k]){ck(fusion_c_flux_field_eval(4,1,x,a,d,F,C,q,t,&old)==0);ck(std::memcmp(&v,&old,sizeof v)==0);}
 }
 // Three-mode regular polynomial geometry: exact within each candidate patch.
 double aa[28]{},dd[28]{};
 for(int i=0;i<4;++i){double q=x[i];aa[7*i]=2+.01*q*q+.02*std::pow(q,4);dd[7*i]=.02*q+.08*std::pow(q,3);
  for(int m=1;m<=3;++m){aa[7*i+m]=.02*std::pow(q,m)+.01*std::pow(q,m+2);dd[7*i+m]=.02*m*std::pow(q,m-1)+.01*(m+2)*std::pow(q,m+1);aa[7*i+3+m]=aa[7*i+m]*1.3;dd[7*i+3+m]=dd[7*i+m]*1.3;}
 }
 for(int k:{1,2})for(double q:{0.,1e-8,.01,.07,.15}){
  double t=.8;ck(fusion_c_flux_core_eval(4,3,x,aa,dd,F,C,k,q,t,&v)==0);
  double R=2+.01*q*q+.02*std::pow(q,4),Z=0,Rx=.02*q+.08*std::pow(q,3);
  for(int m=1;m<=3;++m){R+=(.02*std::pow(q,m)+.01*std::pow(q,m+2))*std::cos(m*t);Z+=1.3*(.02*std::pow(q,m)+.01*std::pow(q,m+2))*std::sin(m*t);Rx+=(.02*m*std::pow(q,m-1)+.01*(m+2)*std::pow(q,m+1))*std::cos(m*t);}
  near(v.R_m,R);near(v.Z_m,Z);near(v.R_xi_m,Rx);
 }
 // Local continuity at both matching surfaces, including geometric derivatives.
 for(int k:{1,2}){double q=x[k];ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,k,std::nextafter(q,0.),.7,&v)==0);ck(fusion_c_flux_field_eval(4,1,x,a,d,F,C,q,.7,&old)==0);near(v.R_xi_m,old.R_xi_m);near(v.Z_xi_m,old.Z_xi_m);near(v.BR_T,old.BR_T);near(v.BZ_T,old.BZ_T);near(v.Bphi_T,old.Bphi_T);}
 ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,0,0,0,&v)==PB11_STATUS_OUT_OF_RANGE);ck(v.R_m==0);
 C[0]=1;ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,1,0,0,&v)==PB11_STATUS_INVALID_ARGUMENT);C[0]=0;
 a[1]=.01;ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,1,0,0,&v)==PB11_STATUS_INVALID_ARGUMENT);a[1]=0;
 d[1]=0;ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,1,0,0,&v)==PB11_STATUS_NUMERICAL_FAILURE);ck(v.Bphi_T==0);
 std::cout<<"regular-core analytic geometry, axis, matching, outer parity and failure checks passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
