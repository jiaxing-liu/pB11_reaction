#include "fusion_flux_inverse.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void ck(bool x){if(!x)throw std::runtime_error("inverse check failed");}
int main(){try{
 double x[]={0,.2,.5,1},a[12],d[12],F[]={1.6,1.6,1.6,1.6},C[]={0,.016,.04,.08};
 for(int i=0;i<4;++i){a[3*i]=2+.01*x[i]*x[i];a[3*i+1]=.4*x[i];a[3*i+2]=.6*x[i];d[3*i]=.02*x[i];d[3*i+1]=.4;d[3*i+2]=.6;}
 fusion_flux_inverse_value_v1 out{};fusion_flux_field_value_v1 p{};
 for(int k:{1,2})for(double q:{0.,1e-12,1e-10,1e-8,.01,.3,.8,1.-1e-10,1.})for(double t:{-3.,-.3,0.,.7,2.,3.141592653589793})for(int mode:{0,1}){
  ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,k,q,t,&p)==0);
  ck(fusion_c_flux_inverse(4,1,x,a,d,F,C,k,p.R_m,p.Z_m,2e-14,60,mode,0,0,&out)==0);
  ck(std::abs(out.u-q*std::cos(t))<2e-13&&std::abs(out.v-q*std::sin(t))<2e-13);
  ck(out.residual_m<=2e-14);
 }
 ck(fusion_c_flux_inverse(4,1,x,a,d,F,C,1,5,0,1e-13,40,0,0,0,&out)==FUSION_FIELD_NOT_CONVERGED);ck(out.xi==0&&out.residual_m==0);
 ck(fusion_c_flux_core_eval(4,1,x,a,d,F,C,1,.6,.7,&p)==0);
 ck(fusion_c_flux_inverse(4,1,x,a,d,F,C,1,p.R_m,p.Z_m,1e-14,1,1,0,0,&out)==FUSION_FIELD_NOT_CONVERGED);
 ck(fusion_c_flux_inverse(4,1,x,a,d,F,C,1,2,0,0,40,0,0,0,&out)==PB11_STATUS_INVALID_ARGUMENT);
 ck(fusion_c_flux_inverse(4,1,x,a,d,F,C,1,2,0,1e-14,40,0,0,0,nullptr)==PB11_STATUS_NULL_OUTPUT);
 d[1]=0;
 ck(fusion_c_flux_inverse(4,1,x,a,d,F,C,1,2.1,0,1e-14,40,0,0,0,&out)==PB11_STATUS_NUMERICAL_FAILURE);
 ck(out.u==0&&out.v==0&&out.xi==0&&out.theta_rad==0&&out.residual_m==0&&out.iterations==0);
 std::cout<<"inverse axis/core/boundary roundtrips, seed agreement and failure separation passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
