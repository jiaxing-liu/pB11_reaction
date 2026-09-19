#include "fusion_magnetic_push.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
void ck(bool b){if(!b)throw std::runtime_error("magnetic push check failed");}
struct Context {double B[3];int calls=0;double at[3]={};int status=0;};
int field(void*p,const double*x,double*B){auto&c=*static_cast<Context*>(p);++c.calls;for(int j=0;j<3;++j){c.at[j]=x[j];B[j]=c.B[j];}return c.status;}
double norm(const double*v){return std::hypot(std::hypot(v[0],v[1]),v[2]);}
int main(){try{
 constexpr double c=299792458.;double mid[3];
 for(double U:{1e5,2e8})for(double sign:{-1.,1.}){
  const double gamma=std::sqrt(1+U*U/(c*c)),T=1.3*gamma;
  double previous=0;
  for(int N:{20,40,80}){
   fusion_orbit_state_v1 s{{0,0,0},{U,0,0}},next{};Context f{{0,0,1}};
   for(int i=0;i<N;++i){ck(fusion_c_magnetic_push(1,sign,T/N,&s,field,&f,&next,mid)==0);s=next;}
   const double exactx=U*std::sin(1.3),exacty=-sign*U*(1-std::cos(1.3));
   double error=std::hypot(s.x_m[0]-exactx,s.x_m[1]-exacty)/U;
   ck(std::abs(norm(s.u_m_s)/U-1)<3e-14);
   ck(std::abs(s.u_m_s[0]/U-std::cos(1.3))<3e-14&&std::abs(s.u_m_s[1]/U+sign*std::sin(1.3))<3e-14);
   if(previous)ck(previous/error>3.99&&previous/error<4.01);
   std::cout<<"U="<<U<<" sign="<<sign<<" N="<<N<<" position_error/U="<<error<<'\n';previous=error;
  }
 }
 Context f{{.2,-.1,.4}};fusion_orbit_state_v1 s{{.1,.2,.3},{1e7,2e7,-3e7}},b{},back{};
 ck(fusion_c_magnetic_push(1,1,.001,&s,field,&f,&b,mid)==0);
 ck(fusion_c_magnetic_push(1,1,-.001,&b,field,&f,&back,mid)==0);
 for(int j=0;j<3;++j){ck(std::abs(back.x_m[j]-s.x_m[j])<1e-10);ck(std::abs(back.u_m_s[j]-s.u_m_s[j])<1e-8);}
 Context zero{{0,0,0}};ck(fusion_c_magnetic_push(1,1,.001,&s,field,&zero,&b,mid)==0);
 double g=std::sqrt(1+std::pow(norm(s.u_m_s)/c,2));
 for(int j=0;j<3;++j){ck(b.u_m_s[j]==s.u_m_s[j]);ck(std::abs(b.x_m[j]-s.x_m[j]-.001*s.u_m_s[j]/g)<1e-10);ck(zero.at[j]==mid[j]);}
 f.status=1234;ck(fusion_c_magnetic_push(1,1,.001,&s,field,&f,&b,mid)==1234);ck(norm(b.x_m)==0&&norm(b.u_m_s)==0&&norm(mid)==0);
 int old=f.calls;ck(fusion_c_magnetic_push(1,0,.001,&s,field,&f,&b,mid)==0&&f.calls==old);
 ck(fusion_c_magnetic_push(1,1,0,&s,field,&f,&b,mid)==0&&f.calls==old);
 ck(fusion_c_magnetic_push(0,1,1,&s,field,&f,&b,mid)==PB11_STATUS_OUT_OF_RANGE);
 ck(fusion_c_magnetic_push(1,1,1,&s,nullptr,&f,&b,mid)==PB11_STATUS_INVALID_ARGUMENT);
 // Parallel component is preserved for an oblique magnetic field.
 f.status=0;double dot0=0,dot1=0;
 ck(fusion_c_magnetic_push(1,1,.3,&s,field,&f,&b,mid)==0);
 for(int j=0;j<3;++j){dot0+=s.u_m_s[j]*f.B[j];dot1+=b.u_m_s[j]*f.B[j];}
 ck(std::abs(dot1-dot0)<1e-8);
 double pmid[3];ck(fusion_c_magnetic_midpoint(&s,.3,pmid)==0);
 for(int j=0;j<3;++j)ck(pmid[j]==f.at[j]);
 for(double angle:{1e-12,1.,1000.}){
  Context axial{{0,0,1}};fusion_orbit_state_v1 a{{0,0,0},{1e7,0,2e7}};
  double ga=std::sqrt(1+std::pow(norm(a.u_m_s)/c,2));
  ck(fusion_c_magnetic_push(1,1,angle*ga,&a,field,&axial,&b,mid)==0);
  ck(std::abs(b.u_m_s[0]/1e7-std::cos(angle))<3e-13);
  ck(std::abs(b.u_m_s[1]/1e7+std::sin(angle))<3e-13);
  ck(b.u_m_s[2]==2e7&&std::abs(norm(b.u_m_s)/norm(a.u_m_s)-1)<3e-15);
 }
 s.u_m_s[0]=std::numeric_limits<double>::quiet_NaN();
 ck(fusion_c_magnetic_push(1,1,1,&s,field,&f,&b,mid)==PB11_STATUS_INVALID_ARGUMENT);
 ck(norm(b.x_m)==0&&norm(b.u_m_s)==0&&norm(mid)==0);
 ck(fusion_c_magnetic_midpoint(&s,1,mid)==PB11_STATUS_INVALID_ARGUMENT&&norm(mid)==0);
 std::cout<<"uniform-B sign, relativistic energy, position convergence, reversal, zero field and errors passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
