#include "fusion_magnetic_push.h"
#include <cmath>
#include <limits>
namespace {
using L=long double;
constexpr L c=299792458.L;
bool put(L x,double&v){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;v=double(x);return x==0||v!=0;}
L gamma(const L*u){return std::hypot(1.L,std::hypot(std::hypot(u[0]/c,u[1]/c),u[2]/c));}
}
extern "C" int fusion_c_magnetic_midpoint(const fusion_orbit_state_v1*in,double dt,double*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 for(int j=0;j<3;++j)out[j]=0;
 if(!in||!std::isfinite(dt))return PB11_STATUS_INVALID_ARGUMENT;
 L u[3];for(int j=0;j<3;++j){if(!std::isfinite(in->x_m[j])||!std::isfinite(in->u_m_s[j]))return PB11_STATUS_INVALID_ARGUMENT;u[j]=in->u_m_s[j];}
 L g=gamma(u);double mid[3];
 for(int j=0;j<3;++j)if(!put(L(in->x_m[j])+L(dt)*u[j]/(2*g),mid[j]))return PB11_STATUS_NUMERICAL_FAILURE;
 for(int j=0;j<3;++j)out[j]=mid[j];
 return 0;
}
extern "C" int fusion_c_magnetic_push(double mass,double charge,double dt,
 const fusion_orbit_state_v1*in,fusion_magnetic_field_callback_v1 field,void*context,
 fusion_orbit_state_v1*out,double*mid){
 if(out)*out={};
 if(mid)for(int j=0;j<3;++j)mid[j]=0;
 if(!out||!mid)return PB11_STATUS_NULL_OUTPUT;
 if(!in||!field||!std::isfinite(mass)||!std::isfinite(charge)||!std::isfinite(dt))return PB11_STATUS_INVALID_ARGUMENT;
 if(mass<=0)return PB11_STATUS_OUT_OF_RANGE;
 L u[3];for(int j=0;j<3;++j){if(!std::isfinite(in->x_m[j])||!std::isfinite(in->u_m_s[j]))return PB11_STATUS_INVALID_ARGUMENT;u[j]=in->u_m_s[j];}
 L g=gamma(u);double midpoint[3];
 int midpoint_status=fusion_c_magnetic_midpoint(in,dt,midpoint);
 if(midpoint_status)return midpoint_status;
 if(dt==0){*out=*in;for(int j=0;j<3;++j)mid[j]=midpoint[j];return 0;}
 double B[3]={};
 if(charge!=0){for(double&v:B)v=std::numeric_limits<double>::quiet_NaN();int rc=field(context,midpoint,B);if(rc)return rc;}
 for(double v:B)if(!std::isfinite(v))return PB11_STATUS_NUMERICAL_FAILURE;
 L norm=std::hypot(std::hypot(L(B[0]),L(B[1])),L(B[2]));
 if(norm!=0&&charge!=0){
  L b[3]={B[0]/norm,B[1]/norm,B[2]/norm};
  L angle=(L(charge)/mass)*L(dt)*norm/g;
  if(!std::isfinite(angle))return PB11_STATUS_NUMERICAL_FAILURE;
  L co=std::cos(angle),si=std::sin(angle),dot=u[0]*b[0]+u[1]*b[1]+u[2]*b[2];
  L cross[3]={u[1]*b[2]-u[2]*b[1],u[2]*b[0]-u[0]*b[2],u[0]*b[1]-u[1]*b[0]};
  for(int j=0;j<3;++j)u[j]=dot*b[j]+(u[j]-dot*b[j])*co+cross[j]*si;
 }
 fusion_orbit_state_v1 trial{};
 for(int j=0;j<3;++j)if(!put(u[j],trial.u_m_s[j]))return PB11_STATUS_NUMERICAL_FAILURE;
 // The second drift uses the actual rounded proper velocity returned to caller.
 for(int j=0;j<3;++j)u[j]=trial.u_m_s[j];
 g=gamma(u);
 for(int j=0;j<3;++j)if(!put(L(midpoint[j])+L(dt)*u[j]/(2*g),trial.x_m[j]))return PB11_STATUS_NUMERICAL_FAILURE;
 *out=trial;for(int j=0;j<3;++j)mid[j]=midpoint[j];return 0;
}
