#include "fusion_fast_moments.h"
#include "fusion_nuclear_data.h"
#include <cmath>
#include <limits>
namespace {
bool store(long double v,double& x) {
 if(!std::isfinite(v)||v>std::numeric_limits<double>::max())return false;
 x=static_cast<double>(v);return v==0||x!=0;
}
}
extern "C" int fusion_c_fast_moments(int cells,const double* edges,const double* s,
 const double* t,fusion_fast_moments_v1* out) {
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(cells<1||cells>100000||!edges||!s||!t)return PB11_STATUS_INVALID_ARGUMENT;
 for(int i=0;i<=cells;++i)
  if(!std::isfinite(edges[i])||edges[i]<0||(i&&edges[i]<=edges[i-1]))return PB11_STATUS_INVALID_ARGUMENT;
 fusion_fast_moments_v1 candidate{};long double q=0,q2=0,u=0;
 for(int id=0;id<6;++id) {
  fusion_nuclear_mass_v1 mass{};
  int status=fusion_c_nuclear_mass(id,&mass);if(status)return status;
  long double number=0,energy=0;
  for(int i=0;i<cells;++i) {
   int k=id*cells+i;
   if(!std::isfinite(s[k])||!std::isfinite(t[k])||s[k]<0||t[k]<0)return PB11_STATUS_INVALID_ARGUMENT;
   long double pop=static_cast<long double>(s[k])+t[k];
   number+=pop;energy+=pop*(static_cast<long double>(edges[i])+edges[i+1])/2;
  }
  if(!store(number,candidate.number_m3[id])||!store(energy,candidate.energy_J_m3[id]))return PB11_STATUS_NUMERICAL_FAILURE;
  q+=number*mass.nuclear_charge;q2+=number*mass.nuclear_charge*mass.nuclear_charge;u+=energy;
 }
 if(!store(q,candidate.charge_number_m3)||!store(q2,candidate.charge_squared_number_m3)||
    !store((2.L/3)*u,candidate.pressure_Pa))return PB11_STATUS_NUMERICAL_FAILURE;
 *out=candidate;return PB11_STATUS_OK;
}
