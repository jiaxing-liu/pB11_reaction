#ifndef FUSION_BIRTH_FLOOR_PRECISE_INTERNAL_H
#define FUSION_BIRTH_FLOOR_PRECISE_INTERNAL_H
#include "fusion_birth_floor.h"
#include <array>
#include <cmath>
namespace fusion_detail {
struct precise_floor_ledger {
 std::array<long double,6> number{},physical_energy{},mapped_energy{},correction{};
 long double borrowed=0;
};
// Private extended-precision amounts. No platform-dependent type crosses C ABI.
// No rate or packet is discarded here; caller performs source conversion once.
inline int project_floor_precise(long double center,long double temperature,
 long double reservoir,double max_center,double max_fraction,
 const std::array<long double,6>& number,
 const std::array<long double,6>& energy,precise_floor_ledger& output){
 output={};
 if(!std::isfinite(center)||!std::isfinite(temperature)||!std::isfinite(reservoir)||
    !std::isfinite(max_center)||!std::isfinite(max_fraction))return PB11_STATUS_INVALID_ARGUMENT;
 if(center<=0||temperature<=0||reservoir<=0||max_center<0||max_center>1||
    max_fraction<0||max_fraction>1)return PB11_STATUS_OUT_OF_RANGE;
 precise_floor_ledger candidate{};bool active=false;
 for(int i=0;i<6;++i){
  if(!std::isfinite(number[i])||!std::isfinite(energy[i])||number[i]<0||energy[i]<0)
   return PB11_STATUS_INVALID_ARGUMENT;
  long double mapped=number[i]*center;
  if(!std::isfinite(mapped))return PB11_STATUS_NUMERICAL_FAILURE;
  if((number[i]==0&&energy[i]!=0)||energy[i]>mapped)return PB11_STATUS_OUT_OF_RANGE;
  active=active||number[i]>0;
  candidate.number[i]=number[i];candidate.physical_energy[i]=energy[i];
  candidate.mapped_energy[i]=mapped;candidate.correction[i]=energy[i]-mapped;
  candidate.borrowed+=mapped-energy[i];
 }
 if(active&&center>static_cast<long double>(max_center)*temperature)
  return PB11_STATUS_OUT_OF_RANGE;
 if(!std::isfinite(candidate.borrowed)||candidate.borrowed>max_fraction*reservoir||
    candidate.borrowed>=reservoir)return PB11_STATUS_OUT_OF_RANGE;
 output=candidate;return PB11_STATUS_OK;
}
}
#endif
