#include "fusion_birth_floor.h"
#include <cmath>
#include <limits>
#include <initializer_list>
namespace {
using R=long double;
bool finite_value(double x){return std::isfinite(x);}
bool put(R x,double& y){
 if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;
 y=static_cast<double>(x);return true;
}
}
extern "C" int fusion_c_birth_floor_project(const fusion_birth_floor_options_v1*o,
 const double*N,const double*U,fusion_birth_floor_ledger_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(!o||!N||!U)return PB11_STATUS_INVALID_ARGUMENT;
 for(double x:{o->first_center_J,o->ion_kT_J,o->ion_energy_J_m3,
               o->max_center_over_ion_kT,o->max_ion_energy_fraction})
  if(!finite_value(x))return PB11_STATUS_INVALID_ARGUMENT;
 if(o->first_center_J<=0||o->ion_kT_J<=0||o->ion_energy_J_m3<=0||
    o->max_center_over_ion_kT<0||o->max_center_over_ion_kT>1||
    o->max_ion_energy_fraction<0||o->max_ion_energy_fraction>1)
  return PB11_STATUS_OUT_OF_RANGE;
 fusion_birth_floor_ledger_v1 result{};
 R borrowed=0,returned_borrowed=0;bool active=false;
 for(int i=0;i<6;++i){
  if(!finite_value(N[i])||!finite_value(U[i])||N[i]<0||U[i]<0)return PB11_STATUS_INVALID_ARGUMENT;
  R mapped=R(N[i])*o->first_center_J;
  if((N[i]==0&&U[i]!=0)||R(U[i])>mapped)return PB11_STATUS_OUT_OF_RANGE;
  active=active||N[i]>0;
  if(!put(mapped,result.mapped_energy_J_m3[i]))return PB11_STATUS_NUMERICAL_FAILURE;
  result.born_number_m3[i]=result.mapped_number_m3[i]=N[i];
  result.born_energy_J_m3[i]=U[i];
  borrowed+=mapped-U[i];
  const R correction=R(U[i])-result.mapped_energy_J_m3[i];
  if(!put(correction,result.ion_energy_correction_J_m3[i]))return PB11_STATUS_NUMERICAL_FAILURE;
  returned_borrowed-=result.ion_energy_correction_J_m3[i];
  R residual=R(result.mapped_energy_J_m3[i])+result.ion_energy_correction_J_m3[i]-U[i];
  R scale=std::abs(R(result.mapped_energy_J_m3[i]))+std::abs(R(result.ion_energy_correction_J_m3[i]))+U[i];
  R allowance=4*std::numeric_limits<double>::epsilon()*scale+2*R(std::numeric_limits<double>::denorm_min());
  if(std::abs(residual)>allowance||!put(residual,result.energy_residual_J_m3[i]))return PB11_STATUS_NUMERICAL_FAILURE;
 }
 if(active&&R(o->first_center_J)>R(o->max_center_over_ion_kT)*o->ion_kT_J)return PB11_STATUS_OUT_OF_RANGE;
 const R limit=R(o->max_ion_energy_fraction)*o->ion_energy_J_m3;
 if(!std::isfinite(borrowed)||!std::isfinite(returned_borrowed)||borrowed>limit||returned_borrowed>limit||
    borrowed>=R(o->ion_energy_J_m3))
  return PB11_STATUS_OUT_OF_RANGE;
 const R remaining=R(o->ion_energy_J_m3)-returned_borrowed;
 if(remaining<=0||!put(remaining,result.remaining_ion_energy_J_m3)||result.remaining_ion_energy_J_m3<=0)
  return PB11_STATUS_OUT_OF_RANGE;
 *out=result;return PB11_STATUS_OK;
}
