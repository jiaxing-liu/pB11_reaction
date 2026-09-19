#include "fusion_prompt_reduce.h"
#include <cmath>
#include <limits>
namespace {
using L=long double;
struct Sum {L s=0,c=0;void add(L x){L t=s+x;if(std::abs(s)>=std::abs(x))c+=(s-t)+x;else c+=(x-t)+s;s=t;}L value()const{return s+c;}};
bool put(L x,double&v){if(!std::isfinite(x)||x<0||x>std::numeric_limits<double>::max())return false;v=double(x);return x==0||v>0;}
}
extern "C" int fusion_c_prompt_reduce(int n,const fusion_prompt_weight_v1*items,fusion_prompt_totals_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(n<0||n>1000000)return PB11_STATUS_OUT_OF_RANGE;
 if(n&&!items)return PB11_STATUS_INVALID_ARGUMENT;
 Sum ns[3],es[3];
 for(int i=0;i<n;++i){
  const auto&r=items[i];if(r.api_status)return r.api_status;
  if(r.outcome<0||r.outcome>2||!std::isfinite(r.number_weight)||!std::isfinite(r.kinetic_energy_J)||r.number_weight<0||r.kinetic_energy_J<0)return PB11_STATUS_INVALID_ARGUMENT;
  ns[r.outcome].add(r.number_weight);es[r.outcome].add(L(r.number_weight)*r.kinetic_energy_J);
 }
 fusion_prompt_totals_v1 result{};Sum nt,et;
 for(int k=0;k<3;++k){nt.add(ns[k].value());et.add(es[k].value());if(!put(ns[k].value(),result.number[k])||!put(es[k].value(),result.energy_J[k]))return PB11_STATUS_NUMERICAL_FAILURE;}
 L N=nt.value(),E=et.value();
 if(!put(N,result.total_number)||!put(E,result.total_energy_J))return PB11_STATUS_NUMERICAL_FAILURE;
 if(N>0){
  if(!put(ns[2].value()/N,result.number_fraction_low)||!put((ns[2].value()+ns[0].value())/N,result.number_fraction_high))return PB11_STATUS_NUMERICAL_FAILURE;
  result.number_fraction_defined=1;
 }
 if(E>0){
  if(!put(es[2].value()/E,result.energy_fraction_low)||!put((es[2].value()+es[0].value())/E,result.energy_fraction_high))return PB11_STATUS_NUMERICAL_FAILURE;
  result.energy_fraction_defined=1;
 }
 *out=result;return 0;
}
