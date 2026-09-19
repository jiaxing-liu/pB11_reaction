#include "fusion_thermal_markers.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using L=long double;
struct Sum {L s=0,c=0;void add(L x){L t=s+x;if(std::abs(s)>=std::abs(x))c+=(s-t)+x;else c+=(x-t)+s;s=t;}L value()const{return s+c;}};
bool valid_number(double x){return std::isfinite(x);}
bool put(L x,double&v){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;v=double(x);return x==0||v!=0;}
bool near(L a,L b,L tol){return std::abs(a-b)<=tol;}
}
extern "C" int fusion_c_thermal_markers(int model,double N,double K,double mass,double V,double tol,
 int nx,const fusion_birth_spatial_node_v1*x,int na,const fusion_birth_direction_v1*a,int capacity,fusion_thermal_marker_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 if(capacity<0||capacity>1000000)return PB11_STATUS_OUT_OF_RANGE;
 std::fill_n(out,capacity,fusion_thermal_marker_v1{});
 if(model!=FUSION_THERMAL_ZERO_DRIFT_UNIFORM_VOLUME||!valid_number(N)||N<0||!valid_number(K)||K<0||
    !valid_number(mass)||mass<=0||!valid_number(V)||V<=0||!valid_number(tol)||tol<=0||tol>1e-3||!x||!a)return PB11_STATUS_INVALID_ARGUMENT;
 if(nx<1||nx>1000000||na<1||na>1000000||static_cast<long long>(nx)*na!=capacity)return PB11_STATUS_OUT_OF_RANGE;
 try{
  Sum volume,weight,m1[3],m2[3][3];
  for(int i=0;i<nx;++i){
   if(!valid_number(x[i].volume_weight_m3)||x[i].volume_weight_m3<=0)return PB11_STATUS_INVALID_ARGUMENT;
   for(double q:x[i].x_m)if(!valid_number(q))return PB11_STATUS_INVALID_ARGUMENT;
   volume.add(x[i].volume_weight_m3);
  }
  if(!near(volume.value(),V,L(tol)*V))return PB11_STATUS_INVALID_ARGUMENT;
  for(int j=0;j<na;++j){
   const auto&d=a[j];if(!valid_number(d.probability_weight)||d.probability_weight<=0)return PB11_STATUS_INVALID_ARGUMENT;
   L norm=0;for(double q:d.direction){if(!valid_number(q))return PB11_STATUS_INVALID_ARGUMENT;norm+=L(q)*q;}
   if(!near(norm,1,tol))return PB11_STATUS_INVALID_ARGUMENT;
   weight.add(d.probability_weight);
   for(int k=0;k<3;++k){m1[k].add(L(d.probability_weight)*d.direction[k]);for(int l=0;l<3;++l)m2[k][l].add(L(d.probability_weight)*d.direction[k]*d.direction[l]);}
  }
  if(!near(weight.value(),1,tol))return PB11_STATUS_INVALID_ARGUMENT;
  for(int k=0;k<3;++k){if(!near(m1[k].value(),0,tol))return PB11_STATUS_INVALID_ARGUMENT;for(int l=0;l<3;++l)if(!near(m2[k][l].value(),k==l?1.L/3:0,tol))return PB11_STATUS_INVALID_ARGUMENT;}
  constexpr L c=299792458.L;
  L proper=std::sqrt(L(K)/mass)*std::sqrt(2+L(K)/(L(mass)*c*c));
  if(!std::isfinite(proper))return PB11_STATUS_NUMERICAL_FAILURE;
  std::vector<fusion_thermal_marker_v1> result(capacity);Sum ns,es;
  for(int i=0;i<nx;++i)for(int j=0;j<na;++j){
   auto&r=result[i*na+j];L n=L(N)*(L(x[i].volume_weight_m3)/V)*a[j].probability_weight;
   if(!put(n,r.number_weight))return PB11_STATUS_NUMERICAL_FAILURE;
   r.kinetic_energy_J=K;
   for(int k=0;k<3;++k){r.initial.x_m[k]=x[i].x_m[k];if(!put(proper*a[j].direction[k],r.initial.u_m_s[k]))return PB11_STATUS_NUMERICAL_FAILURE;}
   ns.add(r.number_weight);es.add(L(r.number_weight)*K);
  }
  double energy_check=0;if(!put(L(N)*K,energy_check))return PB11_STATUS_NUMERICAL_FAILURE;
  if(!near(ns.value(),N,L(tol)*N)||!near(es.value(),L(N)*K,L(tol)*N*K))return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(result.begin(),result.end(),out);return 0;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
