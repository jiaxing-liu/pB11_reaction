// Test-only GNU --wrap fault injection. This does not describe a natural
// physical failure domain and introduces no production hook or library state.
// Link this executable alone with -Wl,--wrap=fusion_c_three_equal_sequential_cm.
#include "fusion_thermal_birth_budget.h"
#include "fusion_reaction_event.h"
#include "fusion_products.h" // Authoritative sequential function declaration.
#include "fusion_rate_model.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>

namespace {
unsigned calls=0,successful_real_calls=0;
}
extern "C" int __real_fusion_c_three_equal_sequential_cm(
 double product_mass_kg,double available_energy_J,
 double intermediate_relative_energy_J,double cos_theta,
 fusion_three_body_cm_v1*out);
extern "C" int __wrap_fusion_c_three_equal_sequential_cm(
 double product_mass_kg,double available_energy_J,
 double intermediate_relative_energy_J,double cos_theta,
 fusion_three_body_cm_v1*out){
 ++calls;
 if(calls==5)return PB11_STATUS_NUMERICAL_FAILURE;
 int status=__real_fusion_c_three_equal_sequential_cm(product_mass_kg,
  available_energy_J,intermediate_relative_energy_J,cos_theta,out);
 if(status==PB11_STATUS_OK)++successful_real_calls;
 return status;
}
template<class T>bool cleared(const T&value){
 const auto*p=reinterpret_cast<const unsigned char*>(&value);
 return std::all_of(p,p+sizeof(value),[](unsigned char byte){return byte==0;});
}
int main(){
 constexpr double kev=1.602176634e-16,mev=1.602176634e-13;
 constexpr int cells=8;
 const std::array<double,cells+1> edges{{0,.1*mev,.2*mev,.4*mev,.8*mev,1.2*mev,1.6*mev,2*mev,2.4*mev}};
 fusion_thermal_birth_options_v1 o{};
 o.relative_max_J=2*mev;o.cm_max_kT=40;o.ground_state_q_J=91.84*kev;
 o.cutoff_J=.001*mev;o.l1_fraction=.76;o.narrow_peak_fraction=.051;
 o.continuum_peak_scale=1;o.continuation=FUSION_ENDPOINT_S;
 o.pb_low=FUSION_PB_LOW_TB;o.remainder_policy=0;o.broad_mode=13;
 o.relative_order=4;o.cm_order=4;o.nq=4;o.ncos=4;
 fusion_thermal_birth_budget_policy_v1 policy{1,4,4,0,.01,.01};
 std::array<double,7*cells> birth;birth.fill(123.5);
 fusion_thermal_birth_v1 out;fusion_thermal_birth_budget_report_v1 report;
 std::memset(&out,0x5a,sizeof(out));std::memset(&report,0x5a,sizeof(report));
 const int status=fusion_c_pb_thermal_birth_budget_grid(.5*kev,&o,&policy,
  cells,edges.data(),birth.data(),&out,&report);
 // Four angular nodes of the first parent execute real sequential kinematics;
 // emit_products maps each returned event immediately. Failure occurs during
 // the next parent's emission, after real prepare_weights and internal mapping.
 if(status!=PB11_STATUS_NUMERICAL_FAILURE||calls!=5||successful_real_calls!=4){
  std::cerr<<"fault path mismatch: status="<<status<<" calls="<<calls
   <<" successful_real_calls="<<successful_real_calls<<'\n';return 1;
 }
 if(!cleared(birth)||!cleared(out)||!cleared(report)){
  std::cerr<<"partial emission failure published external output\n";return 2;
 }
 std::cout<<"PASS test-only fifth-call emission failure: four real sequential calls, all external outputs cleared; not natural physical-domain coverage\n";
 return 0;
}
