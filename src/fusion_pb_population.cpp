#include "fusion_pb_population.h"
#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"
#include "fusion_pb_population_data.h"
#include <algorithm>
#include <cmath>
#include <limits>
extern "C" int fusion_c_pb_population(int policy,int low,double E,double peak,
 double scale,fusion_pb_population_v1 *out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(!std::isfinite(E)||!std::isfinite(peak)||!std::isfinite(scale))return PB11_STATUS_INVALID_ARGUMENT;
 if(E<0||peak<0||peak>1||scale<0||scale>2)return PB11_STATUS_OUT_OF_RANGE;
 fusion_pb_population_v1 result{};
 int status=fusion_c_cross_section_model(FUSION_PB11_3ALPHA,policy,low,E,&result.total_cross_section_m2);if(status)return status;
 constexpr double kev=1.602176634e-16;
 const long double e=static_cast<long double>(E)/kev;
 long double ratio=0;
 if(E<=400*kev){
  const long double d=e-148,narrow=1.82e4L/(d*d+2.35L*2.35L);
  const long double c0=197+(low==FUSION_PB_LOW_C0_MINUS12?-12:low==FUSION_PB_LOW_C0_PLUS12?12:0);
  const long double c1=low==FUSION_PB_LOW_NS?.240L:.269L,c2=low==FUSION_PB_LOW_NS?2.31e-4L:2.54e-4L;
  ratio=narrow/(c0+c1*e+c2*e*e+narrow);
 }
 fusion_nuclear_mass_v1 p{},boron{};
 status=fusion_c_nuclear_mass(FUSION_PROTON,&p);if(status)return status;
 status=fusion_c_nuclear_mass(FUSION_BORON11,&boron);if(status)return status;
 const long double lab_J=static_cast<long double>(E)*(1+static_cast<long double>(p.mass_kg)/boron.mass_kg);
 if(lab_J>std::numeric_limits<double>::max())return PB11_STATUS_NUMERICAL_FAILURE;
 result.equivalent_proton_lab_energy_J=static_cast<double>(lab_J);
 const long double lab=lab_J/kev;
 using namespace pb_population_data;
 auto value=[](std::size_t i){return static_cast<long double>(alpha0_mb[i])/(alpha0_mb[i]+static_cast<long double>(alpha1_mb[i]));};
 long double continuum=0;
 if(lab<=lab_keV.front())continuum=value(0);
 else if(lab>=lab_keV.back())continuum=value(19);
 else {
  auto right=std::upper_bound(lab_keV.begin(),lab_keV.end(),lab);auto j=right-lab_keV.begin();
  long double f=(lab-lab_keV[j-1])/(lab_keV[j]-lab_keV[j-1]);continuum=value(j-1)*(1-f)+value(j)*f;
 }
 continuum*=scale;
 result.narrow_fit_fraction=static_cast<double>(ratio);
 result.narrow_fit_cross_section_m2=static_cast<double>(ratio*result.total_cross_section_m2);
 result.continuum_peak_fraction=static_cast<double>(continuum);
 result.alpha0_peak_fraction=static_cast<double>(ratio*peak+(1-ratio)*continuum);
 result.narrow_remainder_fraction=static_cast<double>(ratio*(1-peak));
 result.other_remainder_fraction=static_cast<double>((1-ratio)*(1-continuum));
 result.continuum_extrapolated_fraction=(lab<lab_keV.front()||lab>lab_keV.back())?static_cast<double>(1-ratio):0.;
 *out=result;return PB11_STATUS_OK;
}
