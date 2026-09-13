#include "fusion_alpha_spectrum.h"
#include "fusion_alpha_events.h"
#include "fusion_products.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
constexpr double mev=1.602176634e-13;
static int spectrum_impl(int mode,int policy,double A,double cutoff,
 double k,double phase,int nq,int nc,int cells,const double *edges,double *birth,
 fusion_alpha_spectrum_v1 *out){
 if(out)*out={};
 if(cells>0 && birth)std::fill(birth,birth+cells,0.);
 if(!birth || !out)return PB11_STATUS_NULL_OUTPUT;
 if(!edges || cells<1 || nq<4 || nc<4 || nq>1024 || nc>1024 ||
  (policy!=0 && policy!=1) || (mode!=1 && mode!=2 && mode!=3 && mode!=13) || !std::isfinite(A) ||
  !std::isfinite(k) || !std::isfinite(phase) || !std::isfinite(cutoff))return PB11_STATUS_INVALID_ARGUMENT;
 if(A<=0 || A>12*mev || k<0 || k>1 || cutoff<.001*mev || cutoff>.01*mev)return PB11_STATUS_OUT_OF_RANGE;
 for(int i=0;i<=cells;++i)
  if(!std::isfinite(edges[i]) || edges[i]<0 || (i>0 && edges[i]<=edges[i-1]))
   return PB11_STATUS_INVALID_ARGUMENT;
 try{
  fusion_detail::AlphaEvents events;
  int status=fusion_detail::alpha_events(mode,policy,A,cutoff,k,phase,nq,nc,events);
  if(status)return status;
  std::vector<double> packet_e,packet_w;packet_e.reserve(3*events.events.size());packet_w.reserve(3*events.events.size());
  for(const auto &event:events.events)for(double e:event.energy_J){packet_e.push_back(e);packet_w.push_back(event.weight);}
  std::vector<double> mapped(cells);fusion_birth_mapping_v1 ledger{};
  // Unit reaction rate1 m^-3 s^-1; divide all outputs by that unit rate.
  status=fusion_c_map_birth_packets(cells,edges,static_cast<int>(packet_e.size()),packet_e.data(),packet_w.data(),mapped.data(),&ledger);if(status)return status;
  fusion_alpha_spectrum_v1 result{};
  result.mapped_number=ledger.mapped_number_m3_s;result.mapped_energy_J=ledger.mapped_energy_W_m3;
  result.below_number=ledger.below_number_m3_s;result.below_energy_J=ledger.below_energy_W_m3;
  result.above_number=ledger.above_number_m3_s;result.above_energy_J=ledger.above_energy_W_m3;
  result.number_residual=ledger.input_number_m3_s-3;result.energy_residual_J=ledger.input_energy_W_m3-A;
  result.normalization_J2=events.normalization_J2;
  result.l1_normalization_J2=mode==13?events.l1_normalization_J2:0;
  result.l3_normalization_J2=mode==13?events.l3_normalization_J2:0;
  result.quadrature_events=nq*nc;result.pruned_events=events.pruned_events;
  if(std::abs(result.number_residual)>3e-12 || std::abs(result.energy_residual_J)>1e-12*A ||
   !std::isfinite(result.normalization_J2) || result.normalization_J2==0)return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(mapped.begin(),mapped.end(),birth);*out=result;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}

} // namespace

extern "C" int fusion_c_alpha_spectrum_grid(int mode,double A,double cutoff,
 double k,double phase,int nq,int nc,int cells,const double *edges,double *birth,
 fusion_alpha_spectrum_v1 *out){
 return spectrum_impl(mode,0,A,cutoff,k,phase,nq,nc,cells,edges,birth,out);
}
extern "C" int fusion_c_alpha_spectrum_model_grid(int mode,int policy,double A,
 double cutoff,double k,double phase,int nq,int nc,int cells,const double *edges,
 double *birth,fusion_alpha_spectrum_v1 *out){
 return spectrum_impl(mode,policy,A,cutoff,k,phase,nq,nc,cells,edges,birth,out);
}
