#include "fusion_coupled_thermal.h"
#include "fusion_thermal_burn.h"
#include "fusion_nuclear_data.h"
#include "fusion_two_component.h"
#include "fusion_handoff.h"
#include "fusion_rate_model.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using R=long double;
constexpr int BAD=PB11_STATUS_INVALID_ARGUMENT,NUM=PB11_STATUS_NUMERICAL_FAILURE;
constexpr double MeV=1.602176634e-13;
bool finite_value(double x){return std::isfinite(x);}
bool nonnegative(double x){return finite_value(x)&&x>=0;}
bool put(R x,double& d){if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;d=double(x);return true;}
bool close(std::initializer_list<R> terms,R roundoff=0){R value=0,scale=0;for(R x:terms){if(!std::isfinite(x))return false;value+=x;scale+=std::abs(x);}return std::abs(value)<=1e-10L*scale+roundoff;}
template<class F>void ledger_fields(fusion_source_ledger_v1&l,F f){
 for(double&x:l.events_m3)f(x);
 for(double&x:l.nuclear_born_number_m3)f(x);
 for(double&x:l.nuclear_born_energy_J_m3)f(x);
 for(double&x:l.external_born_number_m3)f(x);
 for(double&x:l.external_born_energy_J_m3)f(x);
 for(double&x:l.thermal_consumed_number_m3)f(x);
 for(double&x:l.thermal_consumed_energy_J_m3)f(x);
 for(double&x:l.fast_consumed_number_m3)f(x);
 for(double&x:l.fast_consumed_energy_J_m3)f(x);
 for(double&x:l.escaped_number_m3)f(x);
 for(double&x:l.escaped_energy_J_m3)f(x);
 for(double&x:l.handed_off_number_m3)f(x);
 for(double&x:l.handed_off_energy_J_m3)f(x);
 for(double&x:l.heat_to_bath_J_m3)f(x);
 f(l.neutron_number_m3);f(l.neutron_energy_J_m3);
}
int validate_options(const fusion_coupled_thermal_options_v1&o){
 const auto&b=o.birth;
 for(double x:{b.relative_max_J,b.cm_max_kT,b.ground_state_q_J,b.cutoff_J,b.l1_fraction,b.relative_phase,b.narrow_peak_fraction,b.continuum_peak_scale,o.max_source_rate_error,o.max_source_debit_error,o.handoff_max_L1,o.handoff_max_mean_error})if(!finite_value(x))return BAD;
 if(b.relative_max_J<=0||b.cm_max_kT<8||b.cm_max_kT>80||b.ground_state_q_J<0||b.cutoff_J<.001*MeV||b.cutoff_J>.01*MeV||b.l1_fraction<0||b.l1_fraction>1||b.narrow_peak_fraction<0||b.narrow_peak_fraction>1||b.continuum_peak_scale<0||b.continuum_peak_scale>2||o.max_source_rate_error<0||o.max_source_rate_error>1||o.max_source_debit_error<0||o.max_source_debit_error>1||o.handoff_max_L1<0||o.handoff_max_L1>2||o.handoff_max_mean_error<0||o.handoff_max_mean_error>1)return PB11_STATUS_OUT_OF_RANGE;
 if(b.relative_order<4||b.relative_order>64||b.cm_order<4||b.cm_order>32||b.nq<4||b.nq>1024||b.ncos<4||b.ncos>1024||b.remainder_policy<0||b.remainder_policy>2||(b.broad_mode!=1&&b.broad_mode!=3&&b.broad_mode!=13)||b.fsci_policy<0||b.fsci_policy>1||o.handoff_enabled<0||o.handoff_enabled>1)return BAD;
 for(int x:o.channels)if(x!=0&&x!=1)return BAD;
 double sigma=0;int st=fusion_c_cross_section_model(0,b.continuation,b.pb_low,MeV,&sigma);if(st)return st;
 fusion_nuclear_channel_v1 pb{};st=fusion_c_nuclear_channel(0,&pb);if(st)return st;
 return b.ground_state_q_J<=pb.q_J?PB11_STATUS_OK:PB11_STATUS_OUT_OF_RANGE;
}
}
extern "C" int fusion_c_coupled_thermal_trial(double dt,const fusion_coupled_thermal_options_v1*op,
 int n,const double*edges,const double*thermal,double Ue,double Ui,double ne,const double*charge2,
 int ninert,const fusion_inert_ion_v1*inert,const double*logs,const double*old_s,const double*old_t,
 const double*external,const double*escape,double*new_thermal,double*new_s,double*new_t,fusion_coupled_thermal_v1*out){
 if(out)*out={};
 if(new_thermal)std::fill(new_thermal,new_thermal+6,0.);
 if(n>0&&n<=100000){if(new_s)std::fill(new_s,new_s+6*n,0.);if(new_t)std::fill(new_t,new_t+6*n,0.);}
 if(!out||!new_thermal||!new_s||!new_t)return PB11_STATUS_NULL_OUTPUT;
 if(!op||!edges||!thermal||!charge2||!logs||!old_s||!old_t||!external||!escape||n<1||n>100000||ninert<0||ninert>32||(ninert&&!inert))return BAD;
 if(!finite_value(dt)||!finite_value(Ue)||!finite_value(Ui)||!finite_value(ne))return BAD;
 if(dt<=0||Ue<=0||Ui<=0||ne<=0)return PB11_STATUS_OUT_OF_RANGE;
 int st=validate_options(*op);if(st)return st;
 const int nb=7+ninert;
 for(int j=0;j<=n;++j)if(!nonnegative(edges[j])||(j&&edges[j]<=edges[j-1]))return BAD;
 for(int i=0;i<6*n;++i)if(!nonnegative(old_s[i])||!nonnegative(old_t[i])||!nonnegative(external[i])||!nonnegative(escape[i]))return BAD;
 for(int i=0;i<6*nb;++i)if(!finite_value(logs[i])||logs[i]<=0)return PB11_STATUS_OUT_OF_RANGE;
 try{
  std::array<fusion_nuclear_mass_v1,8> mass{};for(int i=0;i<8;++i){st=fusion_c_nuclear_mass(i,&mass[i]);if(st)return st;}
  R Npool=0,Ninert=0;for(int i=0;i<6;++i){if(!nonnegative(thermal[i])||!nonnegative(charge2[i])||charge2[i]>mass[i].nuclear_charge*mass[i].nuclear_charge)return BAD;Npool+=thermal[i];}
  for(int j=0;j<ninert;++j){if(!nonnegative(inert[j].density_m3)||!finite_value(inert[j].mass_kg)||inert[j].mass_kg<=0||!nonnegative(inert[j].mean_charge_squared))return BAD;Ninert+=inert[j].density_m3;}
  Npool+=Ninert;if(Npool<=0)return PB11_STATUS_OUT_OF_RANGE;
  double Ti=0,Te=0;if(!put(R(Ui)/(1.5L*Npool),Ti)||!put(R(Ue)/(1.5L*ne),Te)||Ti<=0||Te<=0)return NUM;
  std::vector<R> centers(n);for(int j=0;j<n;++j)centers[j]=(R(edges[j])+edges[j+1])/2;
  std::array<R,6> oldN{},oldE{};for(int i=0;i<6;++i)for(int j=0;j<n;++j){R z=R(old_s[i*n+j])+old_t[i*n+j];oldN[i]+=z;oldE[i]+=z*centers[j];}
  std::array<fusion_nuclear_channel_v1,5> reactions{};
  std::array<fusion_thermal_birth_v1,5> source{};
  std::array<std::vector<double>,5> grids;
  double rates[5]{},mean_a[5]{},mean_b[5]{},oldUi[6]{},trialNi[6]{},trialUi[6]{};
  fusion_coupled_thermal_v1 result{};auto&l=result.ledger;
  for(int i=0;i<6;++i)if(!put(R(1.5L)*Ti*thermal[i],oldUi[i]))return NUM;
  for(int ch=0;ch<5;++ch){st=fusion_c_nuclear_channel(ch,&reactions[ch]);if(st)return st;
   if(!op->channels[ch]||thermal[reactions[ch].reactant_ids[0]]==0||thermal[reactions[ch].reactant_ids[1]]==0)continue;
   auto options=op->birth;if(ch!=0)options.pb_low=FUSION_PB_LOW_TB;
   grids[ch].resize(7*n);st=fusion_c_thermal_birth_grid(ch,Ti,&options,n,edges,grids[ch].data(),&source[ch]);if(st)return st;
   const auto&r=source[ch];result.max_source_rate_discrepancy=std::max(result.max_source_rate_discrepancy,std::abs(r.relative_rate_discrepancy));
   result.max_source_debit_discrepancy=std::max(result.max_source_debit_discrepancy,r.relative_reactant_energy_discrepancy);
   if(std::abs(r.relative_rate_discrepancy)>op->max_source_rate_error||r.relative_reactant_energy_discrepancy>op->max_source_debit_error)return NUM;
   for(int i=0;i<6;++i)if(r.below_number_m3_s[i]>0||r.above_number_m3_s[i]>0)return PB11_STATUS_OUT_OF_RANGE;
   rates[ch]=r.reactivity_m3_s;
   if(rates[ch]>0){mean_a[ch]=r.reactant_energy_moment_J_m3_s[0]/rates[ch];mean_b[ch]=r.reactant_energy_moment_J_m3_s[1]/rates[ch];}
  }
  fusion_thermal_burn_v1 burn{};st=fusion_c_thermal_burn_trial(dt,thermal,oldUi,rates,mean_a,mean_b,trialNi,trialUi,&burn);if(st)return st;
  R ionU=Ui,electronU=Ue;Npool=Ninert;for(int i=0;i<6;++i){l.thermal_consumed_number_m3[i]=burn.reactant_removed_m3[i];l.thermal_consumed_energy_J_m3[i]=burn.reactant_removed_energy_J_m3[i];ionU-=burn.reactant_removed_energy_J_m3[i];Npool+=trialNi[i];}
  if(ionU<=0||Npool<=0||!put(ionU/(1.5L*Npool),Ti)||Ti<=0)return NUM;
  std::vector<double> birth(6*n),s(6*n),t(6*n);R Q=0,neutronN=0,neutronE=0;
  for(int ch=0;ch<5;++ch){l.events_m3[ch]=burn.events_m3[ch];Q+=R(l.events_m3[ch])*reactions[ch].q_J;}
  for(int i=0;i<6;++i){R bornN=0,bornE=0,extN=0,extE=0;
   for(int j=0;j<n;++j){R amount=0;
    for(int ch=0;ch<5;++ch)if(burn.events_m3[ch]>0){if(rates[ch]<=0)return NUM;amount+=R(grids[ch][i*n+j])/rates[ch]*burn.events_m3[ch];}
    bornN+=amount;bornE+=amount*centers[j];R ex=R(dt)*external[i*n+j];extN+=ex;extE+=ex*centers[j];
    if(!put((amount+ex)/dt,birth[i*n+j]))return NUM;
   }
   if(!put(bornN,l.nuclear_born_number_m3[i])||!put(bornE,l.nuclear_born_energy_J_m3[i])||!put(extN,l.external_born_number_m3[i])||!put(extE,l.external_born_energy_J_m3[i]))return NUM;
  }
  for(int ch=0;ch<5;++ch)if(burn.events_m3[ch]>0){R N=R(source[ch].below_number_m3_s[6])+source[ch].above_number_m3_s[6],E=R(source[ch].below_energy_J_m3_s[6])+source[ch].above_energy_J_m3_s[6];for(int j=0;j<n;++j){N+=grids[ch][6*n+j];E+=R(grids[ch][6*n+j])*centers[j];}neutronN+=N/rates[ch]*burn.events_m3[ch];neutronE+=E/rates[ch]*burn.events_m3[ch];}
  if(!put(neutronN,l.neutron_number_m3)||!put(neutronE,l.neutron_energy_J_m3))return NUM;
  std::vector<fusion_maxwellian_bath_v1> baths(nb);
  baths[0]={ne,mass[7].mass_kg,1,Te,1};
  for(int j=0;j<6;++j)baths[j+1]={trialNi[j],mass[j].mass_kg,charge2[j],Ti,1};
  for(int j=0;j<ninert;++j)baths[j+7]={inert[j].density_m3,inert[j].mass_kg,inert[j].mean_charge_squared,Ti,1};
  std::vector<double> temperatures(nb),diffusion(nb*(n-1)),transfer(n),zero(n),heat(nb),after(n);
  for(int b=0;b<nb;++b)temperatures[b]=baths[b].kT_J;
  for(int i=0;i<6;++i){bool active=oldN[i]>0||l.nuclear_born_number_m3[i]>0||l.external_born_number_m3[i]>0;if(!active)continue;
   std::fill(transfer.begin(),transfer.end(),0.);
   for(int b=0;b<nb;++b){baths[b].coulomb_log=logs[i*nb+b];
    for(int j=0;j<n-1;++j){fusion_coulomb_energy_v1 coefficient{};st=fusion_c_coulomb_energy(edges[j+1],mass[i].mass_kg,mass[i].nuclear_charge,&baths[b],&coefficient);if(st)return st;diffusion[b*(n-1)+j]=coefficient.diffusion_J2_s;}
    if(b>0)for(int j=0;j<n;++j){double lambda=0;st=fusion_c_coulomb_transfer_rate(double(centers[j]),mass[i].mass_kg,mass[i].nuclear_charge,&baths[b],&lambda);if(st)return st;if(!put(R(transfer[j])+lambda,transfer[j]))return NUM;}
   }
   fusion_two_component_ledger_v1 fp{};
   st=fusion_c_two_component_trial(n,nb,dt,edges,old_s+i*n,old_t+i*n,temperatures.data(),diffusion.data(),birth.data()+i*n,zero.data(),escape+i*n,transfer.data(),s.data()+i*n,t.data()+i*n,heat.data(),&fp);if(st)return st;
   l.escaped_number_m3[i]=fp.total.escaped_number_m3;l.escaped_energy_J_m3[i]=fp.total.escaped_energy_J_m3;
   electronU+=heat[0];R extra=0;for(int b=0;b<nb;++b){if(b<7)l.heat_to_bath_J_m3[i*7+b]=heat[b];else extra+=heat[b];if(b>0)ionU+=heat[b];}
   if(!put(extra,result.inert_ion_heat_J_m3[i]))return NUM;
  }
  if(electronU<=0||ionU<=0)return NUM;
  if(op->handoff_enabled)for(int i=0;i<6;++i){R N=0,E=0;for(int j=0;j<n;++j){N+=t[i*n+j];E+=R(t[i*n+j])*centers[j];}if(N==0)continue;
   double mixedTi=0;if(!put((ionU+E)/(1.5L*(Npool+N)),mixedTi)||mixedTi<=0)return NUM;
   fusion_handoff_ledger_v1 handoff{};int projected=0;
   st=fusion_c_maxwellian_handoff_trial(n,mixedTi,op->handoff_max_L1,op->handoff_max_mean_error,edges,t.data()+i*n,after.data(),&projected,&handoff);if(st)return st;
   result.handoff_L1[i]=handoff.distribution_L1;result.handoff_mean_error[i]=handoff.relative_mean_energy_error;result.handoff_projected[i]=projected;
   if(projected){l.handed_off_number_m3[i]=handoff.fluid_number_m3;l.handed_off_energy_J_m3[i]=handoff.fluid_energy_J_m3;
    if(!put(R(l.heat_to_bath_J_m3[i*7+i+1])+handoff.bath_energy_correction_J_m3,l.heat_to_bath_J_m3[i*7+i+1])||!put(R(trialNi[i])+handoff.fluid_number_m3,trialNi[i]))return NUM;
    ionU+=R(handoff.fluid_energy_J_m3)+handoff.bath_energy_correction_J_m3;Npool+=handoff.fluid_number_m3;
    std::copy(after.begin(),after.end(),t.begin()+i*n);
   }
  }
  if(!put(electronU,result.electron_energy_J_m3)||!put(ionU,result.ion_energy_J_m3))return NUM;
  // Independent accounting of reaction, kinetic, and complete local inventories.
  R nuclearE=neutronE,removedE=0,totalOld=R(Ue)+Ui,totalNew=R(result.electron_energy_J_m3)+result.ion_energy_J_m3;
  R escapedE=0,externalE=0;
  for(int i=0;i<6;++i){R N=0,E=0,H=result.inert_ion_heat_J_m3[i];for(int b=0;b<7;++b)H+=l.heat_to_bath_J_m3[i*7+b];
   for(int j=0;j<n;++j){R z=R(s[i*n+j])+t[i*n+j];N+=z;E+=z*centers[j];}
   if(!close({N,-oldN[i],-R(l.nuclear_born_number_m3[i]),-R(l.external_born_number_m3[i]),R(l.escaped_number_m3[i]),R(l.handed_off_number_m3[i])})||
      !close({E,-oldE[i],-R(l.nuclear_born_energy_J_m3[i]),-R(l.external_born_energy_J_m3[i]),R(l.escaped_energy_J_m3[i]),R(l.handed_off_energy_J_m3[i]),H}))return NUM;
   R dN=R(trialNi[i])-thermal[i]+N-oldN[i]-l.nuclear_born_number_m3[i]+l.thermal_consumed_number_m3[i]-l.external_born_number_m3[i]+l.escaped_number_m3[i];
   // Difference-of-inventory residuals inherit rounding from BOTH pools,
   // including a seeded kinetic species with no thermal counterpart.
   R round=16*std::numeric_limits<double>::epsilon()*(R(thermal[i])+trialNi[i]+oldN[i]+N);
   if(!put(dN,result.particle_residual_m3[i])||!close({R(trialNi[i])-thermal[i],N-oldN[i],-R(l.nuclear_born_number_m3[i]),R(l.thermal_consumed_number_m3[i]),-R(l.external_born_number_m3[i]),R(l.escaped_number_m3[i])},round))return NUM;
   nuclearE+=l.nuclear_born_energy_J_m3[i];removedE+=l.thermal_consumed_energy_J_m3[i];totalOld+=oldE[i];totalNew+=E;escapedE+=l.escaped_energy_J_m3[i];externalE+=l.external_born_energy_J_m3[i];
  }
  if(!close({nuclearE,-removedE,-Q}))return NUM;
  R residual=totalNew-totalOld+neutronE+escapedE-externalE-Q;
  R round=16*std::numeric_limits<double>::epsilon()*(std::abs(totalOld)+std::abs(totalNew));
  if(!put(residual,result.energy_residual_J_m3)||!close({totalNew-totalOld,neutronE,escapedE,-externalE,-Q},round))return NUM;
  bool valid=true;ledger_fields(l,[&](double x){if(!finite_value(x))valid=false;});if(!valid)return NUM;
  std::copy(trialNi,trialNi+6,new_thermal);std::copy(s.begin(),s.end(),new_s);std::copy(t.begin(),t.end(),new_t);*out=result;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
