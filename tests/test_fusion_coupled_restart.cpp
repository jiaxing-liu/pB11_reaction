#include "fusion_coupled_thermal.h"
#include "fusion_source_state.h"
#include "fusion_nuclear_data.h"
#include <array>
#include <vector>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <algorithm>
namespace {
using Handle=std::unique_ptr<fusion_source_state_v1,decltype(&fusion_c_source_state_destroy)>;
void need(bool good,const char* message){if(!good)throw std::runtime_error(message);}
void ok(int status,const char* message){need(status==PB11_STATUS_OK,message);}
std::vector<unsigned char> pack(fusion_source_state_v1* p){size_t n=0,written=0;ok(fusion_c_source_state_pack_size(p,&n),"pack size");std::vector<unsigned char>b(n);ok(fusion_c_source_state_pack(p,b.data(),b.size(),&written),"pack");need(written==n,"pack length");return b;}
struct Final {std::array<double,6> N;double Ue,Ui;std::vector<unsigned char> bytes;};
Final evolve(bool restart){
 constexpr int n=300,steps=8;constexpr double keV=1.602176634e-16,MeV=1.602176634e-13,ne=1e22,dt=.001/steps;
 std::vector<double> edges(n+1),s(6*n),t(6*n),sn(6*n),tn(6*n),zero(6*n),logs(48,15.);
 const int low=n/3;
 for(int j=1;j<=low;++j)edges[j]=1e-10*keV*std::pow(1e12,double(j-1)/(low-1));
 for(int j=low+1;j<=n;++j)edges[j]=(100.+24900.*(j-low)/(n-low))*keV;
 fusion_inert_ion_v1 carbon{.001*ne,12*1.66053906892e-27,36};
 std::array<double,6>N{0,(ne-6*carbon.density_m3)/2,(ne-6*carbon.density_m3)/2,0,0,0};
 double Ue=1.5*ne*5*keV,Ui=1.5*(N[1]+N[2]+carbon.density_m3)*20*keV;
 const double initialUe=Ue,initialUi=Ui;
 double Z2[6]{1,1,1,4,4,25};fusion_coupled_thermal_options_v1 o{};
 o.birth={2.5*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,8,8};
 o.channels[3]=1;o.max_source_rate_error=o.max_source_debit_error=1e-5;
 o.handoff_enabled=1;o.handoff_max_L1=o.handoff_max_mean_error=.001;
 fusion_source_state_v1* raw=nullptr;constexpr uint64_t tag=0x1213141516;
 ok(fusion_c_source_state_create(n,edges.data(),s.data(),t.data(),0,tag,&raw),"create");Handle h(raw,fusion_c_source_state_destroy);
 fusion_source_ledger_v1 cumulative{};double inert[6]{},time=0;uint64_t epoch=0;
 for(int step=1;step<=steps;++step){
  ok(fusion_c_source_state_snapshot_inert(h.get(),s.data(),t.data(),&cumulative,inert,&time,&epoch),"snapshot accepted");
  need(epoch==uint64_t(step-1),"accepted epoch");
  std::array<double,6>Nnew{};fusion_coupled_thermal_v1 trial{};uint64_t ticket=0;
  ok(fusion_c_source_state_begin(h.get(),dt,&ticket),"begin");
  ok(fusion_c_coupled_thermal_trial(dt,&o,n,edges.data(),N.data(),Ue,Ui,ne,Z2,1,&carbon,logs.data(),s.data(),t.data(),zero.data(),zero.data(),Nnew.data(),sn.data(),tn.data(),&trial),"coupled trial");
  if(step==1){
   need(trial.inert_ion_heat_J_m3[4]>0,"actual alpha carbon heating");
   need(fusion_c_source_state_stage(h.get(),ticket,sn.data(),tn.data(),&trial.ledger)!=0,"missing inert account must reject");
   need(fusion_c_source_state_commit(h.get(),ticket)!=0,"cannot commit omitted heat");
  }
  ok(fusion_c_source_state_stage_inert(h.get(),ticket,sn.data(),tn.data(),&trial.ledger,trial.inert_ion_heat_J_m3),"complete stage");
  ok(fusion_c_source_state_stage_inert(h.get(),ticket,sn.data(),tn.data(),&trial.ledger,trial.inert_ion_heat_J_m3),"replacement stage");
  // No thermal update until the kinetic candidate is known valid and commits.
  ok(fusion_c_source_state_commit(h.get(),ticket),"commit");N=Nnew;Ue=trial.electron_energy_J_m3;Ui=trial.ion_energy_J_m3;
  ok(fusion_c_source_state_snapshot_inert(h.get(),s.data(),t.data(),&cumulative,inert,&time,&epoch),"snapshot new");
  need(epoch==uint64_t(step)&&time>0,"new epoch");
  // Independently recompute complete local energy from persisted accounts.
  long double fast=0,Q=0;
  for(int i=0;i<6;++i)for(int j=0;j<n;++j)fast+=(static_cast<long double>(s[i*n+j])+t[i*n+j])*(static_cast<long double>(edges[j])+edges[j+1])/2;
  for(int ch=0;ch<5;++ch){fusion_nuclear_channel_v1 c{};ok(fusion_c_nuclear_channel(ch,&c),"channel");Q+=static_cast<long double>(cumulative.events_m3[ch])*c.q_J;}
  long double residual=static_cast<long double>(Ue)+Ui+fast+cumulative.neutron_energy_J_m3-initialUe-initialUi-Q;
  need(std::abs(residual)<1e-12L*(initialUe+initialUi+Q),"persisted total energy");
  need(cumulative.events_m3[3]>0&&inert[4]>0,"persisted real reactions and carbon heat");
  if(restart&&step==4){
   // This test's host bundle captures thermal state alongside the same epoch.
   const auto thermalN=N;const double thermalUe=Ue,thermalUi=Ui,savedTime=time;const uint64_t savedEpoch=epoch;
   const auto bytes=pack(h.get());fusion_source_state_v1* restored=nullptr;
   ok(fusion_c_source_state_unpack(bytes.data(),bytes.size(),tag,&restored),"restore");h.reset(restored);
   N=thermalN;Ue=thermalUe;Ui=thermalUi;
   ok(fusion_c_source_state_snapshot_inert(h.get(),s.data(),t.data(),&cumulative,inert,&time,&epoch),"restored snapshot");
   need(time==savedTime&&epoch==savedEpoch&&pack(h.get())==bytes,"same epoch restart bytes");
  }
 }
 std::cout<<"restart="<<restart<<" DT events="<<cumulative.events_m3[3]<<" carbon heat="<<inert[4]<<" Ue="<<Ue<<" Ui="<<Ui<<'\n';
 return {N,Ue,Ui,pack(h.get())};
}
}
int main(){try{const auto a=evolve(false),b=evolve(true);need(a.N==b.N&&a.Ue==b.Ue&&a.Ui==b.Ui&&a.bytes==b.bytes,"restart trajectory parity");std::cout<<"Coupled thermal inert-state restart passed\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
