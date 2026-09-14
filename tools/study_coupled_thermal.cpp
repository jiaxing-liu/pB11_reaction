#include "fusion_coupled_thermal.h"
#include "fusion_nuclear_data.h"
#include <array>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <chrono>
using R=long double;
constexpr double keV=1.602176634e-16,MeV=1.602176634e-13;
void require(bool p,const char*s){if(!p)throw std::runtime_error(s);}
int main(int argc,char**argv){try{
 int steps=argc>1?std::stoi(argv[1]):64,n=argc>2?std::stoi(argv[2]):800;
 std::string fuel=argc>3?argv[3]:"dt";double duration=argc>4?std::stod(argv[4]):.01;
 int handoff=argc>5?std::stoi(argv[5]):1;
 int alpha_order=argc>6?std::stoi(argv[6]):8;
 require(steps>0&&n>=100&&duration>0&&(handoff==0||handoff==1),"arguments");
 double ne=1e22,N[6]{},Z2[6]={1,1,1,4,4,25};
 fusion_inert_ion_v1 carbon{.001*ne,12*1.66053906892e-27,36};
 double fuel_ne=ne-6*carbon.density_m3,Ti=20*keV,Te=5*keV;
 fusion_coupled_thermal_options_v1 o{};
 o.birth={2.5*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,8,8};
 o.birth.nq=o.birth.ncos=alpha_order;
 o.max_source_rate_error=1e-5;o.max_source_debit_error=1e-5;
 o.handoff_max_L1=.001;o.handoff_max_mean_error=.001;o.handoff_enabled=handoff;
 if(fuel=="dt"){N[1]=N[2]=fuel_ne/2;o.channels[3]=1;}
 else if(fuel=="dd"){N[1]=fuel_ne;o.channels[1]=o.channels[2]=o.channels[3]=o.channels[4]=1;}
 else if(fuel=="dhe3"){N[1]=N[3]=fuel_ne/3;o.channels[4]=1;}
 else if(fuel=="pb"){N[0]=N[5]=fuel_ne/6;Ti=100*keV;o.channels[0]=1;}
 else throw std::runtime_error("fuel");
 R Nt=carbon.density_m3;for(double x:N)Nt+=x;
 double Ue=1.5*ne*Te,Ui=double(1.5L*Nt*Ti),initialUe=Ue,initialUi=Ui;
 std::vector<double> edges(n+1),s(6*n),t(6*n),sn(6*n),tn(6*n),external(6*n),escape(6*n),logs(6*8,15.);
 int low=n/3;edges[0]=0;
 const double first_keV=fuel=="pb"?1e-18:1e-10;
 for(int j=1;j<=low;++j)edges[j]=first_keV*keV*std::pow(100/first_keV,double(j-1)/(low-1));
 for(int j=low+1;j<=n;++j)edges[j]=(100.+24900.*(j-low)/(n-low))*keV;
 const double table_tol=argc>7?std::stod(argv[7]):0;
 require(std::isfinite(table_tol)&&table_tol>=0&&table_tol<=1,"table tolerance");
 std::array<std::shared_ptr<fusion_birth_table_v1>,5> table_owner;
 std::array<const fusion_birth_table_v1*,5> tables{};
 if(table_tol>0){
  require(argc>9,"table mode requires lower and upper keV arguments");
  double lower=std::stod(argv[8])*keV,upper=std::stod(argv[9])*keV;
  fusion_birth_table_control_v1 control{table_tol,table_tol,table_tol,table_tol,o.max_source_rate_error,o.max_source_debit_error,128,2048,16};
  for(int ch=0;ch<5;++ch)if(o.channels[ch]){
   auto source=o.birth;if(ch!=0)source.pb_low=0;
   fusion_birth_table_v1*raw=nullptr;auto start=std::chrono::steady_clock::now();
   int status=fusion_c_birth_table_create(ch,lower,upper,&source,&control,n,edges.data(),&raw);
   if(status){std::fprintf(stderr,"table failed channel=%d status=%d\n",ch,status);return status;}
   table_owner[ch]=std::shared_ptr<fusion_birth_table_v1>(raw,fusion_c_birth_table_destroy);tables[ch]=raw;
   fusion_birth_table_info_v1 info{};require(fusion_c_birth_table_info(raw,&info)==0,"table info");
   std::fprintf(stderr,"table channel=%d knots=%d direct_evaluations=%d build_seconds=%.9g max_sampled_rate=%.9g max_sampled_debit=%.9g max_sampled_number_L1=%.9g max_sampled_energy_L1=%.9g\n",ch,info.knots,info.direct_evaluations,std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(),info.max_validated_rate_error,info.max_validated_debit_error,info.max_validated_number_L1,info.max_validated_energy_L1);
  }
 }
 R Q=0,neutronE=0,neutronN=0,maxResidual=0,inertHeat=0,heatE=0,heatI=0,removedU=0;
 std::array<R,5> events{};std::array<R,6> ash{};int projections=0;
 std::puts("step,time_s,Te_keV,Ti_keV,Ue_J_m3,Ui_J_m3,fast_energy_J_m3,neutron_energy_J_m3,Q_J_m3,energy_relative_residual,inert_heat_J_m3,events_pb,events_dd_tp,events_dd_he3n,events_dt,events_dhe3,thermal_H,thermal_D,thermal_T,thermal_He3,thermal_He4,thermal_B11,fast_He4,thermalized_He4,projections,electron_heat_J_m3,network_ion_heat_J_m3,removed_thermal_energy_J_m3");
 for(int step=1;step<=steps;++step){double Nnew[6]{};fusion_coupled_thermal_v1 result{};
  int st=table_tol>0?
   fusion_c_coupled_thermal_table_trial(duration/steps,&o,tables.data(),n,edges.data(),N,Ue,Ui,ne,Z2,1,&carbon,logs.data(),s.data(),t.data(),external.data(),escape.data(),Nnew,sn.data(),tn.data(),&result):
   fusion_c_coupled_thermal_trial(duration/steps,&o,n,edges.data(),N,Ue,Ui,ne,Z2,1,&carbon,logs.data(),s.data(),t.data(),external.data(),escape.data(),Nnew,sn.data(),tn.data(),&result);
  if(st){std::fprintf(stderr,"trial failed fuel=%s step=%d status=%d Ti=%.9gkeV Te=%.9gkeV\n",fuel.c_str(),step,st,Ui/(1.5*double(Nt))/keV,Ue/(1.5*ne)/keV);return st;}
  // Atomic acceptance in this standalone driver only, after successful trial.
  std::copy(Nnew,Nnew+6,N);Ue=result.electron_energy_J_m3;Ui=result.ion_energy_J_m3;s.swap(sn);t.swap(tn);
  Nt=carbon.density_m3;for(double x:N)Nt+=x;
  for(int ch=0;ch<5;++ch){fusion_nuclear_channel_v1 r{};fusion_c_nuclear_channel(ch,&r);events[ch]+=result.ledger.events_m3[ch];Q+=R(result.ledger.events_m3[ch])*r.q_J;}
  neutronE+=result.ledger.neutron_energy_J_m3;neutronN+=result.ledger.neutron_number_m3;
  R fastU=0,fastHe=0;for(int i=0;i<6;++i){inertHeat+=result.inert_ion_heat_J_m3[i];heatE+=result.ledger.heat_to_bath_J_m3[i*7];
   for(int b=1;b<7;++b)heatI+=result.ledger.heat_to_bath_J_m3[i*7+b];
   removedU+=result.ledger.thermal_consumed_energy_J_m3[i];ash[i]+=result.ledger.handed_off_number_m3[i];projections+=result.handoff_projected[i];for(int j=0;j<n;++j){R pop=R(s[i*n+j])+t[i*n+j];fastU+=pop*(R(edges[j])+edges[j+1])/2;if(i==4)fastHe+=pop;}}
  R residual=(R(Ue)+Ui+fastU+neutronE-initialUe-initialUi-Q)/(R(initialUe)+initialUi+Q);maxResidual=std::max(maxResidual,std::abs(residual));
  std::printf("%d,%.12g,%.12g,%.12Lg,%.14g,%.14g,%.14Lg,%.14Lg,%.14Lg,%.12Lg,%.14Lg",step,step*duration/steps,Ue/(1.5*ne)/keV,R(Ui)/(1.5L*Nt)/keV,Ue,Ui,fastU,neutronE,Q,residual,inertHeat);
  for(R x:events)std::printf(",%.14Lg",x);
  for(double x:N)std::printf(",%.14g",x);
  std::printf(",%.14Lg,%.14Lg,%d,%.14Lg,%.14Lg,%.14Lg\n",fastHe,ash[4],projections,heatE,heatI,removedU);
 }
 std::fprintf(stderr,"fuel=%s steps=%d cells=%d duration=%.12g max_total_energy_relative_residual=%.12Lg projections=%d neutron_number=%.12Lg\n",fuel.c_str(),steps,n,duration,maxResidual,projections,neutronN);
 require(maxResidual<1e-9,"cumulative energy budget");return 0;
}catch(const std::exception&e){std::fprintf(stderr,"ERROR %s\n",e.what());return 1;}}
