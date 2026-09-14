#include "fusion_coupled_thermal.h"
#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <limits>
using R=long double;
void need(bool p,const char*m){if(!p)throw std::runtime_error(m);}
void ok(int s,const char*m){need(s==0,m);}
int main(){try{
 constexpr int n=200;constexpr double keV=1.602176634e-16,MeV=1.602176634e-13,ne=1e22,dt=.000125;
 std::vector<double>edges(n+1),zero(6*n),logs(48,15.);int low=n/3;
 for(int j=1;j<=low;++j)edges[j]=1e-10*keV*std::pow(1e12,double(j-1)/(low-1));
 for(int j=low+1;j<=n;++j)edges[j]=(100.+24900.*(j-low)/(n-low))*keV;
 fusion_inert_ion_v1 carbon{.001*ne,12*1.66053906892e-27,36};
 double z2[6]{1,1,1,4,4,25};std::array<double,6>N{0,(ne-6*carbon.density_m3)/2,(ne-6*carbon.density_m3)/2,0,0,0};
 double Ue=1.5*ne*5*keV,Ui=double(1.5L*(R(carbon.density_m3)+N[1]+N[2])*20*keV);
 const double exactTi=double(R(Ui)/(1.5L*(R(carbon.density_m3)+N[1]+N[2])));
 fusion_coupled_thermal_options_v1 o{};o.birth={2.5*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,4,4};
 o.channels[3]=1;o.max_source_rate_error=o.max_source_debit_error=1e-5;o.handoff_enabled=1;o.handoff_max_L1=o.handoff_max_mean_error=.001;
 fusion_birth_table_control_v1 control{.001,.001,.001,.001,1e-5,1e-5,64,512,12};fusion_birth_table_v1*raw=nullptr;
 ok(fusion_c_birth_table_create(3,exactTi,23*keV,&o.birth,&control,n,edges.data(),&raw),"table build");
 std::unique_ptr<fusion_birth_table_v1,decltype(&fusion_c_birth_table_destroy)>owner(raw,fusion_c_birth_table_destroy);
 std::array<const fusion_birth_table_v1*,5>tables{};tables[3]=raw;
 struct State{std::array<double,6>N{};double Ue=0,Ui=0;std::vector<double>s,t;R carbon=0;};
 State direct{N,Ue,Ui,zero,zero,0},table=direct;
 auto trial=[&](State&v,bool tabulated,const double*edge,const fusion_coupled_thermal_options_v1&options){
  std::vector<double>snew(6*n,99),tnew(6*n,99);std::array<double,6>nextN{};nextN.fill(99);fusion_coupled_thermal_v1 result{};
  int status=tabulated?fusion_c_coupled_thermal_table_trial(dt,&options,tables.data(),n,edge,v.N.data(),v.Ue,v.Ui,ne,z2,1,&carbon,logs.data(),v.s.data(),v.t.data(),zero.data(),zero.data(),nextN.data(),snew.data(),tnew.data(),&result):fusion_c_coupled_thermal_trial(dt,&options,n,edge,v.N.data(),v.Ue,v.Ui,ne,z2,1,&carbon,logs.data(),v.s.data(),v.t.data(),zero.data(),zero.data(),nextN.data(),snew.data(),tnew.data(),&result);
  if(status){need(std::all_of(snew.begin(),snew.end(),[](double x){return x==0;})&&std::all_of(tnew.begin(),tnew.end(),[](double x){return x==0;}),"clear failed kinetic output");for(double x:nextN)need(x==0,"clear failed thermal output");need(result.electron_energy_J_m3==0&&result.ion_energy_J_m3==0,"clear result");return status;}
  fusion_thermal_increment_v1 increment{};
  ok(fusion_c_coupled_thermal_increment(&result.ledger,result.inert_ion_heat_J_m3,&increment),"actual coupled thermal increments");
  auto reconstruct=[](double old,double amount,double next){return std::abs(R(old)+amount-next)<=16*std::numeric_limits<double>::epsilon()*(std::abs(R(old))+std::abs(R(next)));};
  for(int i=0;i<6;++i)need(reconstruct(v.N[i],increment.thermal_number_m3[i],nextN[i]),"actual particle reconstruction");
  need(reconstruct(v.Ue,increment.electron_energy_J_m3,result.electron_energy_J_m3)&&reconstruct(v.Ui,increment.ion_energy_J_m3,result.ion_energy_J_m3),"actual energy reconstruction including carbon");
  v.N=nextN;v.Ue=result.electron_energy_J_m3;v.Ui=result.ion_energy_J_m3;v.s.swap(snew);v.t.swap(tnew);for(double x:result.inert_ion_heat_J_m3)v.carbon+=x;return status;
 };
 ok(trial(direct,false,edges.data(),o),"direct first");ok(trial(table,true,edges.data(),o),"table first");
 need(direct.N==table.N&&direct.s==table.s&&direct.t==table.t&&direct.Ue==table.Ue&&direct.Ui==table.Ui&&direct.carbon==table.carbon,"exact endpoint physical trajectory parity");
 // Reaction-conditioned fuel debit initially cools the ion pool below20keV.
 // A table starting at the initial Ti must reject, not clamp that real state.
 need(trial(table,true,edges.data(),o)==PB11_STATUS_OUT_OF_RANGE,"real early cooling exits endpoint table");
 fusion_birth_table_v1*wide=nullptr;
 ok(fusion_c_birth_table_create(3,19.5*keV,23*keV,&o.birth,&control,n,edges.data(),&wide),"broader trajectory table");
 owner.reset(wide);raw=wide;tables[3]=raw;
 for(int step=1;step<8;++step){ok(trial(direct,false,edges.data(),o),"direct evolving");ok(trial(table,true,edges.data(),o),"table evolving");}
 auto close=[](R a,R b){return std::abs(a-b)<=.002L*std::max(std::abs(a),std::abs(b));};
 need(close(direct.Ue,table.Ue)&&close(direct.Ui,table.Ui)&&close(direct.carbon,table.carbon),"evolving bath/weak heat sensitivity");
 R sd=0,st=0,Ed=0,Et=0;for(int j=0;j<6*n;++j){sd+=direct.s[j]+direct.t[j];st+=table.s[j]+table.t[j];R center=(R(edges[j%n])+edges[j%n+1])/2;Ed+=(R(direct.s[j])+direct.t[j])*center;Et+=(R(table.s[j])+table.t[j])*center;}
 need(close(sd,st)&&close(Ed,Et),"evolving fast sensitivity");
 const State saved=table;auto bad=o;bad.birth.cm_order=8;need(trial(table,true,edges.data(),bad)!=0,"model mismatch rejected");
 auto changed=edges;changed[low]*=1.00001;need(trial(table,true,changed.data(),o)!=0,"grid mismatch rejected");
 tables[3]=nullptr;need(trial(table,true,edges.data(),o)!=0,"missing active channel rejected");tables[3]=raw;
 need(table.N==saved.N&&table.s==saved.s&&table.t==saved.t&&table.Ue==saved.Ue&&table.Ui==saved.Ui,"failed table trials preserve accepted state");
 State outside=table;outside.Ui*=2;need(trial(outside,true,edges.data(),o)==PB11_STATUS_OUT_OF_RANGE,"outside table rejects without fallback");
 std::cout<<"Coupled table/direct parity and failure tests passed; carbon relative difference="<<double((table.carbon-direct.carbon)/direct.carbon)<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
