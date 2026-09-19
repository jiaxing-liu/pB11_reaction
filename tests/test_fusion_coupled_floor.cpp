#define main original_fast_tests_main
#include "test_fusion_coupled_fast.cpp"
#undef main
#include <cstring>
#include <iomanip>
int floor_call(const Inputs& a,double dt,Trial& b,fusion_coupled_floor_limits_v1* limits,fusion_birth_floor_ledger_v1& ledger){
 fusion_handoff_diagnostics_v1 d{};fusion_beam_table_usage_v1 u{};
 return fusion_c_coupled_sources_floor_trial(dt,&a.options,&a.fast_options,nullptr,0,nullptr,0,a.grid.cells(),a.grid.edges.data(),a.thermal_number.data(),a.electron_energy_J_m3,a.ion_energy_J_m3,a.electron_density_m3,a.thermal_charge_squared.data(),0,nullptr,a.coulomb_logs.data(),a.old_s.data(),a.old_t.data(),a.external_birth.data(),a.escape.data(),b.thermal_number.data(),b.s.data(),b.t.data(),&b.result,&d,&u,limits,&ledger);
}
std::vector<double> captured_packets;
fusion_birth_packets_v1 captured_meta{};
int floor_packet_call(const Inputs&a,double dt,Trial&b,fusion_coupled_floor_limits_v1*limits,fusion_birth_floor_ledger_v1&ledger){
 fusion_handoff_diagnostics_v1 d{};fusion_beam_table_usage_v1 u{};uint64_t outside=0;
 captured_packets.assign(70*a.grid.cells(),-7);
 return fusion_c_coupled_sources_packets_trial(dt,&a.options,&a.fast_options,nullptr,0,nullptr,0,a.grid.cells(),a.grid.edges.data(),a.thermal_number.data(),a.electron_energy_J_m3,a.ion_energy_J_m3,a.electron_density_m3,a.thermal_charge_squared.data(),0,nullptr,a.coulomb_logs.data(),a.old_s.data(),a.old_t.data(),a.external_birth.data(),a.escape.data(),b.thermal_number.data(),b.s.data(),b.t.data(),&b.result,&d,&u,limits,&ledger,FUSION_BEAM_TABLE_STRICT,&outside,captured_packets.data(),&captured_meta);
}
int main(){try{
 Inputs a(true);a.fast_options=make_fast_options({{0,0,0,1,0}});
 Trial b=make_trial(a),c=make_trial(a);fusion_coupled_floor_limits_v1 limits{.02,.001};fusion_birth_floor_ledger_v1 ledger{};
 int old=call_fast(a,1e-6,b),now=floor_call(a,1e-6,c,&limits,ledger);
 std::cout<<"DT legacy="<<old<<" floor="<<now<<std::endl;
 require(old==0&&now==0,"DT trials");
 require(b.s==c.s&&b.t==c.t&&b.thermal_number==c.thermal_number&&std::memcmp(&b.result,&c.result,sizeof(b.result))==0,"zero-spill parity");
 for(double x:ledger.born_number_m3)require(x==0,"DT no spill");
 auto invalid=limits;invalid.max_ion_energy_fraction=-1;
 require(floor_call(a,1e-6,c,&invalid,ledger)==PB11_STATUS_OUT_OF_RANGE,"invalid gate");
 for(double x:c.s)require(x==0,"failure clears s");for(double x:ledger.born_number_m3)require(x==0,"failure clears floor");
 for(int j=1;j<=a.grid.cells();++j)a.grid.edges[j]=2*kKeVJ*std::pow(12500.,double(j-1)/(a.grid.cells()-1));
 a.grid.edges.back()=25*kMeVJ;
 a.thermal_number.fill(0);a.thermal_number[5]=1e19;a.ion_energy_J_m3=1.5e19*10*kKeVJ;
 a.options=make_options({{0,0,0,0,0}});a.fast_options=make_fast_options({{1,0,0,0,0}});
 std::fill(a.old_s.begin(),a.old_s.end(),0);std::fill(a.old_t.begin(),a.old_t.end(),0);
 a.old_s[nearest_cell(a.grid,600*kKeVJ)]=1e20;
 a.options.birth.relative_max_J=2500*kKeVJ;a.options.birth.cutoff_J=.001*1.602176634e-13;
 a.options.birth.relative_order=16;a.fast_options.angular_order=16;
 limits.max_center_over_ion_kT=.2;
 b=make_trial(a);c=make_trial(a);old=call_fast(a,1e-3,b);now=floor_packet_call(a,1e-3,c,&limits,ledger);
 std::cout<<std::setprecision(17)<<"pB legacy="<<old<<" floor="<<now<<" lowN="<<ledger.born_number_m3[4]<<" lowU="<<ledger.born_energy_J_m3[4]<<" correction="<<ledger.ion_energy_correction_J_m3[4]<<std::endl;
 require(old==PB11_STATUS_OUT_OF_RANGE,"legacy rejects real low pB");require(now==0,"floor pB passes");require(ledger.born_number_m3[4]>0,"real low source");
 require(captured_meta.events_m3[0][0]==0&&captured_meta.events_m3[1][0]>0,"pB fast-only event packet");
 require(close_scaled(captured_meta.below_number_m3[1][0][4],ledger.born_number_m3[4],1e-12)&&close_scaled(captured_meta.below_energy_J_m3[1][0][4],ledger.born_energy_J_m3[4],1e-12),"physical below packet matches original ledger");
 require(captured_meta.below_energy_J_m3[1][0][4]!=ledger.mapped_energy_J_m3[4],"packet does not replace physical below energy by floor energy");
 long double expected=a.ion_energy_J_m3,delta=0;
 for(int i=0;i<6;++i){expected-=c.result.ledger.thermal_consumed_energy_J_m3[i];
 expected+=c.result.inert_ion_heat_J_m3[i]+c.result.ledger.handed_off_energy_J_m3[i];
 for(int j=1;j<7;++j)expected+=c.result.ledger.heat_to_bath_J_m3[7*i+j];
 delta+=ledger.ion_energy_correction_J_m3[i];}
 std::cout<<"ion_before_numerical="<<double(expected)<<" ion_expected="<<double(expected+delta)<<" ion_actual="<<c.result.ion_energy_J_m3<<" delta="<<double(delta)<<std::endl;
 require(double(expected)!=double(expected+delta),"correction changes representable ion reservoir");
 require(std::abs((long double)c.result.ion_energy_J_m3-(expected+delta))<=2*std::numeric_limits<double>::epsilon()*std::abs(expected+delta),"separate ion ledger reconstruction");
 Trial accepted=c;
 limits.max_center_over_ion_kT=0;
 require(floor_packet_call(a,1e-3,c,&limits,ledger)==PB11_STATUS_OUT_OF_RANGE,"real source rejected by explicit zero displacement gate");
 for(double x:c.s)require(x==0,"late reject clears source");
 for(double x:captured_packets)require(x==0,"late reject clears packets");
 require(captured_meta.events_m3[1][0]==0&&captured_meta.below_number_m3[1][0][4]==0,"late reject clears packet metadata");
 for(double x:ledger.born_number_m3)require(x==0,"late reject clears ledger");
 limits.max_center_over_ion_kT=.2;
 require(floor_call(a,1e-3,c,&limits,ledger)==0,"retry succeeds");
 require(c.s==accepted.s&&c.t==accepted.t&&c.thermal_number==accepted.thermal_number&&std::memcmp(&c.result,&accepted.result,sizeof(c.result))==0,"rejected call leaves no hidden state");
 
 std::cout<<"PASS diagnostic coupled floor"<<std::endl;return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
