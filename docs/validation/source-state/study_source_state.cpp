#include "fusion_source_state.h"
#include "fusion_coulomb.h"
#include "fusion_kinetics.h"
#include "fusion_nuclear_data.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <algorithm>
constexpr double ev=1.602176634e-19;
void ok(int rc){if(rc)throw std::runtime_error("API status "+std::to_string(rc));}
int main(){try{
 const int n=80,steps=400;const double dt=.001,birthrate=1e13,escape_rate=.3;
 std::vector<double> edges(n+1),zero(6*n),s(6*n),t(6*n),birth(n),escape(n,escape_rate),diff(2*(n-1)),next(6*n);
 for(int i=1;i<=n;++i)edges[i]=ev*std::exp((i-1)*std::log(4e6)/(n-1));
 int birthbin=0;double distance=1e99;
 for(int i=0;i<n;++i){double d=std::abs((edges[i]+edges[i+1])/2-3e6*ev);if(d<distance){distance=d;birthbin=i;}}
 birth[birthbin]=birthrate;
 fusion_source_state_v1 *reference=nullptr,*restarted=nullptr;
 ok(fusion_c_source_state_create(n,edges.data(),zero.data(),zero.data(),0,54321,&reference));
 ok(fusion_c_source_state_create(n,edges.data(),zero.data(),zero.data(),0,54321,&restarted));
 fusion_nuclear_mass_v1 alpha{},electron{},proton{};
 ok(fusion_c_nuclear_mass(4,&alpha));ok(fusion_c_nuclear_mass(7,&electron));ok(fusion_c_nuclear_mass(0,&proton));
 std::cout<<std::setprecision(17)<<"time_s,fast_N_m3,fast_U_J_m3,born_N_m3,born_U_J_m3,escaped_N_m3,escaped_U_J_m3,heat_e_J_m3,heat_p_J_m3,particle_relative_residual,energy_relative_residual\n";
 double max_n=0,max_u=0;
 for(int step=1;step<=steps;++step){double tm=(step-.5)*dt;
  double kt[2]={(1+2*tm/(steps*dt))*1e3*ev,(.5+tm/(steps*dt))*1e3*ev};
  fusion_maxwellian_bath_v1 baths[2]={{5.9683e19,electron.mass_kg,1,kt[0],15},{5.9683e19,proton.mass_kg,1,kt[1],15}};
  for(int b=0;b<2;++b)for(int i=0;i<n-1;++i){fusion_coulomb_energy_v1 coefficient{};
   ok(fusion_c_coulomb_energy(edges[i+1],alpha.mass_kg,2,&baths[b],&coefficient));diff[b*(n-1)+i]=coefficient.diffusion_J2_s;}
  for(auto state:{reference,restarted}){
   fusion_source_ledger_v1 cumulative{};double time;uint64_t epoch,ticket;
   ok(fusion_c_source_state_snapshot(state,s.data(),t.data(),&cumulative,&time,&epoch));
   if(step%17==0){ok(fusion_c_source_state_begin(state,2*dt,&ticket));auto invalid=s;invalid[0]=-1;
    fusion_source_ledger_v1 empty{};
    if(fusion_c_source_state_stage(state,ticket,invalid.data(),t.data(),&empty)==0)throw std::runtime_error("bad trial accepted");
    if(fusion_c_source_state_commit(state,ticket)==0)throw std::runtime_error("bad commit accepted");
    ok(fusion_c_source_state_discard(state,ticket));}
   double heat[2];fusion_kinetic_ledger_v1 fp{};next=zero;
   ok(fusion_c_energy_fp_trial(n,2,dt,edges.data(),s.data()+4*n,kt,diff.data(),birth.data(),escape.data(),0,next.data()+4*n,heat,&fp));
   fusion_source_ledger_v1 ledger{};
   ledger.external_born_number_m3[4]=fp.born_number_m3;ledger.external_born_energy_J_m3[4]=fp.born_energy_J_m3;
   ledger.escaped_number_m3[4]=fp.escaped_number_m3;ledger.escaped_energy_J_m3[4]=fp.escaped_energy_J_m3;
   ledger.heat_to_bath_J_m3[4*7]=heat[0];ledger.heat_to_bath_J_m3[4*7+1]=heat[1];
   ok(fusion_c_source_state_begin(state,dt,&ticket));
   ok(fusion_c_source_state_stage(state,ticket,next.data(),t.data(),&ledger));
   if(step%5==0)ok(fusion_c_source_state_stage(state,ticket,next.data(),t.data(),&ledger));
   ok(fusion_c_source_state_commit(state,ticket));
  }
  if(step==200){size_t length,written;ok(fusion_c_source_state_pack_size(restarted,&length));std::vector<unsigned char> packed(length);
   ok(fusion_c_source_state_pack(restarted,packed.data(),length,&written));fusion_c_source_state_destroy(restarted);restarted=nullptr;
   ok(fusion_c_source_state_unpack(packed.data(),written,54321,&restarted));}
  fusion_source_ledger_v1 ledger{};double time;uint64_t epoch;
  ok(fusion_c_source_state_snapshot(restarted,s.data(),t.data(),&ledger,&time,&epoch));
  long double number=0,energy=0;for(int i=0;i<n;++i){number+=s[4*n+i];energy+=static_cast<long double>(s[4*n+i])*(edges[i]+edges[i+1])/2;}
  double nr=static_cast<double>(std::abs(number+ledger.escaped_number_m3[4]-ledger.external_born_number_m3[4])/ledger.external_born_number_m3[4]);
  double ur=static_cast<double>(std::abs(energy+ledger.escaped_energy_J_m3[4]+ledger.heat_to_bath_J_m3[28]+ledger.heat_to_bath_J_m3[29]-ledger.external_born_energy_J_m3[4])/ledger.external_born_energy_J_m3[4]);
  max_n=std::max(max_n,nr);max_u=std::max(max_u,ur);
  std::cout<<time<<','<<double(number)<<','<<double(energy)<<','<<ledger.external_born_number_m3[4]<<','<<ledger.external_born_energy_J_m3[4]<<','<<ledger.escaped_number_m3[4]<<','<<ledger.escaped_energy_J_m3[4]<<','<<ledger.heat_to_bath_J_m3[28]<<','<<ledger.heat_to_bath_J_m3[29]<<','<<nr<<','<<ur<<'\n';
  long double exact=birthrate/escape_rate*(-std::expm1(-step*std::log1p(escape_rate*dt)));
  if(std::abs(number-exact)>1e-10L*exact)throw std::runtime_error("independent BE number reference");
 }
 size_t size,written;ok(fusion_c_source_state_pack_size(reference,&size));std::vector<unsigned char> a(size),b(size);
 ok(fusion_c_source_state_pack(reference,a.data(),size,&written));ok(fusion_c_source_state_pack(restarted,b.data(),size,&written));
 if(a!=b)throw std::runtime_error("restart diverged");
 fusion_c_source_state_destroy(reference);fusion_c_source_state_destroy(restarted);
 std::cerr<<std::setprecision(17)<<"400 steps completed; restart identical including ledgers/tickets; max relative N="<<max_n<<" U="<<max_u<<"; birth cell center MeV="<<(edges[birthbin]+edges[birthbin+1])/(2e6*ev)<<'\n';
 return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
