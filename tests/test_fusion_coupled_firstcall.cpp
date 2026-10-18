#include "fusion_coupled_thermal.h"
#include <array>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>
namespace {
constexpr int cells=64,workers=8;
constexpr double keV=1.602176634e-16,MeV=1.602176634e-13;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
template<class T>bool equal_bytes(const T& a,const T& b){return a.size()==b.size()&&std::memcmp(a.data(),b.data(),a.size()*sizeof(a[0]))==0;}
struct Input {
 std::array<double,cells+1> edges{};
 std::array<double,6> number{},charge2{{1,1,1,4,4,25}};
 std::array<double,6*cells> s{},t{},birth{},escape{};
 std::array<double,48> logs{};
 fusion_inert_ion_v1 carbon{7e16,12*1.66053906892e-27,36};
 double ne=7e19,ue=0,ui=0;
 fusion_coupled_thermal_options_v1 options{};
 Input(){
  edges[0]=0;for(int j=1;j<=cells;++j)edges[j]=.01*keV*std::pow(2e6,double(j-1)/(cells-1));
  number[4]=1e-5*ne;number[1]=number[3]=(ne-6*carbon.density_m3-2*number[4])/3;
  double pool=carbon.density_m3;for(double x:number)pool+=x;
  ue=1.5*ne*.5*keV;ui=1.5*pool*.1*keV;
  options.birth={10*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,8,8};
  options.max_source_rate_error=options.max_source_debit_error=1e-5;
  options.handoff_max_L1=options.handoff_max_mean_error=.001;options.handoff_enabled=1;
  logs.fill(10.);
  for(int species:{0,4}){
   double desired=(species==0?14700:3600)*keV;int j=0;while(j+1<cells&&(edges[j]+edges[j+1])/2<desired)++j;
   s[species*cells+j]=1e11;t[species*cells+j]=2e10;
  }
 }
};
struct Output {
 std::array<double,6> number{};
 std::array<double,6*cells> s{},t{};
 fusion_coupled_thermal_v1 result{};
 int status=-1;
};
void trial(const Input& a,Output& b){
 b.status=fusion_c_coupled_thermal_trial(1.5625e-6,&a.options,cells,a.edges.data(),a.number.data(),a.ue,a.ui,a.ne,a.charge2.data(),1,&a.carbon,a.logs.data(),a.s.data(),a.t.data(),a.birth.data(),a.escape.data(),b.number.data(),b.s.data(),b.t.data(),&b.result);
}
// Copy defined fields only: no comparison of structure tail padding.
std::vector<unsigned char> input_bytes(const Input& a){
 std::vector<unsigned char> bytes;
 auto add=[&](const void* p,size_t n){auto q=static_cast<const unsigned char*>(p);bytes.insert(bytes.end(),q,q+n);};
 for(auto p:{&a.s,&a.t,&a.birth,&a.escape})add(p->data(),sizeof(*p));
 add(a.edges.data(),sizeof(a.edges));add(a.number.data(),sizeof(a.number));add(a.charge2.data(),sizeof(a.charge2));add(a.logs.data(),sizeof(a.logs));
 for(double x:{a.carbon.density_m3,a.carbon.mass_kg,a.carbon.mean_charge_squared,a.ne,a.ue,a.ui})add(&x,sizeof(x));
 // Options are snapshotted as their existing object representation; this object
 // is never assigned after snapshot, so even padding must stay unchanged.
 add(&a.options,sizeof(a.options));return bytes;
}
}
int main(){try{
 const Input input;const auto original=input_bytes(input);
 std::array<Output,workers> outputs;
 std::mutex mutex;std::condition_variable ready;
 int arrived=0;bool release=false;
 std::vector<std::thread> threads;
 // No coupled API call occurs before every worker reaches this gate.
 for(int i=0;i<workers;++i)threads.emplace_back([&,i](){
  {std::unique_lock<std::mutex> lock(mutex);++arrived;ready.notify_all();ready.wait(lock,[&](){return release;});}
  trial(input,outputs[i]);
 });
 {std::unique_lock<std::mutex> lock(mutex);ready.wait(lock,[&](){return arrived==workers;});release=true;}ready.notify_all();
 for(auto& thread:threads)thread.join();
 for(const auto& output:outputs)require(output.status==0,"concurrent trial status");
 require(input_bytes(input)==original,"concurrent calls mutated input");
 Output reference;trial(input,reference);require(reference.status==0,"serial trial status");
 require(input_bytes(input)==original,"serial call mutated input");
 double electron_heat=0,ion_heat=0;
 for(int i=0;i<6;++i){electron_heat+=reference.result.ledger.heat_to_bath_J_m3[7*i];for(int b=1;b<7;++b)ion_heat+=reference.result.ledger.heat_to_bath_J_m3[7*i+b];}
 require(std::isfinite(electron_heat)&&electron_heat!=0,"electron collision heat absent");
 require(std::isfinite(ion_heat)&&ion_heat!=0,"ion collision heat absent");
 require(!equal_bytes(reference.s,input.s)||!equal_bytes(reference.t,input.t),"FP populations did not evolve");
 for(const auto& output:outputs){
  require(equal_bytes(output.number,reference.number),"thermal output differs");
  require(equal_bytes(output.s,reference.s)&&equal_bytes(output.t,reference.t),"kinetic output differs");
  require(std::memcmp(&output.result,&reference.result,offsetof(fusion_coupled_thermal_v1,handoff_projected))==0,"result double prefix differs");
  require(std::memcmp(output.result.handoff_projected,reference.result.handoff_projected,6*sizeof(int))==0,"handoff output differs");
 }
 std::printf("FIRSTCALL_PASS workers=%d cells=%d serial_calls_after_join=1 electron_heat=%.17g ion_heat=%.17g input_unchanged=1 bitwise_equal=1\n",workers,cells,electron_heat,ion_heat);
 return 0;
}catch(const std::exception& error){std::fprintf(stderr,"FIRSTCALL_FAIL %s\n",error.what());return 1;}}
