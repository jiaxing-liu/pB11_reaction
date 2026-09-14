#include "fusion_kinetic_geometry.h"
#include "fusion_radial_transport.h"
#include "fusion_energy_work.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
namespace {
void check(bool x,const char* m){if(!x)throw std::runtime_error(m);}
void ok(int status,const char* m){if(status)throw std::runtime_error(std::string(m)+" status="+std::to_string(status));}
void near(long double x,long double y,const char* m){check(std::abs(x-y)<=1e-11L*std::max(std::abs(x),std::abs(y)),m);}
struct Fixture {
 static constexpr int z=3,c=3;static constexpr size_t components=6*c,total=z*components;
 std::array<double,c+1> edges{0,1,3,7};std::array<double,z> va{2,3,4},vb{2.1,2.9,4.3},C{-.3,0,.2};
 std::array<double,z+1> a{0,-.1,.2,-.15};
 std::vector<double> k,bs,bt,s,t,ns,nt;std::array<fusion_transport_ledger_v1,z> ledgers{};
 double dt=.2;
 Fixture(double scale=1):k(components*(z+1)),bs(2*components),bt(2*components),s(total),t(total),ns(total),nt(total){
  for(double& e:edges)e*=scale;
  for(size_t j=0;j<components;++j){
   for(int face=1;face<=z;++face)k[j*(z+1)+face]=.01*(1+j%3);
   bs[2*j]=.25*(j+1);bs[2*j+1]=.2*(j+1);bt[2*j]=.15*(j+1);bt[2*j+1]=.1*(j+1);
   for(int zone=0;zone<z;++zone){s[zone*components+j]=.3*(zone+1)*(j+1);t[zone*components+j]=.13*(z-zone)*(j+2);}
  }
 }
 int run(){return fusion_c_kinetic_geometry_trial(z,c,dt,edges.data(),va.data(),vb.data(),a.data(),C.data(),
  k.data(),bs.data(),bt.data(),s.data(),t.data(),ns.data(),nt.data(),ledgers.data());}
};
void compare_sequential_and_state(Fixture& f){
 ok(f.run(),"geometry composition");
 // Independent single-component calls check all species/energy/S/T layouts.
 for(int part=0;part<2;++part){
  const auto& input=part==0?f.s:f.t;const auto& boundary=part==0?f.bs:f.bt;
  const auto& actual=part==0?f.ns:f.nt;std::vector<double> transported(Fixture::total);
  for(size_t j=0;j<Fixture::components;++j){double old[3],radial[3],face[4];fusion_radial_ledger_v1 r{};
   for(int z=0;z<3;++z)old[z]=input[z*Fixture::components+j];
   ok(fusion_c_radial_transport_trial(3,1,f.dt,f.va.data(),f.vb.data(),f.a.data(),f.k.data()+4*j,
    old,boundary.data()+2*j,radial,face,&r),"reference single component radial");
   for(int z=0;z<3;++z)transported[z*Fixture::components+j]=radial[z];
  }
  for(int z=0;z<3;++z)for(int sp=0;sp<6;++sp){double worked[3],face[4];fusion_energy_work_ledger_v1 r{};
   const size_t offset=z*Fixture::components+3*sp;
   ok(fusion_c_energy_work_trial(3,f.dt,f.C[z],f.edges.data(),transported.data()+offset,worked,face,&r),"reference energy work");
   for(int e=0;e<3;++e)near(actual[offset+e],worked[e],"all-component layout parity");
  }
 }
 // Independent moments from physical populations; validate output ledger signs/units.
 for(int z=0;z<3;++z)for(int sp=0;sp<6;++sp){
  long double n0=0,u0=0,n1=0,u1=0;
  for(int e=0;e<3;++e){size_t i=z*Fixture::components+3*sp+e;const long double energy=(static_cast<long double>(f.edges[e])+f.edges[e+1])/2;
   n0+=(static_cast<long double>(f.s[i])+f.t[i])*f.va[z];u0+=(static_cast<long double>(f.s[i])+f.t[i])*f.va[z]*energy;
   n1+=(static_cast<long double>(f.ns[i])+f.nt[i])*f.vb[z];u1+=(static_cast<long double>(f.ns[i])+f.nt[i])*f.vb[z]*energy;
  }
  const auto& l=f.ledgers[z];near(n1+l.lower_number[sp]+l.upper_number[sp],n0+l.spatial_number[sp],"extensive particle balance");
  near(u1+l.lower_energy_J[sp]+l.upper_energy_J[sp],u0+l.spatial_energy_J[sp]+l.work_J[sp],"extensive energy balance");
 }
 fusion_source_state_v1* states[3]{};uint64_t tickets[3]{};fusion_source_ledger_v1 source{};double inert[6]{};
 for(int z=0;z<3;++z){size_t offset=z*Fixture::components;
  ok(fusion_c_source_state_create_volume(3,f.edges.data(),f.s.data()+offset,f.t.data()+offset,f.va[z],0,17,&states[z]),"create volume context");
  ok(fusion_c_source_state_begin(states[z],f.dt,&tickets[z]),"begin geometry state");
  ok(fusion_c_source_state_stage_volume(states[z],tickets[z],f.ns.data()+offset,f.nt.data()+offset,f.vb[z],f.va[z],
   &source,inert,&f.ledgers[z]),"stage direct composite ledger");
 }
 ok(fusion_c_source_state_commit_many(3,states,tickets),"batch composed state");
 for(auto* state:states)fusion_c_source_state_destroy(state);
}
void identity_and_gcl(){
 Fixture f;f.C.fill(0);f.a.fill(0);std::fill(f.k.begin(),f.k.end(),0);f.vb=f.va;
 ok(f.run(),"identity");check(f.s==f.ns&&f.t==f.nt,"identity exact populations");
 for(auto& l:f.ledgers)for(int sp=0;sp<6;++sp)
  check(l.spatial_number[sp]==0&&l.spatial_energy_J[sp]==0&&l.work_J[sp]==0&&l.upper_number[sp]==0,"identity zero ledger");
 f.a={0,-.1,-.3,-.2};
 for(int z=0;z<3;++z){f.vb[z]=f.va[z]+f.dt*(f.a[z]-f.a[z+1]);
  for(size_t j=0;j<Fixture::components;++j){f.s[z*Fixture::components+j]=j+1;f.t[z*Fixture::components+j]=.5*(j+1);}}
 for(size_t j=0;j<Fixture::components;++j){f.bs[2*j]=f.bs[2*j+1]=j+1;f.bt[2*j]=f.bt[2*j+1]=.5*(j+1);}
 ok(f.run(),"pure coordinate motion");
 for(size_t i=0;i<Fixture::total;++i){near(f.ns[i],f.s[i],"GCL S uniform density");near(f.nt[i],f.t[i],"GCL T uniform density");}
 for(auto& l:f.ledgers)for(double work:l.work_J)check(work==0,"pure mesh change does not invent work");
 std::cout<<"PASS identity and pure-coordinate GCL with all species and both components\n";
}
void failures(){
 Fixture f;f.C.back()=std::numeric_limits<double>::quiet_NaN();
 std::fill(f.ns.begin(),f.ns.end(),9);std::fill(f.nt.begin(),f.nt.end(),9);for(auto& l:f.ledgers)l.work_J[5]=9;
 check(f.run()!=0,"nonfinite final-zone compression accepted");
 check(std::all_of(f.ns.begin(),f.ns.end(),[](double v){return v==0;})&&std::all_of(f.nt.begin(),f.nt.end(),[](double v){return v==0;}),"all output populations clear");
 for(auto& l:f.ledgers)for(double work:l.work_J)check(work==0,"ledger clearing");
 f.C.back()=.2;f.k.back()=-1;check(f.run()!=0,"invalid last conductance accepted");
 f.k.back()=.01;f.t.back()=-1;check(f.run()!=0,"invalid last T population accepted");
 f.t.back()=1;f.vb.back()=0;check(f.run()!=0,"invalid last volume accepted");
 f.vb.back()=4.3;f.C.back()=.2;
 std::fill(f.s.begin(),f.s.end(),0);std::fill(f.t.begin(),f.t.end(),0);
 std::fill(f.k.begin(),f.k.end(),0);std::fill(f.bs.begin(),f.bs.end(),0);std::fill(f.bt.begin(),f.bt.end(),0);f.a.fill(0);
 for(double& e:f.edges)e*=1.e100;
 f.t.back()=1.e300; // radial N is representable; final-zone energy moment is not.
 check(f.run()==PB11_STATUS_NUMERICAL_FAILURE,"late final-zone work overflow accepted");
 check(std::all_of(f.ns.begin(),f.ns.end(),[](double v){return v==0;}),"partial result published on later failure");
 check(fusion_c_kinetic_geometry_trial(301,1,1,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr)!=0,"dimension limit");
 std::cout<<"PASS invalid late-zone/component inputs and atomic failure clearing\n";
}
}
int main(){try{
 Fixture ordinary;compare_sequential_and_state(ordinary);Fixture si(1.e-16);compare_sequential_and_state(si);
 for(size_t i=0;i<Fixture::total;++i){near(ordinary.ns[i],si.ns[i],"energy-unit covariance S");near(ordinary.nt[i],si.nt[i],"energy-unit covariance T");}
 std::cout<<"PASS actual operator composition, independent SI balances and volume-state batch acceptance\n";
 identity_and_gcl();failures();return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
