#include "../src/fusion_boundary_energy.hpp"
#include "fusion_kinetic_geometry.h"
#include "fusion_source_state.h"
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>
void check(bool v,const char*m){if(!v)throw std::runtime_error(m);}
void rounding(){
 using fusion_boundary_energy::subnormal_product;
 const double d=std::numeric_limits<double>::denorm_min(),m=std::numeric_limits<double>::min();
 struct C{double n,e,w;};
 const C cases[]={{d,.25,0},{d,.5,0},{d,std::nextafter(.5,1.),d},{3*d,.5,2*d},
 {5*d,.5,2*d},{m,.5,m/2},{std::nextafter(m,0.),1.,std::nextafter(m,0.)},
 {m,std::nextafter(1.,0.),m},{3.2114266979681025e-322,4e-12,0}};
 for(int mode:{FE_TONEAREST,FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO}){
  std::fesetround(mode);
  for(auto c:cases){double got=-1;check(subnormal_product(c.n,c.e,got),"classification");check(got==c.w&&!std::signbit(got),"unique nearest-even bits");}
 }
 std::fesetround(FE_TONEAREST);double out=42;
 check(!subnormal_product(m,1.,out)&&out==42,"normal unchanged");
 check(!subnormal_product(-1,1,out)&&!subnormal_product(INFINITY,1,out),"invalid classification");
}
void scenario(double tail,double volume){
 double edge[5]={0,1e-12,2e-12,3e-12,4e-12},v[1]={volume},a[2]={},c[1]={-1e-100};
 std::array<double,48> k{},bs{},bt{};std::array<double,24>s{},t{},ns{},nt{};
 s[16]=1;s[19]=tail;fusion_transport_ledger_v1 movement{};
 check(!fusion_c_kinetic_geometry_trial(1,4,.1,edge,v,v,a,c,k.data(),bs.data(),bt.data(),s.data(),t.data(),ns.data(),nt.data(),&movement),"geometry");
 fusion_source_state_v1*p=nullptr;check(!fusion_c_source_state_create_volume(4,edge,s.data(),t.data(),volume,0,91,&p),"create");
 uint64_t ticket=0;check(!fusion_c_source_state_begin(p,.1,&ticket),"begin");
 fusion_source_ledger_v1 source{};double inert[6]={};
 auto stage=[&](const fusion_transport_ledger_v1& x){return fusion_c_source_state_stage_volume(p,ticket,ns.data(),nt.data(),volume,volume,&source,inert,&x);};
 auto bad=movement;double rounded;
 if(fusion_boundary_energy::subnormal_product(movement.upper_number[4],edge[4],rounded)){
  check(movement.upper_energy_J[4]==rounded,"producer canonical");
  bad.upper_energy_J[4]=std::nextafter(rounded,INFINITY);check(stage(bad)!=0,"adjacent wrong subnormal rejected");
  if(rounded>0){bad.upper_energy_J[4]=0;check(stage(bad)!=0,"false zero rejected");}
 }else {bad.upper_energy_J[4]*=1.001;check(stage(bad)!=0,"wrong normal energy rejected");}
 bad=movement;bad.upper_energy_J[4]=NAN;check(stage(bad)!=0,"nan rejected");
 bad=movement;bad.upper_number[4]=-1;check(stage(bad)!=0,"negative rejected");
 check(stage(movement)==0,"valid stage");check(!fusion_c_source_state_commit(p,ticket),"commit");
 size_t size=0,written=0;check(!fusion_c_source_state_pack_size(p,&size),"size");std::vector<unsigned char> bytes(size);
 check(!fusion_c_source_state_pack(p,bytes.data(),size,&written),"pack");fusion_source_state_v1*q=nullptr;
 check(!fusion_c_source_state_unpack(bytes.data(),size,91,&q),"unpack");
 std::vector<unsigned char> copy(size);check(!fusion_c_source_state_pack(q,copy.data(),size,&written)&&copy==bytes,"restart exact");
 fusion_c_source_state_destroy(q);fusion_c_source_state_destroy(p);
}
void cumulative(double loss,int steps){
 double edges[]={0,1e-12,2e-12,3e-12,4e-12};std::array<double,24>s{},t{};s[16]=1;
 fusion_source_state_v1*p=nullptr;check(!fusion_c_source_state_create_volume(4,edges,s.data(),t.data(),1,0,91,&p),"cumulative create");
 fusion_source_ledger_v1 source{};double inert[6]={};fusion_transport_ledger_v1 move{};
 move.upper_number[4]=loss;move.upper_energy_J[4]=loss*edges[4];
 fusion_boundary_energy::canonicalize(loss,edges[4],move.upper_energy_J[4]);
 double total=0,expected_u=0;
 for(int i=0;i<steps;++i){
  uint64_t ticket=0;check(!fusion_c_source_state_begin(p,.1,&ticket),"cumulative begin");
  check(!fusion_c_source_state_stage_volume(p,ticket,s.data(),t.data(),1,1,&source,inert,&move),"cumulative stage");
  check(!fusion_c_source_state_commit(p,ticket),"cumulative commit");
  total=double((long double)total+loss);expected_u=double((long double)expected_u+move.upper_energy_J[4]);
  fusion_boundary_energy::canonicalize(total,edges[4],expected_u);
  fusion_transport_ledger_v1 actual{};double reference,volume,time;uint64_t epoch;
  check(!fusion_c_source_state_snapshot_volume(p,s.data(),t.data(),&source,inert,&actual,&reference,&volume,&time,&epoch),"cumulative snapshot");
  check(actual.upper_number[4]==total&&actual.upper_energy_J[4]==expected_u,"canonical cumulative value");
 }
 double final_exact=double((long double)total*edges[4]);check(final_exact>0,"nonzero accumulated energy");
 size_t n,w;check(!fusion_c_source_state_pack_size(p,&n),"cumulative pack size");std::vector<unsigned char>b(n);
 check(!fusion_c_source_state_pack(p,b.data(),n,&w),"cumulative pack");fusion_source_state_v1*q=nullptr;
 check(!fusion_c_source_state_unpack(b.data(),n,91,&q),"cumulative restart");fusion_c_source_state_destroy(q);fusion_c_source_state_destroy(p);
}
int main(){try{rounding();cumulative(std::numeric_limits<double>::denorm_min()/16e-12,5);cumulative(std::numeric_limits<double>::min()*.4/4e-12,3);for(double n:{1e-180,1e-190,1e-200,1e-210,1e-215,1e-220})for(double v:{.08,1.,1e4})scenario(n,v);std::puts("PASS boundary exact rounding and 18 public geometry/stage/restart scenarios");}catch(const std::exception&e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
