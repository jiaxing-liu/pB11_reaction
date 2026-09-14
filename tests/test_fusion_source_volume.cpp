#include "fusion_source_state.h"
#include "fusion_nuclear_data.h"
#include "fusion_radial_transport.h"
#include "fusion_energy_work.h"
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Pop=std::array<double,12>;
void check(bool v,const char* msg){if(!v)throw std::runtime_error(msg);}
void ok(int s,const char* msg){if(s!=0)throw std::runtime_error(std::string(msg)+": status="+std::to_string(s));}
void bad(int s,const char* msg){check(s!=0,msg);}
void near(double a,double b,const char* msg){check(std::abs(a-b)<=1e-11*std::max(std::abs(a),std::abs(b)),msg);}
struct State {
 fusion_source_state_v1* p=nullptr;
 State()=default;State(const State&)=delete;State& operator=(const State&)=delete;
 ~State(){fusion_c_source_state_destroy(p);}
 void create(double v=2,double density=4){double e[]={0,2,4};Pop s{},t{};s[0]=density;
  ok(fusion_c_source_state_create_volume(2,e,s.data(),t.data(),v,0,91,&p),"create volume");}
};
struct Snap {
 Pop s{},t{};fusion_source_ledger_v1 source{};fusion_transport_ledger_v1 transport{};
 std::array<double,6> inert{};double reference=0,volume=0,time=0;uint64_t epoch=0;
 explicit Snap(fusion_source_state_v1* p){ok(fusion_c_source_state_snapshot_volume(p,s.data(),t.data(),&source,
  inert.data(),&transport,&reference,&volume,&time,&epoch),"snapshot volume");}
};
std::vector<unsigned char> pack(fusion_source_state_v1* p){size_t n=0,w=0;ok(fusion_c_source_state_pack_size(p,&n),"pack size");
 std::vector<unsigned char> b(n);ok(fusion_c_source_state_pack(p,b.data(),n,&w),"pack");check(n==w,"packed length");return b;}
void step(State& s,double lo,double hi,double volume=1,double source_volume=1,
 const fusion_source_ledger_v1& source={},const fusion_transport_ledger_v1& movement={},double inert_heat=0){
 uint64_t ticket=0;Pop candidate{},t{};candidate[0]=lo;candidate[1]=hi;
 std::array<double,6> inert{};inert[0]=inert_heat;
 ok(fusion_c_source_state_begin(s.p,.125,&ticket),"begin");
 ok(fusion_c_source_state_stage_volume(s.p,ticket,candidate.data(),t.data(),volume,source_volume,
  &source,inert.data(),&movement),"stage volume");
 ok(fusion_c_source_state_commit(s.p,ticket),"commit");
}
void checksum(std::vector<unsigned char>& b){uint64_t h=UINT64_C(14695981039346656037);
 for(size_t j=0;j<b.size()-8;++j){h^=b[j];h*=UINT64_C(1099511628211);}
 for(size_t j=0;j<8;++j)b[b.size()-8+j]=static_cast<unsigned char>(h>>(8*j));}
void word(std::vector<unsigned char>& b,size_t index,double value){uint64_t bits;std::memcpy(&bits,&value,8);
 for(size_t j=0;j<8;++j)b[8*index+j]=static_cast<unsigned char>(bits>>(8*j));checksum(b);}
void restart_rejected(std::vector<unsigned char> b,size_t offset,double value){word(b,offset,value);State restored;
 bad(fusion_c_source_state_unpack(b.data(),b.size(),91,&restored.p),"semantic restart corruption accepted");
 check(restored.p==nullptr,"failed unpack output");}
void analytic_sequence(){
 State state;state.create();auto fresh=pack(state.p);check(fresh[8]==3,"restart version3");
 State fresh_copy;ok(fusion_c_source_state_unpack(fresh.data(),fresh.size(),91,&fresh_copy.p),"fresh restart");
 check(pack(fresh_copy.p)==fresh,"fresh restart bytes");
 constexpr size_t volume_offset=8+12+121+6;
 restart_rejected(fresh,volume_offset,0);restart_rejected(fresh,volume_offset+1,1);
 restart_rejected(fresh,volume_offset+2,1); // no transport before any accepted epoch
 step(state,8,0); // actual V: 2 -> 1, N=8 and U=8 unchanged
 Snap diluted(state.p);check(diluted.reference==2&&diluted.volume==1&&diluted.s[0]==8,"volume dilution");
 fusion_transport_ledger_v1 movement{};movement.work_J[0]=8;
 step(state,4,4,1,1,{},movement); // U:8 ->16, N unchanged
 movement={};movement.spatial_number[0]=-2;movement.spatial_energy_J[0]=-6;
 step(state,4,2,1,1,{},movement); // outward particles of mean energy3
 movement={};movement.upper_number[0]=1;movement.upper_energy_J[0]=4;movement.work_J[0]=1;
 step(state,4,1,1,1,{},movement); // numerical boundary outflow, not ash or physical escape
 fusion_source_ledger_v1 source{};source.external_born_number_m3[0]=1;source.external_born_energy_J_m3[0]=1;
 step(state,6,1,1,2,source); // density source1 * source-volume2 => extensive2
 step(state,6.5,.5,1,2,{}, {},.5); // inert heat .5 * source-volume2 =1J
 Snap final(state.p);check(final.epoch==6&&final.time==.75,"accepted epoch/time");
 check(final.source.external_born_number_m3[0]==1&&final.source.external_born_energy_J_m3[0]==1,"reference source density");
 check(final.inert[0]==.5&&final.transport.work_J[0]==9&&final.transport.spatial_energy_J[0]==-6,"separate heat/work ledgers");
 check(final.source.escaped_number_m3[0]==0&&final.source.handed_off_number_m3[0]==0&&
  final.source.nuclear_born_number_m3[0]==0,"no implicit loss/ash/birth");
 const auto bytes=pack(state.p);State restored;ok(fusion_c_source_state_unpack(bytes.data(),bytes.size(),91,&restored.p),"evolved restart");
 check(pack(restored.p)==bytes,"evolved restart bytes");
 restart_rejected(bytes,volume_offset+1,2); // volume rescale breaks total inventory
 restart_rejected(bytes,volume_offset+2+12,10); // unmatched extra work
 restart_rejected(bytes,volume_offset+2+36,3); // upper energy no longer equals edge*N
 step(state,6.5,.5);step(restored,6.5,.5);check(pack(state.p)==pack(restored.p),"restart continuation identity");
 // Opposite energy-bin face flows can have exactly zero net number exchange.
 movement={};movement.spatial_energy_J[0]=1;
 step(state,6,1,1,1,{},movement);
 check(Snap(state.p).transport.spatial_number[0]==-2,"zero net number with nonzero spatial energy");
 std::cout<<"PASS analytic volume/work/spatial/domain/source/inert sequence and v3 restart\n";
}
void failure_and_atomicity(){
 State a,b;a.create();b.create();auto initial=pack(a.p);
 Pop candidate{},t{};candidate[0]=8;std::array<double,6> inert{};
 fusion_source_ledger_v1 source{};fusion_transport_ledger_v1 movement{};
 uint64_t tickets[2]{};fusion_source_state_v1* handles[]={a.p,b.p};
 ok(fusion_c_source_state_begin(a.p,.125,&tickets[0]),"begin a");
 ok(fusion_c_source_state_begin(b.p,.125,&tickets[1]),"begin b");
 auto stage=[&](State& p,uint64_t ticket,double volume,double source_volume){return fusion_c_source_state_stage_volume(
  p.p,ticket,candidate.data(),t.data(),volume,source_volume,&source,inert.data(),&movement);};
 ok(stage(a,tickets[0],1,1),"stage a");
 candidate[0]=7;bad(stage(b,tickets[1],1,1),"reject number mismatch");
 bad(fusion_c_source_state_commit_many(2,handles,tickets),"batch accepts failed zone");
 check(Snap(a.p).volume==2&&Snap(a.p).epoch==0&&Snap(a.p).s[0]==4,"no partial batch volume publication");
 candidate[0]=8;ok(stage(b,tickets[1],1,1),"restage b");
 ok(fusion_c_source_state_commit_many(2,handles,tickets),"batch publish");
 check(Snap(a.p).volume==1&&Snap(b.p).volume==1,"all volumes published");
 ok(fusion_c_source_state_begin(a.p,.125,&tickets[0]),"second begin");
 for(double invalid:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){
  ok(stage(a,tickets[0],1,1),"valid candidate before invalid replacement");
  bad(stage(a,tickets[0],invalid,1),"invalid volume");
  bad(fusion_c_source_state_commit(a.p,tickets[0]),"stale stage commits after volume failure");
  ok(stage(a,tickets[0],1,1),"restage");bad(stage(a,tickets[0],1,invalid),"invalid source volume");
  bad(fusion_c_source_state_commit(a.p,tickets[0]),"stale stage commits after source volume failure");
 }
 // An invalid tiny negative source must not disappear during normalization.
 source.external_born_number_m3[0]=-std::numeric_limits<double>::denorm_min();
 bad(stage(a,tickets[0],1,.25),"invalid pre-scaled source accepted");source={};
 source.external_born_number_m3[0]=std::numeric_limits<double>::denorm_min();
 bad(stage(a,tickets[0],1,.25),"nonzero source underflow silently erased");source={};
 movement.upper_number[0]=1;movement.upper_energy_J[0]=3;
 bad(stage(a,tickets[0],1,1),"wrong upper edge energy");movement={};
 movement.lower_number[0]=-1;bad(stage(a,tickets[0],1,1),"negative domain flow");movement={};
 ok(stage(a,tickets[0],1,1),"valid before legacy stage");
 bad(fusion_c_source_state_stage_inert(a.p,tickets[0],candidate.data(),t.data(),&source,inert.data()),"legacy inert stage on volume context");
 bad(fusion_c_source_state_commit(a.p,tickets[0]),"legacy call retains stale stage");
 double time=9;uint64_t epoch=9;
 bad(fusion_c_source_state_snapshot_inert(a.p,candidate.data(),t.data(),&source,inert.data(),&time,&epoch),"ambiguous legacy snapshot");
 check(time==0&&epoch==0,"legacy failure outputs");
 ok(fusion_c_source_state_discard(a.p,tickets[0]),"discard");
 check(Snap(a.p).volume==1&&Snap(a.p).s[0]==8,"failed stages preserve accepted state");
 // A rejected/discarded trial cannot enter restart, and fresh legacy v1 is unchanged.
 State legacy;double edges[]={0,2,4};candidate={};candidate[0]=4;
 ok(fusion_c_source_state_create(2,edges,candidate.data(),t.data(),0,91,&legacy.p),"legacy create");
 check(pack(legacy.p)[8]==1,"legacy format unchanged");
 double ref=9,vol=9;bad(fusion_c_source_state_snapshot_volume(legacy.p,candidate.data(),t.data(),&source,inert.data(),&movement,&ref,&vol,&time,&epoch),"legacy volume snapshot");
 check(ref==0&&vol==0&&time==0&&epoch==0,"volume failure clears scalars");
 for(double x:candidate)check(x==0,"volume failure clears population");
 State invalid;bad(fusion_c_source_state_create_volume(2,edges,candidate.data(),t.data(),0,0,91,&invalid.p),"zero initial volume");
 check(invalid.p==nullptr,"invalid create handle");
 std::cout<<"PASS normalization rejection, invalid-stage invalidation, atomic volume batch and legacy guards\n";
}
void nuclear_source_normalization(){
 fusion_nuclear_channel_v1 dt{};ok(fusion_c_nuclear_channel(FUSION_DT_ALPHAN,&dt),"DT nuclear data");
 const double e=dt.q_J/5,edges[]={0,2*e,4*e};Pop zero{},candidate{};
 State state;ok(fusion_c_source_state_create_volume(2,edges,zero.data(),zero.data(),2,0,91,&state.p),"nuclear volume context");
 fusion_source_ledger_v1 source{};fusion_transport_ledger_v1 movement{};double inert[6]{};
 source.events_m3[FUSION_DT_ALPHAN]=1;source.thermal_consumed_number_m3[FUSION_DEUTERON]=1;
 source.thermal_consumed_number_m3[FUSION_TRITON]=1;
 source.nuclear_born_number_m3[FUSION_HELIUM4]=1;source.nuclear_born_energy_J_m3[FUSION_HELIUM4]=e;
 source.neutron_number_m3=1;source.neutron_energy_J_m3=dt.q_J-e;
 // Synthetic budget-only energies: tests normalization, not DT kinematic partition.
 candidate[2*FUSION_HELIUM4]=3;uint64_t ticket=0;
 ok(fusion_c_source_state_begin(state.p,.125,&ticket),"nuclear begin");
 ok(fusion_c_source_state_stage_volume(state.p,ticket,candidate.data(),zero.data(),1,3,&source,inert,&movement),"nuclear SI scaled stage");
 ok(fusion_c_source_state_commit(state.p,ticket),"nuclear commit");
 Snap result(state.p);near(result.source.events_m3[FUSION_DT_ALPHAN]*result.reference,3,"extensive event count");
 near((result.source.nuclear_born_energy_J_m3[FUSION_HELIUM4]+result.source.neutron_energy_J_m3)*result.reference,
  3*dt.q_J,"extensive nuclear energy");
 check(result.s[2*FUSION_HELIUM4]==3&&result.source.handed_off_number_m3[FUSION_HELIUM4]==0,"fast He4 inventory remains kinetic");
 std::cout<<"PASS SI DT ledger source-volume normalization (synthetic budget-only energy split)\n";
}

void compose_actual_operators(){
 // Two closed radial zones, two energy bins. Actual mesh exchange is shared
 // once at the internal face. Energy drift is a separately prescribed operator.
 const double edges[]={0,2,4},vold[]={2,3},vnew[]={1.9,3.1};
 const double advection[]={0,.2,0},conductance[6]{},boundary[4]{};
 const double old[]={4,1,1,3}; // energy-major, zones contiguous
 double radial[4]{},radial_face[6]{};fusion_radial_ledger_v1 radial_ledger[2]{};
 ok(fusion_c_radial_transport_trial(2,2,.5,vold,vnew,advection,conductance,
  old,boundary,radial,radial_face,radial_ledger),"real radial operator");
 State a,b;State* states[]={&a,&b};uint64_t tickets[2]{};
 double total_initial_n=0,total_initial_u=0,total_final_n=0,total_final_u=0;
 double total_work=0,total_out_n=0,total_out_u=0,spatial_n=0,spatial_u=0;
 for(int z=0;z<2;++z){
  Pop initial{},t{},candidate{};initial[0]=old[z];initial[1]=old[2+z];
  ok(fusion_c_source_state_create_volume(2,edges,initial.data(),t.data(),vold[z],0,91,&states[z]->p),"operator context");
  const double radial_zone[]={radial[z],radial[2+z]};double work_face[3]{};
  fusion_energy_work_ledger_v1 work{};
  ok(fusion_c_energy_work_trial(2,.5,z==0?.1:-.2,edges,radial_zone,candidate.data(),work_face,&work),"real work operator");
  fusion_transport_ledger_v1 movement{};fusion_source_ledger_v1 source{};std::array<double,6> inert{};
  for(int e=0;e<2;++e){const double inward=radial_face[e*3+z]-radial_face[e*3+z+1];
   movement.spatial_number[0]+=inward;movement.spatial_energy_J[0]+=inward*(e==0?1:3);}
  movement.work_J[0]=vnew[z]*work.work_on_particles_J_m3;
  movement.lower_number[0]=vnew[z]*work.lower_number_m3;
  movement.lower_energy_J[0]=vnew[z]*work.lower_energy_J_m3;
  movement.upper_number[0]=vnew[z]*work.upper_number_m3;
  movement.upper_energy_J[0]=vnew[z]*work.upper_energy_J_m3;
  ok(fusion_c_source_state_begin(states[z]->p,.5,&tickets[z]),"operator begin");
  ok(fusion_c_source_state_stage_volume(states[z]->p,tickets[z],candidate.data(),t.data(),vnew[z],vnew[z],
   &source,inert.data(),&movement),"operator ledger stage");
  total_initial_n+=vold[z]*(initial[0]+initial[1]);total_initial_u+=vold[z]*(initial[0]+3*initial[1]);
  total_final_n+=vnew[z]*(candidate[0]+candidate[1]);total_final_u+=vnew[z]*(candidate[0]+3*candidate[1]);
  total_work+=movement.work_J[0];total_out_n+=movement.lower_number[0]+movement.upper_number[0];
  total_out_u+=movement.lower_energy_J[0]+movement.upper_energy_J[0];
  spatial_n+=movement.spatial_number[0];spatial_u+=movement.spatial_energy_J[0];
 }
 fusion_source_state_v1* handles[]={a.p,b.p};ok(fusion_c_source_state_commit_many(2,handles,tickets),"operator batch");
 check(spatial_n==0&&spatial_u==0,"shared internal face cancellation");
 near(total_final_n+total_out_n,total_initial_n,"composed global number");
 near(total_final_u+total_out_u,total_initial_u+total_work,"composed global work budget");
 check(total_out_n>0,"heating produces separately tracked upper outflow");
 std::cout<<"PASS actual radial + energy-work composition: N0="<<total_initial_n<<" N1="<<total_final_n
  <<" grid_out_N="<<total_out_n<<" U0="<<total_initial_u<<" U1="<<total_final_u
  <<" work="<<total_work<<" grid_out_U="<<total_out_u<<'\n';
}

}
int main(){try{static_assert(sizeof(fusion_transport_ledger_v1)==42*sizeof(double),"C layout");
 analytic_sequence();failure_and_atomicity();nuclear_source_normalization();compose_actual_operators();return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
