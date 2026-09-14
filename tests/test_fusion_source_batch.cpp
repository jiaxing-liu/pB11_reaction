#include "fusion_source_state.h"
#include <array>
#include <iostream>
#include <stdexcept>
static void require(bool b){if(!b)throw std::runtime_error("batch transaction check failed");}
struct Context {
 fusion_source_state_v1* p=nullptr;std::array<double,12> s{},t{};int cells;double total=0;uint64_t ticket=0;
 Context(int n,double time=0):cells(n){double edges[]={0,2e-15,4e-15};require(fusion_c_source_state_create(n,edges,s.data(),t.data(),time,n,&p)==0);}
 ~Context(){fusion_c_source_state_destroy(p);}
 void begin(double dt){require(fusion_c_source_state_begin(p,dt,&ticket)==0);}
 void stage(double addition,bool extended=false){s[0]=total+addition;fusion_source_ledger_v1 l{};l.external_born_number_m3[0]=addition;l.external_born_energy_J_m3[0]=addition*1e-15;double inert[6]{};
 require((extended?fusion_c_source_state_stage_inert(p,ticket,s.data(),t.data(),&l,inert):fusion_c_source_state_stage(p,ticket,s.data(),t.data(),&l))==0);}
 void check(double count,double time,uint64_t ep){std::array<double,12>a{},b{};double heat[6],tm;uint64_t e;fusion_source_ledger_v1 l{};require(fusion_c_source_state_snapshot_inert(p,a.data(),b.data(),&l,heat,&tm,&e)==0);require(a[0]==count&&tm==time&&e==ep);require(l.external_born_number_m3[0]==count);}
};
int main(){try{
 Context a(1),b(2);a.begin(.1);b.begin(.1);a.stage(4);b.stage(3,true);
 fusion_source_state_v1* states[]={a.p,b.p};uint64_t tickets[]={a.ticket,b.ticket};
 auto unchanged=[&](){a.check(0,0,0);b.check(0,0,0);};
 tickets[1]++;require(fusion_c_source_state_commit_many(2,states,tickets)!=0);unchanged();tickets[1]=b.ticket;
 fusion_source_state_v1* duplicates[]={a.p,a.p};require(fusion_c_source_state_commit_many(2,duplicates,tickets)!=0);unchanged();
 fusion_source_state_v1* null_last[]={a.p,nullptr};require(fusion_c_source_state_commit_many(2,null_last,tickets)!=0);unchanged();
 require(fusion_c_source_state_commit_many(0,states,tickets)!=0);require(fusion_c_source_state_commit_many(2,nullptr,tickets)!=0);require(fusion_c_source_state_commit_many(2,states,nullptr)!=0);unchanged();
 b.begin(.2);b.stage(3,true);tickets[1]=b.ticket;require(fusion_c_source_state_commit_many(2,states,tickets)!=0);unchanged();
 b.begin(.1);tickets[1]=b.ticket;require(fusion_c_source_state_commit_many(2,states,tickets)!=0);unchanged();
 b.stage(3,true);require(fusion_c_source_state_commit_many(2,states,tickets)==0);a.total=4;b.total=3;a.check(4,.1,1);b.check(3,.1,1);
 require(fusion_c_source_state_commit_many(2,states,tickets)!=0);a.check(4,.1,1);b.check(3,.1,1);
 // Equal accepted/trial times with a different accepted epoch must reject.
 Context c(1,.1);a.begin(.1);a.stage(2);c.begin(.1);c.stage(5);states[1]=c.p;tickets[0]=a.ticket;tickets[1]=c.ticket;
 require(fusion_c_source_state_commit_many(2,states,tickets)!=0);a.check(4,.1,1);c.check(0,.1,0);
 // The rejected group did not invalidate the good candidate.
 require(fusion_c_source_state_commit(a.p,a.ticket)==0);a.check(6,.2,2);
 std::cout<<"PASS: atomic multi-context acceptance and rejection preservation\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
