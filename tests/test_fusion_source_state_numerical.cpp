#ifdef NDEBUG
#undef NDEBUG
#endif
#include "fusion_source_state.h"
#include "pb11_c.h"
#include <array>
#include <vector>
#include <cassert>
#include <cstdio>
#include <limits>
std::vector<unsigned char> pack(fusion_source_state_v1*p){size_t n=0,k=0;assert(fusion_c_source_state_pack_size(p,&n)==0);std::vector<unsigned char>b(n);assert(fusion_c_source_state_pack(p,b.data(),n,&k)==0&&k==n);return b;}
void check(bool moving){
 double edges[]={0,2,4};std::array<double,12>s{},t{},trial{},outS{},outT{};s[0]=1;trial[0]=moving?.25:.5;trial[1]=moving?.25:.5;
 fusion_source_state_v1*p=nullptr;int st=moving?fusion_c_source_state_create_volume(2,edges,s.data(),t.data(),2,0,106,&p):fusion_c_source_state_create(2,edges,s.data(),t.data(),0,106,&p);assert(st==0);
 auto original=pack(p);assert(original[8]==(moving?3:1));
 fusion_source_ledger_v1 l{},out{};fusion_transport_ledger_v1 tr{},trout{};double inert[6]{},num[6]{},cum[6]{},outI[6]{},time=0,ref=0,vol=0;uint64_t ticket=0,epoch=0;
 auto stage=[&](const double*n){return moving?fusion_c_source_state_stage_volume_numerical(p,ticket,trial.data(),t.data(),4,3,&l,inert,&tr,n):fusion_c_source_state_stage_numerical(p,ticket,trial.data(),t.data(),&l,inert,n);};
 assert(fusion_c_source_state_begin(p,1,&ticket)==0);num[0]=moving?-2./3.:-1;
 double wrong[6]{};assert(stage(wrong)!=0);wrong[0]=1;assert(stage(wrong)!=0);
 assert(stage(num)==0);assert(fusion_c_source_state_discard(p,ticket)==0);auto discarded=pack(p);assert(discarded[8]==(moving?3:1)); // counter advances, format does not
 assert(fusion_c_source_state_begin(p,1,&ticket)==0);assert(stage(num)==0);
 assert(stage(nullptr)!=0);assert(fusion_c_source_state_commit(p,ticket)!=0);assert(stage(num)==0);assert(fusion_c_source_state_commit(p,ticket)==0);assert(fusion_c_source_state_commit(p,ticket)!=0);
 auto accepted=pack(p);assert(accepted[8]==(moving?5:4));
 outS.fill(-7);outT.fill(-7);
 st=moving?fusion_c_source_state_snapshot_volume_numerical(p,outS.data(),outT.data(),&out,outI,&trout,&ref,&vol,&time,&epoch,nullptr):fusion_c_source_state_snapshot_numerical(p,outS.data(),outT.data(),&out,outI,nullptr,&time,&epoch);
 assert(st!=0);for(double x:outS)assert(x==(moving?0.:-7.));assert(pack(p)==accepted);

 auto snapshot=[&](fusion_source_state_v1*q){return moving?fusion_c_source_state_snapshot_volume_numerical(q,outS.data(),outT.data(),&out,outI,&trout,&ref,&vol,&time,&epoch,cum):fusion_c_source_state_snapshot_numerical(q,outS.data(),outT.data(),&out,outI,cum,&time,&epoch);};
 assert(snapshot(p)==0&&cum[0]==-1&&epoch==1&&time==1&&outS==trial);if(moving)assert(ref==2&&vol==4);
 st=moving?fusion_c_source_state_snapshot_volume(p,outS.data(),outT.data(),&out,outI,&trout,&ref,&vol,&time,&epoch):fusion_c_source_state_snapshot_inert(p,outS.data(),outT.data(),&out,outI,&time,&epoch);assert(st!=0);
 fusion_source_state_v1*q=nullptr;assert(fusion_c_source_state_unpack(accepted.data(),accepted.size(),106,&q)==0);assert(pack(q)==accepted);assert(snapshot(q)==0&&cum[0]==-1);
 assert(fusion_c_source_state_begin(p,1,&ticket)==0);st=moving?fusion_c_source_state_stage_volume(p,ticket,trial.data(),t.data(),4,3,&l,inert,&tr):fusion_c_source_state_stage_inert(p,ticket,trial.data(),t.data(),&l,inert);assert(st!=0);assert(fusion_c_source_state_commit(p,ticket)!=0);
 num[0]=std::numeric_limits<double>::quiet_NaN();assert(stage(num)!=0);assert(fusion_c_source_state_discard(p,ticket)==0);assert(snapshot(p)==0&&cum[0]==-1&&epoch==1);
 auto corrupt=accepted;corrupt[corrupt.size()/2]^=1;fusion_source_state_v1*bad=nullptr;assert(fusion_c_source_state_unpack(corrupt.data(),corrupt.size(),106,&bad)!=0&&bad==nullptr);
 fusion_c_source_state_destroy(q);fusion_c_source_state_destroy(p);
}
int main(){check(false);check(true);puts("PASS numerical state: stationary/moving accounting, reject/discard/commit, legacy API exclusion, restart v4/v5");}
