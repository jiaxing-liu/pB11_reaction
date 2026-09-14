#include "fusion_source_state.h"
#include "fusion_nuclear_data.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
constexpr int OK=PB11_STATUS_OK, BAD=PB11_STATUS_INVALID_ARGUMENT;
constexpr int NUM=PB11_STATUS_NUMERICAL_FAILURE, EXC=PB11_STATUS_EXCEPTION;
constexpr size_t ledger_count=121;
using Moments=std::array<long double,6>;
// Field traversal is explicit: no assumptions about C struct padding.
template<class L,class F> void fields(L& x,F f) {
 for(auto& v:x.events_m3)f(v);
 for(auto& v:x.nuclear_born_number_m3)f(v);
 for(auto& v:x.nuclear_born_energy_J_m3)f(v);
 for(auto& v:x.external_born_number_m3)f(v);
 for(auto& v:x.external_born_energy_J_m3)f(v);
 for(auto& v:x.thermal_consumed_number_m3)f(v);
 for(auto& v:x.thermal_consumed_energy_J_m3)f(v);
 for(auto& v:x.fast_consumed_number_m3)f(v);
 for(auto& v:x.fast_consumed_energy_J_m3)f(v);
 for(auto& v:x.escaped_number_m3)f(v);
 for(auto& v:x.escaped_energy_J_m3)f(v);
 for(auto& v:x.handed_off_number_m3)f(v);
 for(auto& v:x.handed_off_energy_J_m3)f(v);
 for(auto& v:x.heat_to_bath_J_m3)f(v);
 f(x.neutron_number_m3);f(x.neutron_energy_J_m3);
}
bool representable(long double x) {
 return std::isfinite(x)&&std::abs(x)<=std::numeric_limits<double>::max();
}
bool balance(std::initializer_list<long double> terms) {
 long double sum=0,scale=0;
 for(auto x:terms){if(!std::isfinite(x))return false;sum+=x;scale+=std::abs(x);}
 return std::isfinite(sum)&&std::isfinite(scale)&&std::abs(sum)<=1e-10L*scale;
}
bool ledger_valid(const fusion_source_ledger_v1& x) {
 size_t k=0; bool valid=true;
 fields(x,[&](double v){const bool heat=k>=77&&k<119; ++k;
  if(!std::isfinite(v)||(!heat&&v<0))valid=false;});
 if(!valid)return false;
 auto pair=[](double n,double u){return n!=0||u==0;};
 for(int s=0;s<6;++s)
  if(!pair(x.nuclear_born_number_m3[s],x.nuclear_born_energy_J_m3[s])||
     !pair(x.external_born_number_m3[s],x.external_born_energy_J_m3[s])||
     !pair(x.thermal_consumed_number_m3[s],x.thermal_consumed_energy_J_m3[s])||
     !pair(x.fast_consumed_number_m3[s],x.fast_consumed_energy_J_m3[s])||
     !pair(x.escaped_number_m3[s],x.escaped_energy_J_m3[s])||
     !pair(x.handed_off_number_m3[s],x.handed_off_energy_J_m3[s]))return false;
 if(!pair(x.neutron_number_m3,x.neutron_energy_J_m3))return false;
 Moments born{},consumed{}; long double neutrons=0,q=0;
 for(int c=0;c<5;++c){fusion_nuclear_channel_v1 channel{};
  if(fusion_c_nuclear_channel(c,&channel)!=OK)return false;
  long double n=x.events_m3[c];q+=n*channel.q_J;
  for(int id:channel.reactant_ids)consumed[id]+=n;
  for(int j=0;j<channel.product_count;++j){int id=channel.product_ids[j];
   if(id==FUSION_MASS_NEUTRON)neutrons+=n;else born[id]+=n;}
 }
 long double born_u=x.neutron_energy_J_m3,removed_u=0;
 for(int s=0;s<6;++s){
  if(!balance({born[s],-static_cast<long double>(x.nuclear_born_number_m3[s])})||
     !balance({consumed[s],-static_cast<long double>(x.thermal_consumed_number_m3[s]),
               -static_cast<long double>(x.fast_consumed_number_m3[s])}))return false;
  born_u+=x.nuclear_born_energy_J_m3[s];
  removed_u+=static_cast<long double>(x.thermal_consumed_energy_J_m3[s])+x.fast_consumed_energy_J_m3[s];
 }
 return balance({neutrons,-static_cast<long double>(x.neutron_number_m3)})&&
        balance({born_u,-removed_u,-q});
}
bool add_ledger(const fusion_source_ledger_v1& a,const fusion_source_ledger_v1& b,
 fusion_source_ledger_v1& out){
 std::array<double,ledger_count> av{},bv{}; size_t k=0;
 fields(a,[&](double x){av[k++]=x;});k=0;fields(b,[&](double x){bv[k++]=x;});
 k=0;bool valid=true;fields(out,[&](double& x){long double v=static_cast<long double>(av[k])+bv[k];++k;
  if(!representable(v)){valid=false;x=0;}else x=static_cast<double>(v);});
 return valid&&ledger_valid(out);
}
bool inventory(const std::vector<double>& edges,const double* s,const double* t,
 Moments& n,Moments& u) {
 if(!s||!t)return false;
 size_t cells=edges.size()-1;n={};u={};
 for(size_t species=0;species<6;++species)for(size_t i=0;i<cells;++i){
  size_t j=species*cells+i;
  if(!std::isfinite(s[j])||s[j]<0||!std::isfinite(t[j])||t[j]<0)return false;
  long double pop=static_cast<long double>(s[j])+t[j];
  long double center=(static_cast<long double>(edges[i])+edges[i+1])/2;
  n[species]+=pop;u[species]+=pop*center;
 }
 for(int j=0;j<6;++j)if(!representable(n[j])||!representable(u[j]))return false;
 return true;
}
bool inventory_balance(const Moments& old_n,const Moments& old_u,
 const Moments& n,const Moments& u,const fusion_source_ledger_v1& x,
 const std::array<double,6>& inert){
 for(int s=0;s<6;++s){long double heat=inert[s],heat_scale=std::abs(inert[s]);
  for(int b=0;b<7;++b){heat+=x.heat_to_bath_J_m3[7*s+b];heat_scale+=std::abs(x.heat_to_bath_J_m3[7*s+b]);}
  if(!balance({n[s],-old_n[s],-static_cast<long double>(x.nuclear_born_number_m3[s]),
    -static_cast<long double>(x.external_born_number_m3[s]),x.fast_consumed_number_m3[s],
    x.escaped_number_m3[s],x.handed_off_number_m3[s]}))return false;
  // Keep the absolute bath terms in the scale even if different baths cancel.
  if(!balance({u[s],-old_u[s],-static_cast<long double>(x.nuclear_born_energy_J_m3[s]),
    -static_cast<long double>(x.external_born_energy_J_m3[s]),x.fast_consumed_energy_J_m3[s],
    x.escaped_energy_J_m3[s],x.handed_off_energy_J_m3[s],(heat_scale+heat)/2,-(heat_scale-heat)/2}))return false;
 }
 return true;
}
}
struct fusion_source_state_v1 {
 std::vector<double> edges,s,t,trial_s,trial_t;
 Moments initial_n{},initial_u{};
 fusion_source_ledger_v1 cumulative{},staged_cumulative{};
 std::array<double,6> inert{},staged_inert{};
 bool extended=false,staged_extended=false;
 double time=0,initial_time=0,trial_time=0;
 uint64_t tag=0,epoch=0,counter=0;
 bool pending=false,staged=false;
};

extern "C" int fusion_c_source_state_create(int cells,const double* edges,
 const double* s,const double* t,double time,uint64_t tag,fusion_source_state_v1** out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out=nullptr;
 try{
  if(cells<1||cells>1000000||!edges||!s||!t||!std::isfinite(time)||time<0)return BAD;
  for(int i=0;i<=cells;++i)if(!std::isfinite(edges[i])||edges[i]<0||(i&&edges[i]<=edges[i-1]))return BAD;
  auto p=std::make_unique<fusion_source_state_v1>();p->edges.assign(edges,edges+cells+1);
  if(!inventory(p->edges,s,t,p->initial_n,p->initial_u))return BAD;
  // Restart anchors use the same portable double representation from creation.
  for(int j=0;j<6;++j){
   const double anchor_n=static_cast<double>(p->initial_n[j]);
   const double anchor_u=static_cast<double>(p->initial_u[j]);
   if(!balance({p->initial_n[j],-static_cast<long double>(anchor_n)})||
      !balance({p->initial_u[j],-static_cast<long double>(anchor_u)}))return NUM;
   p->initial_n[j]=anchor_n;p->initial_u[j]=anchor_u;
  }
  p->s.assign(s,s+6*size_t(cells));p->t.assign(t,t+6*size_t(cells));
  p->time=p->initial_time=time;p->tag=tag;*out=p.release();return OK;
 }catch(...){return EXC;}
}
extern "C" void fusion_c_source_state_destroy(fusion_source_state_v1* p){delete p;}
extern "C" int fusion_c_source_state_cells(const fusion_source_state_v1* p,int* cells){
 if(!cells)return PB11_STATUS_NULL_OUTPUT;
 *cells=0;if(!p)return BAD;
 *cells=static_cast<int>(p->edges.size()-1);return OK;
}
namespace {
int snapshot(const fusion_source_state_v1* p,double* s,double* t,
 fusion_source_ledger_v1* ledger,double* inert,double* time,uint64_t* epoch,
 bool extended_output){
 if(ledger)*ledger={};
 if(inert)std::fill(inert,inert+6,0.);
 if(time)*time=0;
 if(epoch)*epoch=0;
 if(!p||!s||!t||!ledger||!time||!epoch||
    (extended_output&&!inert)||(!extended_output&&p->extended))return BAD;
 std::copy(p->s.begin(),p->s.end(),s);std::copy(p->t.begin(),p->t.end(),t);
 if(inert)std::copy(p->inert.begin(),p->inert.end(),inert);
 *ledger=p->cumulative;*time=p->time;*epoch=p->epoch;return OK;
}
}
extern "C" int fusion_c_source_state_snapshot(const fusion_source_state_v1* p,
 double* s,double* t,fusion_source_ledger_v1* ledger,double* time,uint64_t* epoch){
 return snapshot(p,s,t,ledger,nullptr,time,epoch,false);
}
extern "C" int fusion_c_source_state_snapshot_inert(const fusion_source_state_v1* p,
 double* s,double* t,fusion_source_ledger_v1* ledger,double* inert,double* time,uint64_t* epoch){
 return snapshot(p,s,t,ledger,inert,time,epoch,true);
}
extern "C" int fusion_c_source_state_begin(fusion_source_state_v1* p,double dt,uint64_t* ticket){
 if(!ticket)return PB11_STATUS_NULL_OUTPUT;
 *ticket=0;if(!p)return BAD;
 // Even an invalid replacement must not leave an older staged result committable.
 p->pending=false;p->staged=false;
 if(!std::isfinite(dt)||dt<=0)return BAD;
 double next=p->time+dt;
 if(!std::isfinite(next)||next<=p->time||p->counter==UINT64_MAX||p->epoch==UINT64_MAX)return NUM;
 p->trial_time=next;++p->counter;p->pending=true;*ticket=p->counter;return OK;
}
namespace {
int stage(fusion_source_state_v1* p,uint64_t ticket,
 const double* s,const double* t,const fusion_source_ledger_v1* step,
 const double* inert,bool extended_step){
 if(!p||!p->pending||ticket!=p->counter)return BAD;
 p->staged=false;
 try{
  if(!step||!ledger_valid(*step)||(extended_step&&!inert))return BAD;
  std::array<double,6> step_inert{},cumulative_inert{};
  for(int i=0;i<6;++i){
   if(inert)step_inert[i]=inert[i];
   if(!std::isfinite(step_inert[i]))return BAD;
   const long double amount=static_cast<long double>(p->inert[i])+step_inert[i];
   if(!representable(amount))return NUM;
   cumulative_inert[i]=static_cast<double>(amount);
  }
  Moments n{},u{},old_n{},old_u{};
  if(!inventory(p->edges,s,t,n,u)||!inventory(p->edges,p->s.data(),p->t.data(),old_n,old_u))return BAD;
  if(!inventory_balance(old_n,old_u,n,u,*step,step_inert))return NUM;
  fusion_source_ledger_v1 cumulative{};
  if(!add_ledger(p->cumulative,*step,cumulative)||
     !inventory_balance(p->initial_n,p->initial_u,n,u,cumulative,cumulative_inert))return NUM;
  // Allocations complete before the valid-stage flag changes.
  std::vector<double> ts(s,s+p->s.size()),tt(t,t+p->t.size());
  p->trial_s.swap(ts);p->trial_t.swap(tt);p->staged_cumulative=cumulative;
  p->staged_inert=cumulative_inert;p->staged_extended=p->extended||extended_step;
  p->staged=true;return OK;
 }catch(...){return EXC;}
}
}
extern "C" int fusion_c_source_state_stage(fusion_source_state_v1* p,uint64_t ticket,
 const double* s,const double* t,const fusion_source_ledger_v1* step){
 return stage(p,ticket,s,t,step,nullptr,false);
}
extern "C" int fusion_c_source_state_stage_inert(fusion_source_state_v1* p,uint64_t ticket,
 const double* s,const double* t,const fusion_source_ledger_v1* step,const double* inert){
 return stage(p,ticket,s,t,step,inert,true);
}
extern "C" int fusion_c_source_state_commit(fusion_source_state_v1* p,uint64_t ticket){
 if(!p||!p->pending||!p->staged||ticket!=p->counter)return BAD;
 p->s.swap(p->trial_s);p->t.swap(p->trial_t);p->cumulative=p->staged_cumulative;
 p->inert=p->staged_inert;p->extended=p->staged_extended;
 p->time=p->trial_time;++p->epoch;p->pending=false;p->staged=false;return OK;
}
extern "C" int fusion_c_source_state_discard(fusion_source_state_v1* p,uint64_t ticket){
 if(!p||!p->pending||ticket!=p->counter)return BAD;
 p->pending=false;p->staged=false;return OK;
}
namespace {
constexpr uint64_t magic=UINT64_C(0x3154535346554250); // PB UFSST1, opaque schema magic
constexpr size_t fixed_words=8+12+ledger_count+1; // header8, anchors12, ledger, checksum
size_t packed_size(size_t cells,bool extended=false){return 8*(fixed_words+13*cells+1+(extended?6:0));}
uint64_t checksum(const unsigned char* p,size_t n){uint64_t h=UINT64_C(14695981039346656037);
 for(size_t i=0;i<n;++i){h^=p[i];h*=UINT64_C(1099511628211);}return h;}
void write_u64(unsigned char*& p,uint64_t x){for(int i=0;i<8;++i){*p++=static_cast<unsigned char>(x&255);x>>=8;}}
uint64_t read_u64(const unsigned char*& p){uint64_t x=0;for(int i=0;i<8;++i)x|=uint64_t(*p++)<<(8*i);return x;}
void write_double(unsigned char*& p,double x){uint64_t bits;std::memcpy(&bits,&x,8);write_u64(p,bits);}
double read_double(const unsigned char*& p){uint64_t bits=read_u64(p);double x;std::memcpy(&x,&bits,8);return x;}
static_assert(sizeof(double)==8&&std::numeric_limits<double>::is_iec559,"restart requires IEEE754 binary64");
}
extern "C" int fusion_c_source_state_pack_size(const fusion_source_state_v1* p,size_t* required){
 if(!required)return PB11_STATUS_NULL_OUTPUT;
 *required=0;if(!p||p->pending)return BAD;
 *required=packed_size(p->edges.size()-1,p->extended);return OK;
}
extern "C" int fusion_c_source_state_pack(const fusion_source_state_v1* p,unsigned char* buffer,
 size_t capacity,size_t* written){
 if(!written)return PB11_STATUS_NULL_OUTPUT;
 *written=0;
 if(!p||p->pending||!buffer)return BAD;
 size_t size=packed_size(p->edges.size()-1,p->extended);
 if(capacity<size)return BAD;
 auto cursor=buffer;
 write_u64(cursor,magic);write_u64(cursor,p->extended?2:1);write_u64(cursor,p->edges.size()-1);
 write_u64(cursor,p->tag);write_u64(cursor,p->epoch);write_u64(cursor,p->counter);
 write_double(cursor,p->time);write_double(cursor,p->initial_time);
 for(auto v:p->initial_n)write_double(cursor,static_cast<double>(v));
 for(auto v:p->initial_u)write_double(cursor,static_cast<double>(v));
 fields(p->cumulative,[&](double v){write_double(cursor,v);});
 if(p->extended)for(double v:p->inert)write_double(cursor,v);
 for(auto v:p->edges)write_double(cursor,v);
 for(auto v:p->s)write_double(cursor,v);
 for(auto v:p->t)write_double(cursor,v);
 write_u64(cursor,checksum(buffer,size-8));*written=size;return OK;
}
extern "C" int fusion_c_source_state_unpack(const unsigned char* buffer,size_t length,
 uint64_t tag,fusion_source_state_v1** out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out=nullptr;
 try{
  if(!buffer||length<packed_size(1))return BAD;
  const unsigned char* cursor=buffer;
  if(read_u64(cursor)!=magic)return BAD;
  const uint64_t version=read_u64(cursor);
  if(version!=1&&version!=2)return BAD;
  uint64_t cells=read_u64(cursor);
  if(cells<1||cells>1000000||length!=packed_size(size_t(cells),version==2))return BAD;
  const unsigned char* end=buffer+length-8;
  if(read_u64(end)!=checksum(buffer,length-8))return BAD;
  auto p=std::make_unique<fusion_source_state_v1>();p->extended=version==2;p->tag=read_u64(cursor);
  if(p->tag!=tag)return BAD;
 p->epoch=read_u64(cursor);p->counter=read_u64(cursor);
  p->time=read_double(cursor);p->initial_time=read_double(cursor);
  if(!std::isfinite(p->time)||!std::isfinite(p->initial_time)||p->initial_time<0||
     p->time<p->initial_time||p->epoch>p->counter||
     (p->epoch==0&&p->time!=p->initial_time)||(p->epoch>0&&p->time<=p->initial_time))return BAD;
  for(auto& v:p->initial_n){v=read_double(cursor);if(!std::isfinite(v)||v<0)return BAD;}
  for(auto& v:p->initial_u){v=read_double(cursor);if(!std::isfinite(v)||v<0)return BAD;}
  fields(p->cumulative,[&](double& v){v=read_double(cursor);});
  if(!ledger_valid(p->cumulative))return BAD;
  if(p->extended){
   // Extended format can only arise from an accepted extended stage.
   if(p->epoch==0)return BAD;
   for(double& v:p->inert){v=read_double(cursor);if(!std::isfinite(v))return BAD;}
  }
  if(p->epoch==0){bool zero=true;fields(p->cumulative,[&](double v){if(v!=0)zero=false;});if(!zero)return BAD;}
  p->edges.resize(size_t(cells)+1);p->s.resize(6*size_t(cells));p->t.resize(p->s.size());
  for(size_t i=0;i<p->edges.size();++i){double v=read_double(cursor);p->edges[i]=v;
   if(!std::isfinite(v)||v<0||(i&&v<=p->edges[i-1]))return BAD;}
  for(auto& v:p->s)v=read_double(cursor);
  for(auto& v:p->t)v=read_double(cursor);
  Moments n{},u{};
  for(int j=0;j<6;++j){
   if(p->initial_n[j]==0&&p->initial_u[j]!=0)return BAD;
   const long double minimum=p->initial_n[j]*(static_cast<long double>(p->edges[0])+p->edges[1])/2;
   const long double maximum=p->initial_n[j]*(static_cast<long double>(p->edges[cells-1])+p->edges[cells])/2;
   if(p->initial_u[j]<minimum&&!balance({p->initial_u[j],-minimum}))return BAD;
   if(p->initial_u[j]>maximum&&!balance({p->initial_u[j],-maximum}))return BAD;
  }
  if(!inventory(p->edges,p->s.data(),p->t.data(),n,u)||
     !inventory_balance(p->initial_n,p->initial_u,n,u,p->cumulative,p->inert))return BAD;
  *out=p.release();return OK;
 }catch(...){return EXC;}
}
