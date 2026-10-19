// Standalone: link with --wrap=fusion_c_birth_table_destroy and
// --wrap=fusion_c_beam_birth_table_destroy. Supply THERMAL_PACKED BEAM_PACKED.
#include "fusion_capture_registry.h"
#include <array>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>
static unsigned thermal_destroys=0,beam_destroys=0;
static bool reject_content_digest=false;
extern "C" int __real_fusion_c_birth_table_content_digest(const fusion_birth_table_v1*,unsigned char*,size_t*);
extern "C" int __real_fusion_c_beam_birth_table_content_digest(const fusion_beam_birth_table_v1*,unsigned char*,size_t*);
extern "C" int __wrap_fusion_c_birth_table_content_digest(const fusion_birth_table_v1* p,unsigned char* d,size_t* n) {
 if(reject_content_digest){std::memset(d,0,32);*n=0;return 3;}
 return __real_fusion_c_birth_table_content_digest(p,d,n);
}
extern "C" int __wrap_fusion_c_beam_birth_table_content_digest(const fusion_beam_birth_table_v1* p,unsigned char* d,size_t* n) {
 if(reject_content_digest){std::memset(d,0,32);*n=0;return 3;}
 return __real_fusion_c_beam_birth_table_content_digest(p,d,n);
}
extern "C" void __real_fusion_c_birth_table_destroy(fusion_birth_table_v1 *);
extern "C" void __real_fusion_c_beam_birth_table_destroy(fusion_beam_birth_table_v1 *);
extern "C" void __wrap_fusion_c_birth_table_destroy(fusion_birth_table_v1 *p) {
 ++thermal_destroys;__real_fusion_c_birth_table_destroy(p);
}
extern "C" void __wrap_fusion_c_beam_birth_table_destroy(fusion_beam_birth_table_v1 *p) {
 ++beam_destroys;__real_fusion_c_beam_birth_table_destroy(p);
}
static void require(bool ok,const char *message) {if(!ok)throw std::runtime_error(message);}
static std::vector<char> read(const char *path) {
 std::ifstream f(path,std::ios::binary);require(bool(f),"open fixture");
 return {std::istreambuf_iterator<char>(f),std::istreambuf_iterator<char>()};
}
template<class T,class Unpack,class Adopt,class Destroy,class RawDestroy,class Digest,class Pack>
static void exercise(const char *path,int kind,unsigned &released,Unpack unpack,
 Adopt adopt,Destroy destroy,RawDestroy raw_destroy,Digest digest,Pack pack) {
 const auto bytes=read(path);require(!bytes.empty(),"empty fixture");
 auto raw=[&](){T *p=nullptr;require(unpack(bytes.data(),bytes.size(),&p)==0&&p,"core unpack");return p;};
 auto unchanged=[&](T *p){std::vector<char> packed(bytes.size());size_t n=0;
  require(pack(p,packed.data(),packed.size(),&n)==0&&n==bytes.size()&&packed==bytes,"raw remains packable");};
 fusion_capture_context_v1 *ctx=nullptr;require(fusion_capture_context_create_v1(4,&ctx)==0,"context");
 T *p=raw();T *same=p;std::array<unsigned char,32> before{},after{};size_t n=0,m=0;
 require(digest(p,before.data(),&n)==0,"initial digest");int status=77;
 unsigned base=released;
 require(adopt(ctx,p,bytes.data(),bytes.size(),path,&status)==0&&status==0&&p==same,"adoption identity");
 require(digest(p,after.data(),&m)==0&&before==after&&n==m,"digest preserved");
 fusion_capture_table_identity_v1 identity{};
 require(fusion_capture_table_identity_v1_get(ctx,kind,p,&identity)==0,"registered identity");
 require(identity.packed_bytes==bytes.size()&&std::string(identity.qualified_path)==path,"metadata");
 constexpr char hex[]="0123456789abcdef";std::string expected;
 for(auto b:before){expected+=hex[b>>4];expected+=hex[b&15];}
 require(expected==identity.content_sha256,"registered digest");
 status=77;require(adopt(ctx,p,bytes.data(),bytes.size(),path,&status)==100&&status==-1&&released==base,"duplicate preserves ownership");
 // Wrong-kind destroy searches the registry before dereferencing the opaque
 // pointer. Never pass a wrong typed pointer into adoption/core operations.
 if(kind==FUSION_CAPTURE_THERMAL_TABLE){auto wrong=reinterpret_cast<fusion_beam_birth_table_v1 *>(p);
  require(fusion_capture_beam_destroy_v1(ctx,&wrong)==103&&wrong==reinterpret_cast<fusion_beam_birth_table_v1 *>(p),"wrong kind destroy");}
 else {auto wrong=reinterpret_cast<fusion_birth_table_v1 *>(p);
  require(fusion_capture_thermal_destroy_v1(ctx,&wrong)==103&&wrong==reinterpret_cast<fusion_birth_table_v1 *>(p),"wrong kind destroy");}
 require(released==base,"wrong kind retains");
 require(destroy(ctx,&p)==0&&!p&&released==base+1,"registry destroy clears");
 p=raw();require(adopt(ctx,p,bytes.data(),bytes.size(),path,&status)==0,"context adoption");
 fusion_capture_context_destroy_v1(ctx);require(released==base+2,"context destroys exactly once");
 require(fusion_capture_context_create_v1(2,&ctx)==0,"failure context");
 T *anchor=raw();require(adopt(ctx,anchor,bytes.data(),bytes.size(),path,&status)==0,"anchor");
 p=raw();base=released;
 auto fail=[&](fusion_capture_context_v1 *c,const void *b,size_t len,const char *name,int *out,bool called){
  if(out)*out=77;
  require(adopt(c,p,b,len,name,out)==100,"invalid diagnostic");
  if(out)require(*out==(called?0:-1),"core status separation");
  require(released==base,"failure must not destroy");unchanged(p);
  fusion_capture_table_identity_v1 now{};
  require(fusion_capture_table_identity_v1_get(ctx,kind,anchor,&now)==0&&
   std::memcmp(&identity,&now,sizeof now)==0,"existing entry unchanged");
  require(fusion_capture_table_identity_v1_get(ctx,kind,p,&now)==103,"failure unregistered");
 };
 auto corrupt=bytes;corrupt.back()^=1;
 fail(ctx,corrupt.data(),corrupt.size(),path,&status,true);
 require(corrupt.back()==char(bytes.back()^1),"input unchanged");
 fail(ctx,bytes.data(),bytes.size()-1,path,&status,true);
 fail(ctx,bytes.data(),0,path,&status,false);
 fail(ctx,nullptr,bytes.size(),path,&status,false);
 fail(nullptr,bytes.data(),bytes.size(),path,&status,false);
 fail(ctx,bytes.data(),bytes.size(),nullptr,&status,false);
 fail(ctx,bytes.data(),bytes.size(),"relative",&status,false);
 fail(ctx,bytes.data(),bytes.size(),"/",&status,false);
 const std::string longpath="/"+std::string(4096,'x');
 fail(ctx,bytes.data(),bytes.size(),longpath.c_str(),&status,false);
 fail(ctx,bytes.data(),bytes.size(),path,nullptr,false);
 fail(ctx,bytes.data(),(kind==1?512ULL:256ULL)*1024*1024+1,path,&status,false);
 status=77;require(adopt(ctx,nullptr,bytes.data(),bytes.size(),path,&status)==100&&status==-1,"null table");
 status=77;reject_content_digest=true;
 const int rejected=adopt(ctx,p,bytes.data(),bytes.size(),path,&status);
 reject_content_digest=false;
 require(rejected==100&&status==3&&released==base,"simulated serializer status preserved; no destruction");
 unchanged(p);
 fusion_capture_table_identity_v1 absent{};
 require(fusion_capture_table_identity_v1_get(ctx,kind,p,&absent)==103,"serializer failure not registered");
 raw_destroy(p);require(released==base+1,"caller destroys rejected raw");
 fusion_capture_context_destroy_v1(ctx);require(released==base+2,"anchor retained");
 require(fusion_capture_context_create_v1(1,&ctx)==0,"capacity context");
 anchor=raw();require(adopt(ctx,anchor,bytes.data(),bytes.size(),path,&status)==0,"fill capacity");
 p=raw();base=released;status=77;
 require(adopt(ctx,p,bytes.data(),bytes.size(),path,&status)==102&&status==-1&&released==base,"capacity rejection");
 unchanged(p);raw_destroy(p);fusion_capture_context_destroy_v1(ctx);
 require(released==base+2,"capacity ownership");
 require(read(path)==bytes,"fixture file unchanged");
}
int main(int argc,char **argv) {
 try {
  require(argc==3,"usage: test_created_adoption THERMAL_PACKED BEAM_PACKED (absolute paths)");
  exercise<fusion_birth_table_v1>(argv[1],1,thermal_destroys,
   fusion_c_birth_table_unpack,fusion_capture_thermal_adopt_v1,fusion_capture_thermal_destroy_v1,
   fusion_c_birth_table_destroy,fusion_c_birth_table_content_digest,fusion_c_birth_table_pack);
  exercise<fusion_beam_birth_table_v1>(argv[2],2,beam_destroys,
   fusion_c_beam_birth_table_unpack,fusion_capture_beam_adopt_v1,fusion_capture_beam_destroy_v1,
   fusion_c_beam_birth_table_destroy,fusion_c_beam_birth_table_content_digest,fusion_c_beam_birth_table_pack);
  std::cout<<"created adoption checks passed\n";return 0;
 }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
