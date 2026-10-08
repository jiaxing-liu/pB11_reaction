#include "fusion_capture_registry.h"
#include "detail/sha256_stream.hpp"
#include <algorithm>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <stdexcept>
#include <string>

namespace {
struct Entry { int kind; uint64_t length; std::string path, digest, kernel; };
void release(int kind,const void *handle) noexcept {
 if(kind==FUSION_CAPTURE_THERMAL_TABLE)
  fusion_c_birth_table_destroy(const_cast<fusion_birth_table_v1 *>(
   static_cast<const fusion_birth_table_v1 *>(handle)));
 else fusion_c_beam_birth_table_destroy(const_cast<fusion_beam_birth_table_v1 *>(
   static_cast<const fusion_beam_birth_table_v1 *>(handle)));
}
std::string hex(const std::array<unsigned char,32> &bytes) {
 constexpr char digits[]="0123456789abcdef";
 std::string out(64,'0');
 for(size_t i=0;i<bytes.size();++i) {
  out[2*i]=digits[bytes[i]>>4];out[2*i+1]=digits[bytes[i]&15];
 }
 return out;
}
bool path_valid(const char *path) {
 if(!path||path[0]!='/')return false;
 size_t n=0;while(n<=4096&&path[n])++n;
 return n>1&&n<=4096;
}
}
struct fusion_capture_context_v1 {
 uint64_t capacity;
 std::mutex mutex;
 std::map<const void *,Entry> entries;
 explicit fusion_capture_context_v1(uint64_t limit):capacity(limit){}
 ~fusion_capture_context_v1() {
  for(const auto &entry:entries)release(entry.second.kind,entry.first);
 }
};
namespace {
template<typename T,typename Unpack>
int unpack(fusion_capture_context_v1 *context,const void *bytes,size_t length,
 const char *path,T **out,int *unpack_status,int kind,Unpack original) noexcept {
 if(out)*out=nullptr;
 if(unpack_status)*unpack_status=-1;
 if(!context||!out||!unpack_status||!bytes||!length||!path_valid(path)||
  length>(kind==FUSION_CAPTURE_THERMAL_TABLE?512ULL:256ULL)*1024*1024)
  return FUSION_CAPTURE_INVALID;
 T *table=nullptr;
 try {
  std::lock_guard<std::mutex> lock(context->mutex);
  if(context->entries.size()>=context->capacity)return FUSION_CAPTURE_CAPACITY;
  // Capture a stable reference from the caller's loader buffer. The original
  // unpack imports this same buffer, without an additional pack allocation.
  pb11_diagnostics::detail::Sha256 digest;
  digest.update(bytes,length);
  Entry entry{kind,static_cast<uint64_t>(length),path,hex(digest.digest()),
   fusion_c_beam_birth_table_kernel_identity()};
  *unpack_status=original(bytes,length,&table);
  if(*unpack_status!=0) {
   if(table)release(kind,table);
   return FUSION_CAPTURE_UNPACK_REJECTED;
  }
  if(!table)return FUSION_CAPTURE_UNPACK_REJECTED;
  if(!context->entries.emplace(table,std::move(entry)).second) {
   release(kind,table);return FUSION_CAPTURE_UNKNOWN_HANDLE;
  }
  *out=table;
  return FUSION_CAPTURE_OK;
 } catch(const std::bad_alloc &) {
  if(table)release(kind,table);
  return FUSION_CAPTURE_ALLOCATION;
 } catch(...) {
  if(table)release(kind,table);
  return FUSION_CAPTURE_INVALID;
 }
}
template<typename T>
int destroy(fusion_capture_context_v1 *context,T **handle,int kind) noexcept {
 if(!context||!handle||!*handle)return FUSION_CAPTURE_INVALID;
 try {
  std::lock_guard<std::mutex> lock(context->mutex);
  auto entry=context->entries.find(*handle);
  if(entry==context->entries.end()||entry->second.kind!=kind)
   return FUSION_CAPTURE_UNKNOWN_HANDLE;
  release(kind,*handle);context->entries.erase(entry);*handle=nullptr;
  return FUSION_CAPTURE_OK;
 } catch(...) {return FUSION_CAPTURE_INVALID;}
}
}
extern "C" int fusion_capture_context_create_v1(uint64_t capacity,
 fusion_capture_context_v1 **out) {
 if(out)*out=nullptr;
 if(!out||capacity==0||capacity>1000000)return FUSION_CAPTURE_INVALID;
 try {*out=new fusion_capture_context_v1(capacity);return FUSION_CAPTURE_OK;}
 catch(const std::bad_alloc &){return FUSION_CAPTURE_ALLOCATION;}
 catch(...){return FUSION_CAPTURE_INVALID;}
}
extern "C" void fusion_capture_context_destroy_v1(fusion_capture_context_v1 *context) {
 delete context;
}
extern "C" int fusion_capture_thermal_unpack_v1(fusion_capture_context_v1 *context,
 const void *bytes,size_t length,const char *path,fusion_birth_table_v1 **out,int *status) {
 return unpack(context,bytes,length,path,out,status,FUSION_CAPTURE_THERMAL_TABLE,
  fusion_c_birth_table_unpack);
}
extern "C" int fusion_capture_beam_unpack_v1(fusion_capture_context_v1 *context,
 const void *bytes,size_t length,const char *path,fusion_beam_birth_table_v1 **out,int *status) {
 return unpack(context,bytes,length,path,out,status,FUSION_CAPTURE_BEAM_TABLE,
  fusion_c_beam_birth_table_unpack);
}
extern "C" int fusion_capture_table_identity_v1_get(fusion_capture_context_v1 *context,
 int kind,const void *handle,fusion_capture_table_identity_v1 *out) {
 if(out)std::memset(out,0,sizeof(*out));
 if(!context||!handle||!out)return FUSION_CAPTURE_INVALID;
 try {
  std::lock_guard<std::mutex> lock(context->mutex);
  auto found=context->entries.find(handle);
  if(found==context->entries.end()||found->second.kind!=kind)
   return FUSION_CAPTURE_UNKNOWN_HANDLE;
  const auto &entry=found->second;
  if(entry.kernel.size()!=64)return FUSION_CAPTURE_INVALID;
  out->kind=kind;out->packed_bytes=entry.length;
  std::memcpy(out->content_sha256,entry.digest.data(),64);
  std::memcpy(out->kernel_identity,entry.kernel.data(),64);
  std::memcpy(out->qualified_path,entry.path.data(),entry.path.size());
  return FUSION_CAPTURE_OK;
 } catch(...){return FUSION_CAPTURE_INVALID;}
}
extern "C" int fusion_capture_thermal_destroy_v1(fusion_capture_context_v1 *context,
 fusion_birth_table_v1 **handle) {
 return destroy(context,handle,FUSION_CAPTURE_THERMAL_TABLE);
}
extern "C" int fusion_capture_beam_destroy_v1(fusion_capture_context_v1 *context,
 fusion_beam_birth_table_v1 **handle) {
 return destroy(context,handle,FUSION_CAPTURE_BEAM_TABLE);
}
