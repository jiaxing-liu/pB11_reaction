#include "fusion_birth_table.h"
#include "fusion_beam_birth_table.h"
#include "sha256_internal.hpp"
#include <algorithm>
#include <vector>

namespace {
template<class Table, class Size, class Pack>
int digest_table(const Table* table, unsigned char* digest, size_t* length,
                 size_t cap, Size size_fn, Pack pack_fn) noexcept {
 if(digest)std::fill_n(digest,32,0);
 if(length)*length=0;
 if(!digest||!length)return PB11_STATUS_NULL_OUTPUT;
 if(!table)return PB11_STATUS_INVALID_ARGUMENT;
 try {
  size_t needed=0,written=0;
  int status=size_fn(table,&needed);
  if(status)return status;
  if(!needed||needed>cap)return PB11_STATUS_INVALID_ARGUMENT;
  std::vector<unsigned char> bytes(needed);
  status=pack_fn(table,bytes.data(),bytes.size(),&written);
  if(status)return status;
  if(written!=needed)return PB11_STATUS_INVALID_ARGUMENT;
  const auto result=fusion_detail::sha256(bytes.data(),bytes.size());
  std::copy(result.begin(),result.end(),digest);*length=written;
  return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
}
extern "C" int fusion_c_birth_table_content_digest(
 const fusion_birth_table_v1* table,unsigned char* digest,size_t* length){
 return digest_table(table,digest,length,512ULL*1024*1024,
                     fusion_c_birth_table_pack_size,fusion_c_birth_table_pack);
}
extern "C" int fusion_c_beam_birth_table_content_digest(
 const fusion_beam_birth_table_v1* table,unsigned char* digest,size_t* length){
 return digest_table(table,digest,length,256ULL*1024*1024,
                     fusion_c_beam_birth_table_pack_size,fusion_c_beam_birth_table_pack);
}
