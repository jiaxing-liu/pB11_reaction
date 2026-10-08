#pragma once
#include "fusion_capture_trial.h"
#include "sha256_stream.hpp"
#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

namespace pb11_diagnostics { namespace detail {
static_assert(CHAR_BIT==8&&sizeof(double)==8&&std::numeric_limits<double>::is_iec559,
 "Trial records require IEEE754 binary64 and 8-bit bytes");
static_assert(sizeof(int)==4,"Trial C ABI requires 32-bit int");
constexpr uint64_t maximum_record=128ULL*1024*1024;
struct Error:std::runtime_error {
 int status;
 Error(int code,const char *message):std::runtime_error(message),status(code){}
};
inline void need(bool ok,int status,const char *message) {
 if(!ok)throw Error(status,message);
}
inline std::string hex(const std::array<unsigned char,32> &digest) {
 constexpr char digits[]="0123456789abcdef";std::string out(64,'0');
 for(size_t i=0;i<digest.size();++i){out[2*i]=digits[digest[i]>>4];out[2*i+1]=digits[digest[i]&15];}
 return out;
}
inline bool valid_hex(const std::string &value) {
 if(value.size()!=64)return false;
 for(char c:value)if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false;
 return true;
}
inline std::array<unsigned char,8> word(uint64_t value) {
 std::array<unsigned char,8> bytes{};
 for(unsigned i=0;i<8;++i)bytes[i]=static_cast<unsigned char>(value>>(8*i));
 return bytes;
}
inline uint64_t bits(double value) {uint64_t out;std::memcpy(&out,&value,8);return out;}
inline double number(uint64_t value) {double out;std::memcpy(&out,&value,8);return out;}
// Hash-only output sink uses the same explicit word encoding as the wire file.
struct HashSink {
 Sha256 hash;
 void raw(const void *p,size_t n){hash.update(p,n);}
 void u64(uint64_t v){auto b=word(v);raw(b.data(),8);}
 void u32(uint32_t v){auto b=word(v);raw(b.data(),4);}
 void i32(int v){u32(static_cast<uint32_t>(v));}
 void f64(double v){u64(bits(v));}
 void doubles(const double *v,size_t n){for(size_t i=0;i<n;++i)f64(v[i]);}
};
class Writer {
 std::string temporary_,target_;
 int fd_=-1;
 uint64_t limit_,position_=0,field_end_=0;
 Sha256 hash_;
public:
 Writer(const char *path,uint64_t maximum):limit_(maximum) {
  need(path&&*path&&maximum>0&&maximum<=maximum_record,FUSION_CAPTURE_INVALID,"invalid record path/budget");
  target_=path;
#ifdef _WIN32
  throw Error(FUSION_CAPTURE_UNSUPPORTED,"exclusive atomic writer requires POSIX");
#else
  auto target=std::filesystem::path(path);
  auto parent=target.parent_path();if(parent.empty())parent=".";
  temporary_=(parent/("."+target.filename().string()+".capture.XXXXXX")).string();
  std::vector<char> name(temporary_.begin(),temporary_.end());name.push_back(0);
  fd_=::mkstemp(name.data());
  need(fd_>=0,FUSION_CAPTURE_IO,"cannot create exclusive temporary record");
  std::memcpy(&temporary_[0],name.data(),temporary_.size());
#endif
 }
 ~Writer() {
#ifndef _WIN32
  if(fd_>=0)::close(fd_);
  if(!temporary_.empty())::unlink(temporary_.c_str());
#endif
 }
 Writer(const Writer&)=delete;Writer& operator=(const Writer&)=delete;
 void raw(const void *data,size_t size,bool hashed=true) {
  need(size<=limit_-position_,FUSION_CAPTURE_CAPACITY,"record budget exceeded");
  if(field_end_)need(size<=field_end_-position_,FUSION_CAPTURE_FORMAT,"writer field length exceeded");
#ifndef _WIN32
  auto *p=static_cast<const unsigned char *>(data);size_t done=0;
  while(done<size){auto n=::write(fd_,p+done,size-done);if(n<0&&errno==EINTR)continue;
   need(n>0,FUSION_CAPTURE_IO,"record write failed");done+=static_cast<size_t>(n);}
#endif
  if(hashed)hash_.update(data,size);
  position_+=size;
 }
 void u64(uint64_t v){auto b=word(v);raw(b.data(),8);}
 void u32(uint32_t v){auto b=word(v);raw(b.data(),4);}
 void i32(int v){u32(static_cast<uint32_t>(v));}
 void f64(double v){u64(bits(v));}
 void doubles(const double *v,size_t n) {
  std::array<unsigned char,8192> chunk{};size_t consumed=0;
  while(consumed<n){size_t count=std::min(n-consumed,chunk.size()/8);
   for(size_t j=0;j<count;++j){auto b=word(bits(v[consumed+j]));std::memcpy(chunk.data()+8*j,b.data(),8);}
   raw(chunk.data(),8*count);consumed+=count;}
 }
 void field(uint32_t tag,uint64_t size) {
  need(field_end_==0,FUSION_CAPTURE_FORMAT,"unclosed writer field");
  u32(tag);u64(size);
  need(size<=limit_-position_,FUSION_CAPTURE_CAPACITY,"field exceeds record budget");
  field_end_=position_+size;
 }
 void end(){need(position_==field_end_,FUSION_CAPTURE_FORMAT,"writer field length mismatch");field_end_=0;}
 void publish() {
  need(field_end_==0,FUSION_CAPTURE_FORMAT,"incomplete record field");
  auto digest=hash_.digest();raw(digest.data(),digest.size(),false);
#ifndef _WIN32
  need(::fsync(fd_)==0,FUSION_CAPTURE_IO,"record flush failed");
  int fd=fd_;fd_=-1;need(::close(fd)==0,FUSION_CAPTURE_IO,"record close failed");
  need(::link(temporary_.c_str(),target_.c_str())==0,FUSION_CAPTURE_IO,"cannot publish exclusive record");
  ::unlink(temporary_.c_str());temporary_.clear();
#endif
 }
};
class Reader {
 std::ifstream file_;
 uint64_t size_,position_=0,field_end_=0;
 Sha256 hash_;
public:
 Reader(const char *path,uint64_t maximum) {
  need(path&&*path&&maximum>0&&maximum<=maximum_record,FUSION_CAPTURE_INVALID,"invalid read path/budget");
  auto status=std::filesystem::symlink_status(path);
  need(std::filesystem::is_regular_file(status),FUSION_CAPTURE_IO,"regular nonsymlink record required");
  size_=std::filesystem::file_size(path);
  need(size_>=44&&size_<=maximum,FUSION_CAPTURE_CAPACITY,"record size outside budget");
  file_.open(path,std::ios::binary);need(bool(file_),FUSION_CAPTURE_IO,"cannot open record");
 }
 void raw(void *data,size_t n,bool hashed=true) {
  need(n<=size_-position_,FUSION_CAPTURE_FORMAT,"truncated record");
  if(field_end_)need(n<=field_end_-position_,FUSION_CAPTURE_FORMAT,"reader field length exceeded");
  need(bool(file_.read(static_cast<char *>(data),static_cast<std::streamsize>(n))),FUSION_CAPTURE_FORMAT,"truncated bytes");
  if(hashed)hash_.update(data,n);
  position_+=n;
 }
 uint64_t u64(){std::array<unsigned char,8>b{};raw(b.data(),8);uint64_t v=0;for(unsigned i=0;i<8;++i)v|=uint64_t(b[i])<<(8*i);return v;}
 uint32_t u32(){std::array<unsigned char,4>b{};raw(b.data(),4);uint32_t v=0;for(unsigned i=0;i<4;++i)v|=uint32_t(b[i])<<(8*i);return v;}
 int i32(){uint32_t v=u32();return v<=INT_MAX?static_cast<int>(v):static_cast<int>(static_cast<int64_t>(v)-4294967296LL);}
 double f64(){return number(u64());}
 bool flag(){auto value=u32();need(value<=1,FUSION_CAPTURE_FORMAT,"invalid presence flag");return value!=0;}
 void field(uint32_t expected,uint64_t exact=UINT64_MAX) {
  need(field_end_==0,FUSION_CAPTURE_FORMAT,"unclosed reader field");
  auto tag=u32();auto length=u64();
  need(tag==expected,FUSION_CAPTURE_FORMAT,"unknown/duplicate/reordered field");
  need(size_-position_>=32&&length<=size_-position_-32,FUSION_CAPTURE_FORMAT,"field exceeds payload");
  if(exact!=UINT64_MAX)need(length==exact,FUSION_CAPTURE_FORMAT,"invalid field shape");
  field_end_=position_+length;
 }
 void end(){need(position_==field_end_,FUSION_CAPTURE_FORMAT,"reader field length mismatch");field_end_=0;}
 uint64_t remaining() const {need(field_end_>=position_&&field_end_!=0,FUSION_CAPTURE_FORMAT,"no reader field");return field_end_-position_;}
 std::vector<double> doubles(size_t n) {
  need(field_end_&&n<=(field_end_-position_)/8,FUSION_CAPTURE_FORMAT,"invalid array extent");
  std::vector<double> v(n);for(auto &x:v)x=f64();return v;
 }
 std::string text(size_t n){std::string value(n,'\0');if(n)raw(&value[0],n);return value;}
 void finish() {
  need(field_end_==0&&position_+32==size_,FUSION_CAPTURE_FORMAT,"incomplete/trailing record");
  auto expected=hash_.digest();std::array<unsigned char,32>actual{};raw(actual.data(),32,false);
  need(actual==expected,FUSION_CAPTURE_FORMAT,"snapshot checksum differs");
  char extra;need(!file_.get(extra)&&file_.eof(),FUSION_CAPTURE_FORMAT,"file changed or has trailing bytes");
 }
};
}} // namespace
