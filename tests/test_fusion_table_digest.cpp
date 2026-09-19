#include "fusion_birth_table.h"
#include "fusion_beam_birth_table.h"
#include "../src/sha256_internal.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using Bytes=std::vector<unsigned char>;
using Digest=std::array<unsigned char,32>;
void check(bool ok,const char*label){if(!ok)throw std::runtime_error(label);}
std::string hex(const Digest&d){std::ostringstream s;for(auto b:d)s<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(b);return s.str();}
void known(const std::string&s,const char*expected){check(hex(fusion_detail::sha256(reinterpret_cast<const unsigned char*>(s.data()),s.size()))==expected,"SHA256 known answer");}
template<class T,class Size,class Pack,class Unpack,class Hash,class Destroy>
void table_check(T*t,const char*name,Size size_fn,Pack pack,Unpack unpack,Hash hash,Destroy destroy){
 size_t n=0,w=0,len=0;Digest d{},again{};
 check(size_fn(t,&n)==0,"pack size");Bytes bytes(n),after(n);
 check(pack(t,bytes.data(),n,&w)==0&&w==n,"pack");
 check(hash(t,d.data(),&len)==0&&len==n,"table digest");
 check(d==fusion_detail::sha256(bytes.data(),n),"hash entire canonical pack");
 check(pack(t,after.data(),n,&w)==0&&after==bytes,"digest leaves table unchanged");
 T*copy=nullptr;check(unpack(bytes.data(),n,&copy)==0&&copy,"unpack");
 check(hash(copy,again.data(),&len)==0&&len==n&&again==d,"unpacked identity");destroy(copy);
 again.fill(9);len=9;check(hash(nullptr,again.data(),&len)==PB11_STATUS_INVALID_ARGUMENT,"null table");
 check(len==0&&std::all_of(again.begin(),again.end(),[](auto b){return b==0;}),"failure clears outputs");
 again.fill(9);check(hash(t,again.data(),nullptr)==PB11_STATUS_NULL_OUTPUT,"null length");
 check(std::all_of(again.begin(),again.end(),[](auto b){return b==0;}),"null length clears digest");
 len=9;check(hash(t,nullptr,&len)==PB11_STATUS_NULL_OUTPUT&&len==0,"null digest clears length");
 check(hash(t,nullptr,nullptr)==PB11_STATUS_NULL_OUTPUT,"both outputs null");
 std::ofstream f(std::string(name)+"-digest.bin",std::ios::binary);f.write(reinterpret_cast<const char*>(bytes.data()),n);check(bool(f),"write test bytes");
 std::cout<<name<<' '<<n<<' '<<hex(d)<<'\n';
}
int main(int argc,char**argv){try{
 if(argc==3&&std::string(argv[1])=="--hash"){
  std::ifstream f(argv[2],std::ios::binary);check(bool(f),"hash input");Bytes b((std::istreambuf_iterator<char>(f)),{});
  std::cout<<hex(fusion_detail::sha256(b.data(),b.size()))<<'\n';return 0;
 }
 known("","e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
 known("abc","ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
 known("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq","248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
 known(std::string(1000000,'a'),"cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
 constexpr double kev=1.602176634e-16,mev=1.602176634e-13;constexpr int cells=32;
 double edges[cells+1];for(int j=0;j<=cells;++j)edges[j]=25*mev*j/cells;
 fusion_thermal_birth_options_v1 source{5*mev,40,91.84*kev,.001*mev,.76,0,.051,1,1,0,0,13,0,16,12,4,4};
 fusion_birth_table_control_v1 control{.03,.03,.03,.03,1e-4,1e-4,128,1024,12};
 fusion_birth_table_v1*t=nullptr;
 check(fusion_c_birth_table_create(3,9*kev,11*kev,&source,&control,cells,edges,&t)==0,"thermal create");
 table_check(t,"thermal",fusion_c_birth_table_pack_size,fusion_c_birth_table_pack,fusion_c_birth_table_unpack,fusion_c_birth_table_content_digest,fusion_c_birth_table_destroy);
 Digest first{},second{};size_t n=0;check(fusion_c_birth_table_content_digest(t,first.data(),&n)==0,"first identity");fusion_c_birth_table_destroy(t);
 control.max_rate_error=.04;
 check(fusion_c_birth_table_create(3,9*kev,11*kev,&source,&control,cells,edges,&t)==0,"changed thermal create");
 check(fusion_c_birth_table_content_digest(t,second.data(),&n)==0&&first!=second,"changed table metadata changes identity");fusion_c_birth_table_destroy(t);
 fusion_beam_birth_options_v1 beam{};beam.relative_max_J=2.5*mev;beam.angular_max_exponent=40;
 beam.ground_state_q_J=91.84*kev;beam.cutoff_J=.001*mev;beam.l1_fraction=.76;beam.narrow_peak_fraction=.051;
 beam.continuum_peak_scale=1;beam.continuation=1;beam.broad_mode=13;beam.relative_order=16;beam.angular_order=4;beam.nq=4;beam.ncos=4;
 fusion_beam_birth_table_v1*b=nullptr;
 check(fusion_c_beam_birth_table_create(3,0,100*kev,9*kev,11*kev,&beam,&control,cells,edges,&b)==0,"beam create");
 table_check(b,"beam",fusion_c_beam_birth_table_pack_size,fusion_c_beam_birth_table_pack,fusion_c_beam_birth_table_unpack,fusion_c_beam_birth_table_content_digest,fusion_c_beam_birth_table_destroy);
 fusion_c_beam_birth_table_destroy(b);std::cout<<"TABLE_DIGEST_PASS\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
