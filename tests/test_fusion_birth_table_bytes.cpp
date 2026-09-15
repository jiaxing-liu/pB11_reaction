// Independent round-trip and malformed-buffer tests. Reuse only the existing
// physical table fixture; the original test main remains a separate CTest.
#define main thermal_table_original_test_main
#include "test_fusion_birth_table.cpp"
#undef main
#include <cstring>
#include <memory>

namespace {
using BytesTableOwner=std::unique_ptr<fusion_birth_table_v1,decltype(&fusion_c_birth_table_destroy)>;
std::vector<unsigned char> pack_table(fusion_birth_table_v1* table){
 size_t n=123,written=123;
 expect_ok(fusion_c_birth_table_pack_size(table,&n),"size");check(n>96,"nonempty format");
 std::vector<unsigned char>b(n);
 expect_ok(fusion_c_birth_table_pack(table,b.data(),b.size(),&written),"pack");check(written==n,"exact size");return b;
}
uint64_t checksum(const unsigned char*p,size_t n){
 uint64_t h=14695981039346656037ULL;
 for(size_t i=0;i<n;++i){h^=p[i];h*=1099511628211ULL;}return h;
}
void put_word(std::vector<unsigned char>&b,size_t pos,uint64_t value){
 check(pos<=b.size()&&b.size()-pos>=8,"mutation offset");
 for(int i=0;i<8;++i)b[pos+i]=static_cast<unsigned char>(value>>(8*i));
}
void put_real(std::vector<unsigned char>&b,size_t pos,double value){
 uint64_t bits;std::memcpy(&bits,&value,8);put_word(b,pos,bits);
}
void fix_checksum(std::vector<unsigned char>&b){put_word(b,b.size()-8,checksum(b.data(),b.size()-8));}
void reject_buffer(const std::vector<unsigned char>& b,size_t n){
 auto*p=reinterpret_cast<fusion_birth_table_v1*>(uintptr_t(1));
 expect_rejected(fusion_c_birth_table_unpack(b.data(),n,&p),"invalid buffer");check(p==nullptr,"atomic null on failure");
}
}
int main(){try{
 auto source=source_options();auto control=table_control();auto edges=wide_energy_edges();
 const int channel=FUSION_DT_ALPHAN;const double lo=9*kKeVJ,hi=11*kKeVJ;
 fusion_birth_table_v1*p=nullptr;
 expect_ok(fusion_c_birth_table_create(channel,lo,hi,&source,&control,kCells,edges.data(),&p),"real DT fixture");
 BytesTableOwner original(p,fusion_c_birth_table_destroy);auto bytes=pack_table(p);
 p=nullptr;expect_ok(fusion_c_birth_table_unpack(bytes.data(),bytes.size(),&p),"unpack");
 BytesTableOwner loaded(p,fusion_c_birth_table_destroy);check(pack_table(p)==bytes,"byte-identical repack");
 for(double T:{lo,9.3*kKeVJ,10*kKeVJ,10.7*kKeVJ,hi}){
  std::vector<double>a(7*kCells),b(a.size());fusion_birth_coefficients_v1 ca{},cb{};
  expect_ok(fusion_c_birth_table_evaluate(original.get(),T,kCells,a.data(),&ca),"before");
  expect_ok(fusion_c_birth_table_evaluate(loaded.get(),T,kCells,b.data(),&cb),"after");
  check(std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0,"full spectrum bits");
  // Coefficients contain only double fields: no struct padding comparison.
  check(ca.reactivity_m3_s==cb.reactivity_m3_s,"rate");
  for(int i=0;i<2;++i)check(ca.reactant_energy_moment_J_m3_s[i]==cb.reactant_energy_moment_J_m3_s[i],"reactant debit");
  for(int i=0;i<7;++i){
   check(ca.below_number_m3_s[i]==cb.below_number_m3_s[i]&&ca.above_number_m3_s[i]==cb.above_number_m3_s[i],"number tails");
   check(ca.below_energy_J_m3_s[i]==cb.below_energy_J_m3_s[i]&&ca.above_energy_J_m3_s[i]==cb.above_energy_J_m3_s[i],"energy tails");
  }
 }
 int matches=-1;
 auto match=[&](double lower,const fusion_thermal_birth_options_v1&s,const fusion_birth_table_control_v1&c,const std::vector<double>&e){
  matches=-1;expect_ok(fusion_c_birth_table_matches_request(loaded.get(),channel,lower,hi,&s,&c,kCells,e.data(),&matches),"request compare");return matches;
 };
 check(match(lo,source,control,edges)==1,"full exact request");
 check(match(std::nextafter(lo,0.),source,control,edges)==0,"domain mismatch");
 auto altered=source;altered.cm_max_kT=39;check(match(lo,altered,control,edges)==0,"source mismatch");
 auto gate=control;gate.max_rate_error*=.9;check(match(lo,source,gate,edges)==0,"gate mismatch");
 auto grid=edges;grid[1]=std::nextafter(grid[1],grid[2]);check(match(lo,source,control,grid)==0,"edge mismatch");
 size_t count=9;expect_rejected(fusion_c_birth_table_pack_size(nullptr,&count),"null table");check(count==0,"size reset");
 std::vector<unsigned char>small(bytes.size()-1,0xa5);count=9;
 expect_rejected(fusion_c_birth_table_pack(original.get(),small.data(),small.size(),&count),"capacity");check(count==0,"written reset");
 for(size_t n:{size_t(0),size_t(95),bytes.size()-1})reject_buffer(bytes,n);
 auto bad=bytes;bad[24]^=1;reject_buffer(bad,bad.size());
 bad=bytes;bad.back()^=1;reject_buffer(bad,bad.size());
 // Repair the checksum to exercise semantic guards rather than integrity alone.
 auto reject_word=[&](size_t offset,uint64_t value){auto x=bytes;put_word(x,offset,value);fix_checksum(x);reject_buffer(x,x.size());};
 auto reject_real=[&](size_t offset,double value){auto x=bytes;put_real(x,offset,value);fix_checksum(x);reject_buffer(x,x.size());};
 reject_word(0,0);reject_word(8,999);reject_word(16,bytes.size()+8);
 bad=bytes;bad[24]^=1;fix_checksum(bad);reject_buffer(bad,bad.size());
 constexpr size_t header=88,edge_offset=header+38*8;
 reject_real(header,std::numeric_limits<double>::quiet_NaN());
 reject_real(header+8,lo); // invalid interval
 reject_real(header+2*8,1.); // claimed achieved error exceeds frozen tolerance
 reject_word(header+8*8,99); // unknown reaction
 reject_word(header+9*8,100001); // over-budget cells; reject before allocation
 reject_word(header+10*8,1); // insufficient knots
 reject_word(header+11*8,1000001); // direct evaluation budget
 reject_real(header+13*8,81.); // source CM cutoff outside direct API domain
 reject_word(header+35*8,0); // constructor knot budget
 reject_real(edge_offset+8,edges[0]); // repeated edge
 const size_t first_node=edge_offset+8*(kCells+1);
 reject_real(first_node,lo*.99); // missing exact lower endpoint
 reject_real(first_node+8,-1.); // negative reactivity
 reject_real(first_node+32*8,-1.); // negative first spectral cell
 bad=bytes;bad.insert(bad.end()-8,8,0);put_word(bad,16,bad.size());fix_checksum(bad);
 reject_buffer(bad,bad.size()); // extra unconsumed payload bytes

 // Serialization-only conserved perturbation: retain a subnormal alpha cell
 // too small to alter the existing node's macroscopic conservation check.
 std::vector<double> clean_grid(7*kCells);fusion_birth_coefficients_v1 clean_coeff{};
 expect_ok(fusion_c_birth_table_evaluate(original.get(),lo,kCells,clean_grid.data(),&clean_coeff),"clean knot");
 check(clean_grid[4*kCells]==0.,"synthetic tiny perturbation starts from zero cell");
 auto tiny=bytes;
 const size_t alpha_first_cell=first_node+32*8+4*kCells*8;
 put_real(tiny,alpha_first_cell,std::numeric_limits<double>::denorm_min());fix_checksum(tiny);
 p=nullptr;expect_ok(fusion_c_birth_table_unpack(tiny.data(),tiny.size(),&p),"subnormal node import");
 BytesTableOwner tiny_owner(p,fusion_c_birth_table_destroy);
 check(pack_table(p)==tiny,"subnormal bits preserved in repack");
 std::vector<double> tiny_grid(7*kCells);fusion_birth_coefficients_v1 tiny_coeff{};
 expect_ok(fusion_c_birth_table_evaluate(p,lo,kCells,tiny_grid.data(),&tiny_coeff),"subnormal knot eval");
 check(tiny_grid[4*kCells]==std::numeric_limits<double>::denorm_min(),"subnormal cell not erased");
 // Deliberately coarse, narrow-domain pB fixture: serialization coverage only,
 // not a claim that unit discrepancy gates validate a physical operating case.
 auto ps=source;ps.relative_max_J=2.5*kMeVJ;ps.relative_order=4;ps.cm_order=4;
 auto pc=control;pc.max_rate_error=pc.max_debit_error=pc.max_number_L1=pc.max_energy_L1=1.;
 pc.max_direct_rate_discrepancy=pc.max_direct_debit_discrepancy=1.;
 double pe[3]={0.,12.5*kMeVJ,25*kMeVJ};p=nullptr;
 expect_ok(fusion_c_birth_table_create(0,100*kKeVJ,100.001*kKeVJ,&ps,&pc,2,pe,&p),"real pB byte fixture");
 BytesTableOwner pb(p,fusion_c_birth_table_destroy);auto pbytes=pack_table(p);
 p=nullptr;expect_ok(fusion_c_birth_table_unpack(pbytes.data(),pbytes.size(),&p),"pB unpack");
 BytesTableOwner pbcopy(p,fusion_c_birth_table_destroy);check(pack_table(p)==pbytes,"pB exact repack");
 put_real(pbytes,header+12*8,10*kMeVJ);fix_checksum(pbytes);
 reject_buffer(pbytes,pbytes.size()); // static scalar domains pass; parent energy exceeds12MeV
 std::cout<<"Thermal byte roundtrip and basic rejection PASS\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
