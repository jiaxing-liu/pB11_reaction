#define main budget_fixture_main
#include "test_fusion_coupled_sources_budget.cpp"
#undef main
#include "detail/trial_snapshot.hpp"
#include "detail/wire.hpp"
#include "sha256_internal.hpp"
#include <filesystem>
#include <fstream>
namespace fs=std::filesystem;
constexpr uint64_t cap=16*1024*1024;
fusion_capture_trial_input_v1 view(const Input&i,int entry,const fusion_coupled_floor_limits_v1&limits){
 fusion_capture_trial_input_v1 v{};v.entry_point=entry;v.dt_s=1e-8;v.options=&i.options;v.fast_options=&i.fast;
 v.cells=n;v.edges_J=i.edges.data();v.thermal_number_m3=i.thermal.data();v.electron_energy_J_m3=1.5*i.ne*i.temperature;
 v.ion_energy_J_m3=i.ui();v.electron_density_m3=i.ne;v.thermal_charge_squared=i.charge.data();v.coulomb_logs=i.logs.data();
 v.old_s_m3=i.s.data();v.old_t_m3=i.t.data();v.external_birth_m3_s=i.external.data();v.escape_s_inv=i.escape.data();
 v.floor_limits=&limits;v.table_domain_policy=FUSION_BEAM_TABLE_STRICT;return v;
}
int covered(const Input&i,Output&o,const fusion_thermal_birth_budget_policy_v1&p,const fusion_coupled_floor_limits_v1&l,const fusion_birth_table_v1*const*tables=nullptr){
 return fusion_c_coupled_sources_covered_budget_trial(1e-8,&i.options,&i.fast,tables,0,nullptr,0,n,i.edges.data(),i.thermal.data(),1.5*i.ne*i.temperature,i.ui(),i.ne,i.charge.data(),0,nullptr,i.logs.data(),i.s.data(),i.t.data(),i.external.data(),i.escape.data(),o.thermal.data(),o.s.data(),o.t.data(),&o.result,&o.diagnostics,&o.usage,&l,&o.floor,FUSION_BEAM_TABLE_STRICT,&o.outside,&p,&o.accuracy);
}
void check(const Output&o,const pb11_diagnostics::Replay&r){
 need(std::equal(o.thermal.begin(),o.thermal.end(),r.n.begin())&&std::equal(o.s.begin(),o.s.end(),r.s.begin())&&std::equal(o.t.begin(),o.t.end(),r.t.begin()),"independent physical arrays");
 need(same(o.result,r.result)&&same(o.diagnostics,r.handoff)&&same(o.usage,r.usage)&&same(o.floor,r.floor)&&same(o.accuracy,r.accuracy)&&o.outside==r.outside,"independent full observation fields");
 if(!r.mapped.empty())need(std::equal(o.mapped.begin(),o.mapped.end(),r.mapped.begin())&&same(o.packets,r.packets),"independent packet arrays");
}
std::vector<unsigned char> load(const fs::path&p){std::ifstream f(p,std::ios::binary);return {std::istreambuf_iterator<char>(f),{}};}
void malformed(const fs::path&p,std::vector<unsigned char>b,bool fix=true){
 if(fix){auto h=fusion_detail::sha256(b.data(),b.size()-32);std::copy(h.begin(),h.end(),b.end()-32);}
 std::ofstream f(p,std::ios::binary);f.write(reinterpret_cast<const char*>(b.data()),b.size());f.close();
 fusion_capture_replay_report_v1 r{};need(fusion_capture_replay_file_v1(p.c_str(),cap,&r)==FUSION_CAPTURE_FORMAT,"malformed rejects");
}
void table_snapshot(const fs::path&dir){
 Input i;fusion_coupled_floor_limits_v1 limits{.5,.5};
 fusion_thermal_birth_budget_policy_v1 policy{1,4,4,0,1e-4,1e-4};
 fusion_birth_table_control_v1 control{.01,.01,.01,.01,.05,.05,8,64,3};
 const double lo=i.temperature,hi=lo*1.000001;fusion_birth_table_v1*raw=nullptr;
 need(fusion_c_pb_birth_table_create_budget(lo,hi,&i.options.birth,&policy,&control,n,i.edges.data(),&raw)==0,"capture fixture budget table create");
 using TableOwner=std::unique_ptr<fusion_birth_table_v1,decltype(&fusion_c_birth_table_destroy)>;
 TableOwner original(raw,fusion_c_birth_table_destroy);size_t size=0,written=0;
 need(fusion_c_birth_table_pack_size(raw,&size)==0&&size>0&&size<cap,"budget table packed size");
 std::vector<unsigned char> packed(size);need(fusion_c_birth_table_pack(raw,packed.data(),packed.size(),&written)==0&&written==size,"budget table pack");
 const auto tablepath=fs::absolute(dir/"budget-table.bin");
 auto save_table=[&](const std::vector<unsigned char>&b){std::ofstream f(tablepath,std::ios::binary);f.write(reinterpret_cast<const char*>(b.data()),b.size());need(bool(f),"write budget table");};
 save_table(packed);original.reset();
 fusion_capture_context_v1*context=nullptr;need(fusion_capture_context_create_v1(1,&context)==0,"table fixture registry");
 using ContextOwner=std::unique_ptr<fusion_capture_context_v1,decltype(&fusion_capture_context_destroy_v1)>;
 ContextOwner owner(context,fusion_capture_context_destroy_v1);fusion_birth_table_v1*registered=nullptr;int unpack=-1;
 need(fusion_capture_thermal_unpack_v1(context,packed.data(),packed.size(),tablepath.c_str(),&registered,&unpack)==0&&unpack==0&&registered,"registered budget table unpack");
 fusion_birth_table_budget_info_v1 info{};need(fusion_c_pb_birth_table_budget_info(registered,&info)==0,"imported budget table metadata");
 std::array<const fusion_birth_table_v1*,5> tables{{registered,nullptr,nullptr,nullptr,nullptr}};
 i.temperature=std::sqrt(lo*hi);
 for(int entry:{3,4})for(int wrong=0;wrong<2;++wrong){
  auto p=policy;if(wrong)p.energy_L1_tolerance=std::nextafter(p.energy_L1_tolerance,1.);
  Output direct;std::memset(&direct,0x5a,sizeof direct);
  int status=entry==3?covered(i,direct,p,limits,tables.data()):call(i,direct,&p,tables.data());
  need(wrong?status!=0:status==0,"independent table trial policy status");
  if(!wrong)need(direct.result.ledger.events_m3[0]>0&&direct.accuracy.version==1&&direct.accuracy.mode==2&&zero(direct.accuracy.direct)&&same(direct.accuracy.table,info),"successful sampled mode2 actual pB and no fictitious direct report");
  else need(zero(direct.accuracy),"wrong table policy clears accuracy");
  auto v=view(i,entry,limits);v.thermal_tables=tables.data();v.original_status=status;
  const auto path=dir/("table-"+std::to_string(entry)+"-"+std::to_string(wrong)+".bin");
  need(fusion_capture_write_trial_v2(context,path.c_str(),&v,&p,cap)==0,"budget table snapshot capture");
  auto snapshot=pb11_diagnostics::read_snapshot(path.c_str(),cap);
  need(snapshot.thermal_array&&snapshot.thermal[0].present&&snapshot.thermal[0].bytes==size&&snapshot.thermal[0].path==tablepath.string(),"budget table reference metadata");
  auto replay=pb11_diagnostics::run_snapshot(snapshot);need(replay.status==status,"budget table replay status");check(direct,replay);
  fusion_capture_replay_report_v1 report{};need(fusion_capture_replay_file_v1(path.c_str(),cap,&report)==0&&report.status_matches,"budget table public replay");
  auto changed=packed;changed[changed.size()/2]^=1;save_table(changed);
  need(fusion_capture_replay_file_v1(path.c_str(),cap,&report)==FUSION_CAPTURE_IDENTITY,"changed packed budget bytes reject before dispatch");
  save_table(packed);need(fusion_capture_replay_file_v1(path.c_str(),cap,&report)==0&&report.status_matches,"restored budget bytes replay");
 }
 // Active pB uses direct mode with nullptr, whereas a present all-null table
 // array requests table mode and fails. Preserve that distinction in v2 too.
 for(int entry:{3,4}){
  tables[0]=nullptr;Output direct;
  int status=entry==3?covered(i,direct,policy,limits,tables.data()):call(i,direct,&policy,tables.data());
  need(status!=0&&zero(direct.accuracy),"all-null budget table array rejects active pB");
  auto v=view(i,entry,limits);v.thermal_tables=tables.data();v.original_status=status;
  const auto path=dir/("allnull-table-"+std::to_string(entry)+".bin");
  need(fusion_capture_write_trial_v2(nullptr,path.c_str(),&v,&policy,cap)==0,"all-null table mode capture");
  auto snapshot=pb11_diagnostics::read_snapshot(path.c_str(),cap);need(snapshot.thermal_array&&!snapshot.thermal[0].present,"all-null handle array presence retained");
  auto replay=pb11_diagnostics::run_snapshot(snapshot);need(replay.status==status,"all-null table mode replay status");check(direct,replay);
 }
 std::cout<<"PASS serialized budget table used by entries3/4 mode2, exact policy mismatch, altered-file identity guard, all-null table distinction\n";
}
int main(int argc,char**argv){try{
 need(argc==2&&fs::create_directory(argv[1]),"fresh directory required");fs::path dir=argv[1],good;
 for(int entry:{3,4})for(int mode=0;mode<6;++mode){
  Input i;fusion_coupled_floor_limits_v1 l{.5,.5};fusion_thermal_birth_budget_policy_v1 p{1,4,4,0,0,0};
  if(mode==1)p.number_L1_tolerance=p.energy_L1_tolerance=2;
  if(mode==2)p.number_L1_tolerance=p.energy_L1_tolerance=1e-4;
  if(mode==3){i.s[n/2]=1e10;i.t[n/2]=1e10;p.number_L1_tolerance=p.energy_L1_tolerance=2;}
  if(mode==4)p.version=99;
  if(mode==5){uint64_t nan=0x7ff8000000000042ULL;std::memcpy(&p.number_L1_tolerance,&nan,8);p.reserved=17;}
  Output baseline;std::memset(&baseline,0x5a,sizeof baseline);
  int status=entry==3?covered(i,baseline,p,l):call(i,baseline,&p);
  need((mode<4)==(status==0),"expected numerical status");
  if(status==0){need(baseline.result.ledger.events_m3[0]>0&&baseline.accuracy.mode==1,"actual nonzero pB source");
   if(mode==1)need(baseline.accuracy.direct.low_pairs>0,"actual low parents");
  }else need(zero(baseline.accuracy),"invalid policy accuracy clearing");
  auto v=view(i,entry,l);v.original_status=status;auto path=dir/(std::to_string(entry)+"-"+std::to_string(mode)+".bin");
  need(fusion_capture_write_trial_v2(nullptr,path.c_str(),&v,&p,cap)==0,"v2 capture");
  auto s=pb11_diagnostics::read_snapshot(path.c_str(),cap);need(s.schema==2&&s.policy_present&&same(s.policy,p),"exact invalid and valid policy bits");
  auto r=pb11_diagnostics::run_snapshot(s);need(r.status==status,"replay status");check(baseline,r);
  fusion_capture_replay_report_v1 report{};need(fusion_capture_replay_file_v1(path.c_str(),cap,&report)==0&&report.status_matches,"public v2 replay");
  auto digest=pb11_diagnostics::output_digest(r,entry,n,true);r.accuracy.table.base.control.max_depth++;
  need(digest!=pb11_diagnostics::output_digest(r,entry,n,true),"accuracy included in digest");
  need(fusion_capture_write_trial_v1(nullptr,(dir/"wrong-v1").c_str(),&v,cap)==FUSION_CAPTURE_INVALID,"v1 rejects new entry");
  auto oldentry=v;oldentry.entry_point=1;need(fusion_capture_write_trial_v2(nullptr,(dir/"wrong-v2").c_str(),&oldentry,&p,cap)==FUSION_CAPTURE_INVALID,"v2 rejects old entry");
  need(fusion_capture_write_trial_v2(nullptr,(dir/"null").c_str(),&v,nullptr,cap)==FUSION_CAPTURE_INVALID,"mandatory policy");
  need(fusion_capture_write_trial_v2(nullptr,(dir/"alias").c_str(),&v,reinterpret_cast<const fusion_thermal_birth_budget_policy_v1*>(v.options),cap)==FUSION_CAPTURE_UNSUPPORTED,"overlap policy rejects");
  need(fusion_capture_write_trial_v2(nullptr,path.c_str(),&v,&p,cap)==FUSION_CAPTURE_IO,"exclusive file");
  need(fusion_capture_write_trial_v2(nullptr,(dir/"tiny").c_str(),&v,&p,64)==FUSION_CAPTURE_CAPACITY,"bounded capture");
  if(good.empty())good=path;
 }
 table_snapshot(dir);
 auto original=load(good);
 auto b=original;b[8]=1;malformed(dir/"v1-newentry",b);
 b=original;b[8]=3;malformed(dir/"unknown-version",b);
 b=original;b[100]=1;malformed(dir/"v2-oldentry",b);
 b=original;b[b.size()-76]=18;malformed(dir/"policy-tag",b);
 b=original;b[b.size()-72]=31;malformed(dir/"policy-length",b);
 b=original;b[250]^=1;malformed(dir/"checksum",b,false);
 b=original;b.pop_back();malformed(dir/"truncated",b,false);
 b=original;b.erase(b.end()-76,b.end()-32);malformed(dir/"missing-policy",b);
 std::cout<<"PASS budget snapshot independent covered/packet pB success/failure and malformed records\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
