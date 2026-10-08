// Identity-only real packed tables: every reaction is disabled and fast pools
// are empty. Neither table is used by a nuclear reaction; this is not a
// used-cache physics test or campaign acceptance. No table constructor runs.
// Identity import succeeds, then the original API faithfully rejects disabled
// beam routes (and duplicate routes) as INVALID_ARGUMENT; cleared output parity
// verifies that rejection, not a successful physical trial.
// Reuse the existing numerical fixture helpers without calling its main.
#define main fixture_original_main
#include "test_fusion_coupled_fast.cpp"
#undef main
#include "detail/trial_snapshot.hpp"
#include "detail/wire.hpp"
#include "sha256_internal.hpp"
#include <filesystem>
#include <fstream>
#include <cstring>

namespace {
using pb11_diagnostics::Replay;
using pb11_diagnostics::Snapshot;
constexpr uint64_t budget=16*1024*1024;
namespace fs=std::filesystem;
void bit(double a,double b,const char *label){require(pb11_diagnostics::detail::bits(a)==pb11_diagnostics::detail::bits(b),label);}
template<class A,class B> void array_bits(const A &a,const B &b,const char *label){
 require(std::size(a)==std::size(b),label);for(size_t i=0;i<std::size(a);++i)bit(a[i],b[i],label);
}
void options_equal(const fusion_coupled_thermal_options_v1&a,const fusion_coupled_thermal_options_v1&b){
#define D(x) bit(a.birth.x,b.birth.x,#x)
 D(relative_max_J);D(cm_max_kT);D(ground_state_q_J);D(cutoff_J);D(l1_fraction);D(relative_phase);D(narrow_peak_fraction);D(continuum_peak_scale);
#undef D
#define I(x) require(a.birth.x==b.birth.x,#x)
 I(continuation);I(pb_low);I(remainder_policy);I(broad_mode);I(fsci_policy);I(relative_order);I(cm_order);I(nq);I(ncos);
#undef I
#define D(x) bit(a.x,b.x,#x)
 D(max_source_rate_error);D(max_source_debit_error);D(handoff_max_L1);D(handoff_max_mean_error);
#undef D
 for(int i=0;i<5;++i)require(a.channels[i]==b.channels[i],"thermal mask");
 require(a.handoff_enabled==b.handoff_enabled,"handoff flag");
}
Replay direct(const fusion_capture_trial_input_v1&v){
 Replay o;o.s.resize(6*v.cells);o.t.resize(6*v.cells);if(v.entry_point==2)o.mapped.resize(70*v.cells);
#define ARGS v.dt_s,v.options,v.fast_options,v.thermal_tables,v.beam_table_count,v.beam_tables,v.effective_charge,v.cells,v.edges_J,v.thermal_number_m3,v.electron_energy_J_m3,v.ion_energy_J_m3,v.electron_density_m3,v.thermal_charge_squared,v.inert_count,v.inert,v.coulomb_logs,v.old_s_m3,v.old_t_m3,v.external_birth_m3_s,v.escape_s_inv,o.n.data(),o.s.data(),o.t.data(),&o.result,&o.handoff,&o.usage,v.floor_limits,v.floor_limits?&o.floor:nullptr,v.table_domain_policy,&o.outside
 if(v.entry_point==1)o.status=fusion_c_coupled_sources_covered_trial(ARGS);
 else o.status=fusion_c_coupled_sources_packets_trial(ARGS,o.mapped.data(),&o.packets);
#undef ARGS
 return o;
}
void outputs(const Replay&a,const Replay&b){
 require(a.status==b.status&&a.outside==b.outside,"status/counter");array_bits(a.n,b.n,"thermal");array_bits(a.s,b.s,"S");array_bits(a.t,b.t,"T");array_bits(a.mapped,b.mapped,"mapped");
 // Fixture comparators enumerate all declared result, diagnostic and usage members.
 require(same_result(a.result,b.result)&&same_diagnostics(a.handoff,b.handoff)&&same_usage(a.usage,b.usage),"all ordinary output members");
#define A(x) array_bits(a.floor.x,b.floor.x,#x)
 A(born_number_m3);A(born_energy_J_m3);A(mapped_number_m3);A(mapped_energy_J_m3);A(ion_energy_correction_J_m3);A(energy_residual_J_m3);
#undef A
 bit(a.floor.remaining_ion_energy_J_m3,b.floor.remaining_ion_energy_J_m3,"remaining floor energy");
 for(int src=0;src<2;++src)for(int ch=0;ch<5;++ch){bit(a.packets.events_m3[src][ch],b.packets.events_m3[src][ch],"packet events");
#define A(x) array_bits(a.packets.x[src][ch],b.packets.x[src][ch],#x)
 A(below_number_m3);A(below_energy_J_m3);A(above_number_m3);A(above_energy_J_m3);
#undef A
 }
}
std::vector<unsigned char> bytes(const fs::path&p){std::ifstream f(p,std::ios::binary);return {std::istreambuf_iterator<char>(f),{}};}
void save(const fs::path&p,const std::vector<unsigned char>&b){std::ofstream f(p,std::ios::binary);f.write(reinterpret_cast<const char*>(b.data()),b.size());require(bool(f),"write mutant");}

}

#include <memory>

namespace {
using RegistryOwner=std::unique_ptr<fusion_capture_context_v1,decltype(&fusion_capture_context_destroy_v1)>;
RegistryOwner registry(uint64_t capacity){fusion_capture_context_v1*p=nullptr;
 require(fusion_capture_context_create_v1(capacity,&p)==0&&p,"registry create");
 return RegistryOwner(p,fusion_capture_context_destroy_v1);}
std::string table_sha(const std::vector<unsigned char>&b){return pb11_diagnostics::detail::hex(fusion_detail::sha256(b.data(),b.size()));}
struct PackedInput{fs::path path;std::vector<unsigned char> data;std::string sha;};
PackedInput private_copy(const fs::path&source,const fs::path&dest,uint64_t limit){
 require(fs::is_regular_file(source)&&fs::file_size(source)>0&&fs::file_size(source)<=limit,"source packed size");
 PackedInput p{fs::absolute(dest).lexically_normal(),bytes(source),{}};p.sha=table_sha(p.data);
 require(!fs::exists(p.path),"private packed path fresh");save(p.path,p.data);return p;
}
void reference_equal(const pb11_diagnostics::TableReference&r,int kind,const PackedInput&p){
 require(r.present&&r.kind==kind&&r.bytes==p.data.size()&&r.digest==p.sha&&r.path==p.path.string()&&r.kernel==fusion_c_beam_birth_table_kernel_identity(),"complete packed reference identity");
}
void registry_identity(fusion_capture_context_v1*c,int kind,const void*h,const PackedInput&p){
 fusion_capture_table_identity_v1 r{};require(fusion_capture_table_identity_v1_get(c,kind,h,&r)==0,"registry identity");
 require(r.kind==kind&&r.packed_bytes==p.data.size()&&std::string(r.content_sha256)==p.sha&&std::string(r.qualified_path)==p.path.string()&&std::string(r.kernel_identity)==fusion_c_beam_birth_table_kernel_identity(),"registry metadata/hash/path");
}
void failure_report(const fs::path&p,int expected){fusion_capture_replay_report_v1 r;
 std::memset(&r,0xa5,sizeof(r));require(fusion_capture_replay_file_v1(p.c_str(),budget,&r)==expected,"identity failure before physics");
 const auto*b=reinterpret_cast<const unsigned char*>(&r);for(size_t i=0;i<sizeof(r);++i)require(b[i]==0,"failure report cleared");
}
void unknown_capture(fusion_capture_context_v1*c,const fs::path&p,const fusion_capture_trial_input_v1&v,int expected=FUSION_CAPTURE_UNKNOWN_HANDLE){
 require(fusion_capture_write_trial_v1(c,p.c_str(),&v,budget)==expected&&!fs::exists(p),"unknown registration no publication");
}
}
int main(int argc,char**argv){try{
 require(argc==4,"usage: test_trial_tables fresh-output-directory thermal-packed-file beam-packed-file");
 const fs::path thermal_source=argv[2],beam_source=argv[3];
 fs::path dir=fs::absolute(argv[1]).lexically_normal();require(fs::create_directory(dir),"fresh exclusive directory");
 auto thermal=private_copy(thermal_source,dir/"thermal-private.bin",512ULL*1024*1024);
 auto beam=private_copy(beam_source,dir/"beam-private.bin",256ULL*1024*1024);
 auto ctx=registry(2);fusion_birth_table_v1*th=nullptr;fusion_beam_birth_table_v1*bh=nullptr;int original=-1;
 require(fusion_capture_thermal_unpack_v1(ctx.get(),thermal.data.data(),thermal.data.size(),thermal.path.c_str(),&th,&original)==0&&original==0&&th,"exact loader thermal unpack");
 require(fusion_capture_beam_unpack_v1(ctx.get(),beam.data.data(),beam.data.size(),beam.path.c_str(),&bh,&original)==0&&original==0&&bh,"exact loader beam unpack");
 registry_identity(ctx.get(),1,th,thermal);registry_identity(ctx.get(),2,bh,beam);
 Inputs in(true);in.options=make_options({{0,0,0,0,0}});in.fast_options=make_fast_options({{0,0,0,0,0}});
 // 508 cells makes the recorded D519 route (3,1,507) a valid cell index.
 in.grid.edges.resize(509);in.grid.edges[0]=0;
 for(int i=1;i<=508;++i)in.grid.edges[i]=2*kElectronVoltJ*std::pow(25*kMeVJ/(2*kElectronVoltJ),double(i-1)/507);
 in.old_s.assign(6*508,0);in.old_t.assign(6*508,0);in.external_birth.assign(6*508,0);in.escape.assign(6*508,0);
 const Inputs before=in;
 std::array<const fusion_birth_table_v1*,5> tables{{th,nullptr,nullptr,nullptr,nullptr}};
 // Duplicate reference retains one live registry handle; replay's Tables map
 // deduplicates by complete identity. This is not a table usage counter claim.
 std::array<fusion_beam_table_entry_v1,2> routes{{{3,1,507,bh},{3,1,507,bh}}};
 fusion_capture_trial_input_v1 v{};v.host_zone=73;v.host_time_s=-0.0;v.dt_s=1e-4;v.options=&in.options;v.fast_options=&in.fast_options;
 v.thermal_tables=tables.data();v.beam_table_count=2;v.beam_tables=routes.data();v.cells=in.grid.cells();v.edges_J=in.grid.edges.data();v.thermal_number_m3=in.thermal_number.data();v.electron_energy_J_m3=in.electron_energy_J_m3;v.ion_energy_J_m3=in.ion_energy_J_m3;v.electron_density_m3=in.electron_density_m3;v.thermal_charge_squared=in.thermal_charge_squared.data();v.coulomb_logs=in.coulomb_logs.data();v.old_s_m3=in.old_s.data();v.old_t_m3=in.old_t.data();v.external_birth_m3_s=in.external_birth.data();v.escape_s_inv=in.escape.data();v.table_domain_policy=FUSION_BEAM_TABLE_STRICT;
 require(routes[0].energy_cell<v.cells,"route cell valid");
 fs::path good;
 for(int entry=1;entry<=2;++entry){v.entry_point=entry;auto baseline=direct(v);require(baseline.status==PB11_STATUS_INVALID_ARGUMENT,"original API rejects disabled/duplicate beam routes");v.original_status=baseline.status;
 auto path=dir/("tables-"+std::to_string(entry)+".bin");require(fusion_capture_write_trial_v1(ctx.get(),path.c_str(),&v,budget)==0,"real reference capture");
 auto immutable=bytes(path);auto s=pb11_diagnostics::read_snapshot(path.c_str(),budget);
 require(s.thermal_array&&s.beam_array&&s.cells==508&&s.beam.size()==2&&s.entry==entry&&s.zone==v.host_zone,"snapshot metadata");
 reference_equal(s.thermal[0],1,thermal);for(int ch=1;ch<5;++ch)require(!s.thermal[ch].present,"five table array null slots");
 for(auto&r:s.beam){require(r.channel==3&&r.slot==1&&r.cell==507,"beam route bits");reference_equal(r.table,2,beam);}
 auto replay=pb11_diagnostics::run_snapshot(s);outputs(baseline,replay);fusion_beam_table_usage_v1 zero{};
 require(same_usage(baseline.usage,zero)&&same_usage(replay.usage,zero)&&baseline.outside==0&&replay.outside==0,"identity only zero usage");
 auto digest=pb11_diagnostics::output_digest(baseline,entry,v.cells,false);require(digest==pb11_diagnostics::output_digest(replay,entry,v.cells,false),"all output digest parity");
 fusion_capture_replay_report_v1 r{};require(fusion_capture_replay_file_v1(path.c_str(),budget,&r)==0&&r.original_status==baseline.status&&r.replay_status==baseline.status&&r.status_matches==1&&r.entry_point==entry&&r.host_zone==v.host_zone&&r.cells==v.cells&&std::string(r.output_sha256)==digest&&std::string(r.kernel_identity)==s.kernel,"public report parity");bit(r.host_time_s,v.host_time_s,"report time bits");bit(r.dt_s,v.dt_s,"report dt bits");
 require(fusion_capture_write_trial_v1(ctx.get(),path.c_str(),&v,budget)==106&&bytes(path)==immutable,"exclusive immutable evidence");
 for(const auto*p:{&thermal,&beam}){auto mutant=p->data;mutant[mutant.size()/2]^=1;save(p->path,mutant);require(fs::file_size(p->path)==p->data.size(),"same size mutation");failure_report(path,107);save(p->path,p->data);require(fusion_capture_replay_file_v1(path.c_str(),budget,&r)==0&&std::string(r.output_sha256)==digest,"restore identity restores parity");
 auto moved=p->path.string()+".away";fs::rename(p->path,moved);failure_report(path,106);fs::rename(moved,p->path);require(fusion_capture_replay_file_v1(path.c_str(),budget,&r)==0&&std::string(r.output_sha256)==digest,"restore moved file restores parity");}
 require(bytes(path)==immutable,"snapshot immutable across failures");if(good.empty())good=path;
 }
 auto other=registry(2);unknown_capture(other.get(),dir/"unknown-context.bin",v);unknown_capture(nullptr,dir/"null-context.bin",v,FUSION_CAPTURE_INVALID);
 const auto*stale=th;require(fusion_capture_thermal_destroy_v1(ctx.get(),&th)==0&&!th,"destroy registered thermal");tables[0]=stale; // Address token only: never dereference a destroyed handle.
 unknown_capture(ctx.get(),dir/"stale-handle.bin",v);
 array_bits(in.old_s,before.old_s,"input S immutable");array_bits(in.old_t,before.old_t,"input T immutable");array_bits(in.grid.edges,before.grid.edges,"input edges immutable");array_bits(in.thermal_number,before.thermal_number,"input N immutable");array_bits(in.thermal_charge_squared,before.thermal_charge_squared,"input charges immutable");array_bits(in.coulomb_logs,before.coulomb_logs,"input logs immutable");array_bits(in.external_birth,before.external_birth,"input births immutable");array_bits(in.escape,before.escape,"input escape immutable");options_equal(in.options,before.options);
 require(bytes(thermal_source)==thermal.data&&bytes(beam_source)==beam.data&&bytes(thermal.path)==thermal.data&&bytes(beam.path)==beam.data,"source and private table bytes preserved");
 for(auto&item:fs::directory_iterator(dir))require(item.path().filename().string().find(".capture.")==std::string::npos,"temporary publication cleanup");
 std::cout<<"PASS: real packed-table identity capture/replay; original invalid status preserved, usage zero; no physics acceptance\n";return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
