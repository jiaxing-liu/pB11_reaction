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
void roundtrip(const fusion_capture_trial_input_v1&v,const Snapshot&s){
 require(s.entry==v.entry_point&&s.status==v.original_status&&s.zone==v.host_zone&&s.cells==v.cells&&s.effective==v.effective_charge&&s.domain==v.table_domain_policy&&s.inert_count==v.inert_count,"integer metadata");
 require(s.thermal_array==bool(v.thermal_tables)&&s.beam_array==bool(v.beam_tables)&&s.inert_array==bool(v.inert)&&s.floor_present==bool(v.floor_limits),"pointer presence");
 bit(s.time,v.host_time_s,"host time bits");bit(s.dt,v.dt_s,"dt bits");bit(s.ue,v.electron_energy_J_m3,"ue");bit(s.ui,v.ion_energy_J_m3,"ui");bit(s.ne,v.electron_density_m3,"ne");
 require(s.kernel==fusion_c_beam_birth_table_kernel_identity(),"kernel identity");options_equal(s.options,*v.options);
 for(int i=0;i<5;++i)require(s.fast.channels[i]==v.fast_options->channels[i],"fast mask");
 require(s.fast.angular_order==v.fast_options->angular_order,"angular order");bit(s.fast.angular_max_exponent,v.fast_options->angular_max_exponent,"angular exponent");
 if(v.floor_limits){bit(s.floor.max_center_over_ion_kT,v.floor_limits->max_center_over_ion_kT,"floor center");bit(s.floor.max_ion_energy_fraction,v.floor_limits->max_ion_energy_fraction,"floor fraction");}
 auto values=[](const std::vector<double>&a,const double*b,size_t n){require(a.size()==n,"input dimension");for(size_t i=0;i<n;++i)bit(a[i],b[i],"input element bits");};
 values(s.edges,v.edges_J,v.cells+1);values(s.n,v.thermal_number_m3,6);values(s.z2,v.thermal_charge_squared,6);values(s.logs,v.coulomb_logs,6*(7+v.inert_count));
 values(s.s,v.old_s_m3,6*v.cells);values(s.t,v.old_t_m3,6*v.cells);values(s.birth,v.external_birth_m3_s,6*v.cells);values(s.escape,v.escape_s_inv,6*v.cells);
 require(s.inert.size()==size_t(v.inert_count),"inert dimension");for(int i=0;i<v.inert_count;++i){bit(s.inert[i].density_m3,v.inert[i].density_m3,"inert density");bit(s.inert[i].mass_kg,v.inert[i].mass_kg,"inert mass");bit(s.inert[i].mean_charge_squared,v.inert[i].mean_charge_squared,"inert charge");}
 require(s.beam.size()==size_t(v.beam_table_count),"beam dimension");for(auto&t:s.thermal)require(!t.present,"direct thermal reference");
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
void checksum(std::vector<unsigned char>&b){auto h=fusion_detail::sha256(b.data(),b.size()-32);std::copy(h.begin(),h.end(),b.end()-32);}
void put(std::vector<unsigned char>&b,size_t at,uint64_t n,int width){for(int i=0;i<width;++i)b.at(at+i)=static_cast<unsigned char>(n>>(8*i));}
void bad(const fs::path&p,const std::vector<unsigned char>&b,int expected){save(p,b);try{(void)pb11_diagnostics::read_snapshot(p.c_str(),budget);require(false,"mutant accepted");}catch(const pb11_diagnostics::detail::Error&e){require(e.status==expected,"private reader error code");}fusion_capture_replay_report_v1 r{};require(fusion_capture_replay_file_v1(p.c_str(),budget,&r)==expected,"public reader error code");}
}

int main(int argc,char**argv){try{
 require(argc==2,"usage: test_trial_snapshot fresh-output-directory");fs::path dir=argv[1];require(fs::create_directory(dir),"output directory must be fresh/exclusive");
 fs::path good;
 for(int entry=1;entry<=2;++entry)for(int floor=0;floor<=1;++floor)for(int inert=0;inert<=1;++inert)for(int mode=0;mode<3;++mode){
  Inputs in(true);in.fast_options=make_fast_options({{0,0,0,1,0}});in.inert_count=inert;
  fusion_inert_ion_v1 bath{1e18,6.6446573357e-27,4};in.coulomb_logs.assign(6*(7+inert),15);
  if(mode==1){in.options.max_source_rate_error=0;in.options.max_source_debit_error=0;}
  if(mode==2)for(double&e:in.grid.edges)if(e!=0)e*=.001;
  const Inputs before=in;fusion_coupled_floor_limits_v1 limits{.5,.5};
  std::array<const fusion_birth_table_v1*,5> thermal{};fusion_beam_table_entry_v1 zero_beam{};
  fusion_capture_trial_input_v1 v{};v.entry_point=entry;v.host_zone=48;v.host_time_s=-0.0;v.dt_s=1e-4;v.options=&in.options;v.fast_options=&in.fast_options;
  v.thermal_tables=nullptr;v.beam_tables=&zero_beam;v.effective_charge=inert;v.cells=in.grid.cells();v.edges_J=in.grid.edges.data();v.thermal_number_m3=in.thermal_number.data();v.electron_energy_J_m3=in.electron_energy_J_m3;v.ion_energy_J_m3=in.ion_energy_J_m3;v.electron_density_m3=in.electron_density_m3;v.thermal_charge_squared=in.thermal_charge_squared.data();v.inert_count=inert;v.inert=&bath;v.coulomb_logs=in.coulomb_logs.data();v.old_s_m3=in.old_s.data();v.old_t_m3=in.old_t.data();v.external_birth_m3_s=in.external_birth.data();v.escape_s_inv=in.escape.data();v.floor_limits=floor?&limits:nullptr;v.table_domain_policy=FUSION_BEAM_TABLE_STRICT;
  auto baseline=direct(v);require(baseline.status==(mode==0?0:mode==1?4:3),"expected direct API numerical status: entry="+std::to_string(entry)+" floor="+std::to_string(floor)+" inert="+std::to_string(inert)+" mode="+std::to_string(mode)+" actual="+std::to_string(baseline.status));v.original_status=baseline.status;
  auto path=dir/("trial-"+std::to_string(entry)+std::to_string(floor)+std::to_string(inert)+std::to_string(mode)+".bin");
  require(fusion_capture_write_trial_v1(nullptr,path.c_str(),&v,budget)==0,"capture");auto snapshot=pb11_diagnostics::read_snapshot(path.c_str(),budget);roundtrip(v,snapshot);auto replay=pb11_diagnostics::run_snapshot(snapshot);outputs(baseline,replay);
  auto digest=pb11_diagnostics::output_digest(baseline,entry,v.cells,floor);require(digest==pb11_diagnostics::output_digest(replay,entry,v.cells,floor),"complete output bit digest");
  fusion_capture_replay_report_v1 report{};require(fusion_capture_replay_file_v1(path.c_str(),budget,&report)==0,"public replay");require(report.original_status==baseline.status&&report.replay_status==baseline.status&&report.status_matches==1&&report.entry_point==entry&&report.host_zone==v.host_zone&&report.cells==v.cells&&report.output_sha256==digest&&report.kernel_identity==snapshot.kernel,"public report metadata/digest");bit(report.host_time_s,v.host_time_s,"report signed zero time");bit(report.dt_s,v.dt_s,"report dt");
  auto original=bytes(path);require(fusion_capture_write_trial_v1(nullptr,path.c_str(),&v,budget)==FUSION_CAPTURE_IO&&bytes(path)==original,"exclusive no overwrite");
  auto tiny=dir/"tiny.bin";require(fusion_capture_write_trial_v1(nullptr,tiny.c_str(),&v,64)==FUSION_CAPTURE_CAPACITY&&!fs::exists(tiny),"tiny budget no final file");
  array_bits(in.old_s,before.old_s,"unchanged old S");array_bits(in.old_t,before.old_t,"unchanged old T");array_bits(in.grid.edges,before.grid.edges,"unchanged edges");array_bits(in.thermal_number,before.thermal_number,"unchanged thermal");array_bits(in.thermal_charge_squared,before.thermal_charge_squared,"unchanged charges");array_bits(in.coulomb_logs,before.coulomb_logs,"unchanged logs");array_bits(in.external_birth,before.external_birth,"unchanged birth");array_bits(in.escape,before.escape,"unchanged escape");options_equal(in.options,before.options);bit(in.electron_energy_J_m3,before.electron_energy_J_m3,"unchanged Ue");bit(in.ion_energy_J_m3,before.ion_energy_J_m3,"unchanged Ui");bit(in.electron_density_m3,before.electron_density_m3,"unchanged ne");require(in.fast_options.angular_order==before.fast_options.angular_order,"unchanged angular order");bit(in.fast_options.angular_max_exponent,before.fast_options.angular_max_exponent,"unchanged fast exponent");for(int i=0;i<5;++i)require(in.fast_options.channels[i]==before.fast_options.channels[i],"unchanged fast masks");require(bytes(path)==original,"replay input file immutable");if(good.empty())good=path;
  if(entry==1&&floor==0&&inert==0&&mode==0){
   // A nonnull all-null thermal table array requests table mode; it is
   // not equivalent to direct-mode nullptr when thermal DT is active.
   v.thermal_tables=thermal.data();
   auto invalid=direct(v);require(invalid.status==PB11_STATUS_INVALID_ARGUMENT,"all-null table mode rejects active thermal channel");
   v.original_status=invalid.status;
   auto negative=dir/"all-null-thermal-table-mode.bin";
   require(fusion_capture_write_trial_v1(nullptr,negative.c_str(),&v,budget)==0,"capture invalid table-mode trial");
   auto invalid_snapshot=pb11_diagnostics::read_snapshot(negative.c_str(),budget);roundtrip(v,invalid_snapshot);
   auto invalid_replay=pb11_diagnostics::run_snapshot(invalid_snapshot);outputs(invalid,invalid_replay);
   require(pb11_diagnostics::output_digest(invalid,entry,v.cells,floor)==pb11_diagnostics::output_digest(invalid_replay,entry,v.cells,floor),"invalid table-mode output bits");
   v.thermal_tables=nullptr;v.original_status=baseline.status;
   in.options.channels[3]=0;in.fast_options.channels[3]=0;
   in.options.birth.relative_phase=-0.0;in.options.birth.narrow_peak_fraction=std::numeric_limits<double>::denorm_min();
   in.old_s[0]=-0.0;in.old_t[0]=std::numeric_limits<double>::denorm_min();in.external_birth[0]=-0.0;in.escape[0]=std::numeric_limits<double>::denorm_min();
   v.thermal_tables=nullptr;v.beam_tables=nullptr;v.inert=nullptr;
   auto probe=dir/"binary64-codec-only.bin";require(fusion_capture_write_trial_v1(nullptr,probe.c_str(),&v,budget)==0,"codec-only capture");roundtrip(v,pb11_diagnostics::read_snapshot(probe.c_str(),budget));
  }
 }
 // Codec-only signed-zero/subnormal probe: disabled controls/populations are
 // roundtripped without a physics call or acceptance claim.
 auto original=bytes(good);
 auto mutate=[&](const char*name,auto edit,int code=FUSION_CAPTURE_FORMAT){auto b=original;edit(b);checksum(b);bad(dir/name,b,code);};
 auto corrupt=original;corrupt[250]^=1;bad(dir/"checksum.bin",corrupt,105);
 auto trunc=original;trunc.pop_back();bad(dir/"truncated.bin",trunc,105);
 auto trailing=original;trailing.push_back(0);bad(dir/"trailing.bin",trailing,105);
 mutate("version.bin",[](auto&b){put(b,8,2,4);});
 mutate("duplicate.bin",[](auto&b){put(b,88,1,4);});
 mutate("length.bin",[](auto&b){put(b,16,65,8);});
 mutate("cells.bin",[](auto&b){put(b,100+12,100001,4);});
 mutate("count.bin",[](auto&b){put(b,100+28,0x7fffffff,4);});
 mutate("kernel.bin",[](auto&b){b[24]=b[24]=='a'?'b':'a';},107);
 for(auto&item:fs::directory_iterator(dir))require(item.path().filename().string().find(".capture.")==std::string::npos,"temporary cleanup");
 std::cout<<"PASS: trial codec and original covered/packet status 0/3/4 parity; no campaign acceptance\n";return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
