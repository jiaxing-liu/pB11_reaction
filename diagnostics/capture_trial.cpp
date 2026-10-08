#include "detail/trial_snapshot.hpp"
#include "detail/wire.hpp"
#include <limits>
#include <new>
namespace pb11_diagnostics {
namespace {
using detail::need;
constexpr unsigned char magic[8]={'P','B','T','R','I','A','L','1'};
template<class W> void write_options(W &w,const fusion_coupled_thermal_options_v1 &o) {
 const auto &b=o.birth;
#define D(x) w.f64(b.x)
 D(relative_max_J);D(cm_max_kT);D(ground_state_q_J);D(cutoff_J);
 D(l1_fraction);D(relative_phase);D(narrow_peak_fraction);D(continuum_peak_scale);
#undef D
#define I(x) w.i32(b.x)
 I(continuation);I(pb_low);I(remainder_policy);I(broad_mode);I(fsci_policy);
 I(relative_order);I(cm_order);I(nq);I(ncos);
#undef I
 w.f64(o.max_source_rate_error);w.f64(o.max_source_debit_error);
 w.f64(o.handoff_max_L1);w.f64(o.handoff_max_mean_error);
 for(int value:o.channels)w.i32(value);
 w.i32(o.handoff_enabled);
}
void read_options(detail::Reader &r,fusion_coupled_thermal_options_v1 &o) {
 auto &b=o.birth;
#define D(x) b.x=r.f64()
 D(relative_max_J);D(cm_max_kT);D(ground_state_q_J);D(cutoff_J);
 D(l1_fraction);D(relative_phase);D(narrow_peak_fraction);D(continuum_peak_scale);
#undef D
#define I(x) b.x=r.i32()
 I(continuation);I(pb_low);I(remainder_policy);I(broad_mode);I(fsci_policy);
 I(relative_order);I(cm_order);I(nq);I(ncos);
#undef I
 o.max_source_rate_error=r.f64();o.max_source_debit_error=r.f64();
 o.handoff_max_L1=r.f64();o.handoff_max_mean_error=r.f64();
 for(auto &value:o.channels)value=r.i32();
 o.handoff_enabled=r.i32();
}
TableReference reference(fusion_capture_context_v1 *context,int kind,const void *handle) {
 TableReference out;
 if(!handle)return out;
 fusion_capture_table_identity_v1 identity{};
 int status=fusion_capture_table_identity_v1_get(context,kind,handle,&identity);
 need(status==0,status,"missing live table registration");
 out.present=true;out.kind=kind;out.bytes=identity.packed_bytes;
 out.digest=identity.content_sha256;out.kernel=identity.kernel_identity;out.path=identity.qualified_path;
 need(out.kernel==fusion_c_beam_birth_table_kernel_identity(),FUSION_CAPTURE_IDENTITY,"registered table kernel differs");
 return out;
}
uint64_t ref_size(const TableReference &ref) {return ref.present?148+ref.path.size():4;}
template<class W> void write_reference(W &w,const TableReference &ref) {
 w.u32(ref.present?1:0);if(!ref.present)return;
 w.i32(ref.kind);w.u64(ref.bytes);
 w.raw(ref.digest.data(),64);w.raw(ref.kernel.data(),64);
 w.u32(static_cast<uint32_t>(ref.path.size()));w.raw(ref.path.data(),ref.path.size());
}
TableReference read_reference(detail::Reader &r,int expected,const std::string &kernel) {
 TableReference ref;ref.present=r.flag();if(!ref.present)return ref;
 ref.kind=r.i32();ref.bytes=r.u64();ref.digest=r.text(64);ref.kernel=r.text(64);
 auto length=r.u32();need(length>1&&length<=4096,FUSION_CAPTURE_FORMAT,"invalid table path length");
 ref.path=r.text(length);
 need(ref.kind==expected&&ref.bytes>0&&ref.bytes<=(expected==1?512ULL:256ULL)*1024*1024,
  FUSION_CAPTURE_FORMAT,"invalid table kind/packed extent");
 need(detail::valid_hex(ref.digest)&&ref.kernel==kernel,FUSION_CAPTURE_IDENTITY,"table identity differs");
 need(ref.path[0]=='/'&&ref.path.find('\0')==std::string::npos,FUSION_CAPTURE_FORMAT,"absolute qualified table path required");
 return ref;
}
void write_trial(fusion_capture_context_v1 *context,const char *path,
 const fusion_capture_trial_input_v1 &v,uint64_t budget) {
 need((v.entry_point==1||v.entry_point==2)&&v.cells>=1&&v.cells<=100000&&
  v.inert_count>=0&&v.inert_count<=32&&v.beam_table_count>=0&&v.beam_table_count<=10*v.cells,
  FUSION_CAPTURE_INVALID,"unsupported entry/dimensions");
 need(v.options&&v.fast_options&&v.edges_J&&v.thermal_number_m3&&v.thermal_charge_squared&&
  v.coulomb_logs&&v.old_s_m3&&v.old_t_m3&&v.external_birth_m3_s&&v.escape_s_inv&&
  (!v.inert_count||v.inert)&&(!v.beam_table_count||v.beam_tables),FUSION_CAPTURE_INVALID,"required input pointer absent");
 std::vector<std::pair<uintptr_t,uintptr_t>> spans;
 auto span=[&](const void *pointer,size_t length){if(!length)return;
  auto start=reinterpret_cast<uintptr_t>(pointer);
  need(pointer&&length<=std::numeric_limits<uintptr_t>::max()-start,FUSION_CAPTURE_INVALID,"invalid input span");
  auto end=start+length;for(auto existing:spans)
   need(end<=existing.first||start>=existing.second,FUSION_CAPTURE_UNSUPPORTED,"aliased input layout cannot be recorded as independent arrays");
  spans.emplace_back(start,end);};
 span(v.options,sizeof(*v.options));span(v.fast_options,sizeof(*v.fast_options));
 if(v.thermal_tables)span(v.thermal_tables,5*sizeof(*v.thermal_tables));
 span(v.beam_tables,size_t(v.beam_table_count)*sizeof(*v.beam_tables));
 span(v.edges_J,(size_t(v.cells)+1)*8);span(v.thermal_number_m3,48);span(v.thermal_charge_squared,48);
 span(v.inert,size_t(v.inert_count)*sizeof(*v.inert));span(v.coulomb_logs,size_t(6*(7+v.inert_count))*8);
 for(auto pointer:{v.old_s_m3,v.old_t_m3,v.external_birth_m3_s,v.escape_s_inv})span(pointer,size_t(v.cells)*48);
 if(v.floor_limits)span(v.floor_limits,sizeof(*v.floor_limits));
 std::string kernel=fusion_c_beam_birth_table_kernel_identity();
 need(detail::valid_hex(kernel),FUSION_CAPTURE_IDENTITY,"invalid current kernel identity");
 detail::Writer w(path,budget);w.raw(magic,8);w.u32(1);
 w.field(1,64);w.raw(kernel.data(),64);w.end();
 w.field(2,88);
 for(int i:{v.entry_point,v.original_status,v.host_zone,v.cells,v.effective_charge,v.table_domain_policy,
  v.inert_count,v.beam_table_count,int(v.thermal_tables!=nullptr),int(v.beam_tables!=nullptr),
  int(v.inert!=nullptr),int(v.floor_limits!=nullptr)})w.i32(i);
 for(double d:{v.host_time_s,v.dt_s,v.electron_energy_J_m3,v.ion_energy_J_m3,v.electron_density_m3})w.f64(d);
 w.end();
 w.field(3,156);write_options(w,*v.options);w.end();
 w.field(4,32);for(int i:v.fast_options->channels)w.i32(i);
 w.i32(v.fast_options->angular_order);w.f64(v.fast_options->angular_max_exponent);w.end();
 w.field(5,v.floor_limits?16:0);
 if(v.floor_limits){w.f64(v.floor_limits->max_center_over_ion_kT);w.f64(v.floor_limits->max_ion_energy_fraction);}w.end();
 auto array=[&](uint32_t tag,const double *data,size_t count){w.field(tag,8*count);w.doubles(data,count);w.end();};
 array(6,v.edges_J,size_t(v.cells)+1);array(7,v.thermal_number_m3,6);array(8,v.thermal_charge_squared,6);
 w.field(9,uint64_t(v.inert_count)*24);
 for(int i=0;i<v.inert_count;++i){w.f64(v.inert[i].density_m3);w.f64(v.inert[i].mass_kg);w.f64(v.inert[i].mean_charge_squared);}w.end();
 array(10,v.coulomb_logs,size_t(6*(7+v.inert_count)));
 size_t count=size_t(v.cells)*6;
 array(11,v.old_s_m3,count);array(12,v.old_t_m3,count);
 array(13,v.external_birth_m3_s,count);array(14,v.escape_s_inv,count);
 uint64_t length=0;
 for(int i=0;i<5;++i)length+=ref_size(reference(context,1,v.thermal_tables?v.thermal_tables[i]:nullptr));
 w.field(15,length);
 for(int i=0;i<5;++i)write_reference(w,reference(context,1,v.thermal_tables?v.thermal_tables[i]:nullptr));
 w.end();length=0;
 for(int i=0;i<v.beam_table_count;++i){auto n=12+ref_size(reference(context,2,v.beam_tables[i].table));
  need(n<=budget&&length<=budget-n,FUSION_CAPTURE_CAPACITY,"beam reference field exceeds budget");length+=n;}
 w.field(16,length);
 for(int i=0;i<v.beam_table_count;++i){const auto &entry=v.beam_tables[i];
  w.i32(entry.channel);w.i32(entry.projectile_slot);w.i32(entry.energy_cell);
  write_reference(w,reference(context,2,entry.table));}
 w.end();w.publish();
}
}
Snapshot read_snapshot(const char *path,uint64_t budget) {
 detail::Reader r(path,budget);unsigned char actual[8];r.raw(actual,8);
 need(std::memcmp(actual,magic,8)==0&&r.u32()==1,FUSION_CAPTURE_FORMAT,"unsupported record magic/version");
 Snapshot s;r.field(1,64);s.kernel=r.text(64);r.end();
 need(s.kernel==fusion_c_beam_birth_table_kernel_identity(),FUSION_CAPTURE_IDENTITY,"snapshot kernel differs");
 r.field(2,88);s.entry=r.i32();s.status=r.i32();s.zone=r.i32();s.cells=r.i32();s.effective=r.i32();s.domain=r.i32();
 s.inert_count=r.i32();int beam_count=r.i32();s.thermal_array=r.flag();s.beam_array=r.flag();s.inert_array=r.flag();s.floor_present=r.flag();
 s.time=r.f64();s.dt=r.f64();s.ue=r.f64();s.ui=r.f64();s.ne=r.f64();r.end();
 need((s.entry==1||s.entry==2)&&s.cells>=1&&s.cells<=100000&&s.inert_count>=0&&s.inert_count<=32&&
  beam_count>=0&&beam_count<=10*s.cells&&(!beam_count||s.beam_array)&&(!s.inert_count||s.inert_array),
  FUSION_CAPTURE_FORMAT,"invalid shape/pointer-presence metadata");
 r.field(3,156);read_options(r,s.options);r.end();
 r.field(4,32);for(auto &value:s.fast.channels)value=r.i32();s.fast.angular_order=r.i32();s.fast.angular_max_exponent=r.f64();r.end();
 r.field(5,s.floor_present?16:0);if(s.floor_present){s.floor.max_center_over_ion_kT=r.f64();s.floor.max_ion_energy_fraction=r.f64();}r.end();
 auto array=[&](uint32_t tag,size_t count){r.field(tag,8*count);auto out=r.doubles(count);r.end();return out;};
 s.edges=array(6,size_t(s.cells)+1);s.n=array(7,6);s.z2=array(8,6);
 r.field(9,uint64_t(s.inert_count)*24);s.inert.resize(s.inert_count);
 for(auto &inert:s.inert){inert.density_m3=r.f64();inert.mass_kg=r.f64();inert.mean_charge_squared=r.f64();}r.end();
 s.logs=array(10,size_t(6*(7+s.inert_count)));size_t count=size_t(s.cells)*6;
 s.s=array(11,count);s.t=array(12,count);s.birth=array(13,count);s.escape=array(14,count);
 r.field(15);for(auto &table:s.thermal){table=read_reference(r,1,s.kernel);need(s.thermal_array||!table.present,FUSION_CAPTURE_FORMAT,"table without handle-array presence");}r.end();
 r.field(16);need(uint64_t(beam_count)<=r.remaining()/16,FUSION_CAPTURE_FORMAT,"beam count exceeds field before allocation");s.beam.reserve(beam_count);
 for(int i=0;i<beam_count;++i){BeamReference beam;beam.channel=r.i32();beam.slot=r.i32();beam.cell=r.i32();
  beam.table=read_reference(r,2,s.kernel);s.beam.push_back(std::move(beam));}
 r.end();r.finish();return s;
}
} // namespace
extern "C" int fusion_capture_write_trial_v1(fusion_capture_context_v1 *context,
 const char *path,const fusion_capture_trial_input_v1 *input,uint64_t maximum) {
 if(!input)return FUSION_CAPTURE_INVALID;
 try{pb11_diagnostics::write_trial(context,path,*input,maximum);return FUSION_CAPTURE_OK;}
 catch(const pb11_diagnostics::detail::Error &error){return error.status;}
 catch(const std::bad_alloc &){return FUSION_CAPTURE_ALLOCATION;}
 catch(...){return FUSION_CAPTURE_IO;}
}
