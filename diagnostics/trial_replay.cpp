#include "detail/trial_snapshot.hpp"
#include "detail/wire.hpp"
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <new>
namespace pb11_diagnostics {
namespace {
using detail::need;
struct OwnedTable {
 int kind;void *handle;
 OwnedTable(int k,void *h):kind(k),handle(h){}
 ~OwnedTable(){if(!handle)return;if(kind==1)fusion_c_birth_table_destroy(static_cast<fusion_birth_table_v1 *>(handle));
  else fusion_c_beam_birth_table_destroy(static_cast<fusion_beam_birth_table_v1 *>(handle));}
};
class Tables {
 std::map<std::string,std::unique_ptr<OwnedTable>> tables_;
public:
 const void *load(const TableReference &ref) {
  if(!ref.present)return nullptr;
  need((ref.kind==1||ref.kind==2)&&ref.bytes>0&&
   ref.bytes<=(ref.kind==1?512ULL:256ULL)*1024*1024&&
   detail::valid_hex(ref.digest)&&ref.kernel==fusion_c_beam_birth_table_kernel_identity(),
   FUSION_CAPTURE_IDENTITY,"invalid replay table identity");
  auto key=std::to_string(ref.kind)+":"+std::to_string(ref.bytes)+":"+ref.digest+":"+ref.path;
  auto found=tables_.find(key);if(found!=tables_.end())return found->second->handle;
  std::error_code file_error;
  auto status=std::filesystem::symlink_status(ref.path,file_error);
  need(!file_error&&std::filesystem::exists(status),FUSION_CAPTURE_IO,"cannot locate replay table");
  need(std::filesystem::is_regular_file(status),FUSION_CAPTURE_IDENTITY,"replay table file kind changed");
  auto length=std::filesystem::file_size(ref.path,file_error);
  need(!file_error,FUSION_CAPTURE_IO,"cannot inspect replay table size");
  need(length==ref.bytes,FUSION_CAPTURE_IDENTITY,"replay table file size changed");
  std::ifstream file(ref.path,std::ios::binary);need(bool(file),FUSION_CAPTURE_IO,"cannot open replay table");
  std::vector<unsigned char> bytes(static_cast<size_t>(ref.bytes));
  need(bool(file.read(reinterpret_cast<char *>(bytes.data()),static_cast<std::streamsize>(bytes.size()))),
   FUSION_CAPTURE_IDENTITY,"truncated replay table");
  char extra;need(!file.get(extra)&&file.eof(),FUSION_CAPTURE_IDENTITY,"table file changed/trailing bytes");
  detail::Sha256 hash;hash.update(bytes.data(),bytes.size());
  need(detail::hex(hash.digest())==ref.digest,FUSION_CAPTURE_IDENTITY,"packed table digest differs");
  void *handle=nullptr;
  int code;
  if(ref.kind==1){fusion_birth_table_v1 *table=nullptr;
   code=fusion_c_birth_table_unpack(bytes.data(),bytes.size(),&table);handle=table;}
  else{fusion_beam_birth_table_v1 *table=nullptr;
   code=fusion_c_beam_birth_table_unpack(bytes.data(),bytes.size(),&table);handle=table;}
  // Hold a successfully imported handle before any allocation in the map.
  OwnedTable guard(ref.kind,handle);
  need(code==0&&handle,FUSION_CAPTURE_UNPACK_REJECTED,"original table unpack rejected reference");
  auto owner=std::make_unique<OwnedTable>(ref.kind,handle);
  guard.handle=nullptr;
  auto entry=tables_.emplace(std::move(key),std::move(owner));return entry.first->second->handle;
 }
};
}
Replay run_snapshot(const Snapshot &s) {
 need((s.entry==1||s.entry==2)&&s.cells>=1&&s.cells<=100000&&s.inert_count>=0&&s.inert_count<=32&&
  s.beam.size()<=size_t(10*s.cells)&&(!s.beam.size()||s.beam_array)&&(!s.inert_count||s.inert_array),
  FUSION_CAPTURE_FORMAT,"invalid replay layout");
 size_t count=size_t(s.cells)*6;
 need(s.kernel==fusion_c_beam_birth_table_kernel_identity()&&s.edges.size()==size_t(s.cells)+1&&
  s.n.size()==6&&s.z2.size()==6&&s.logs.size()==size_t(6*(7+s.inert_count))&&
  s.inert.size()==size_t(s.inert_count)&&s.s.size()==count&&s.t.size()==count&&
  s.birth.size()==count&&s.escape.size()==count,FUSION_CAPTURE_FORMAT,"invalid replay array shape/kernel");
 Tables loaded;
 std::array<const fusion_birth_table_v1 *,5> thermal{};
 for(size_t i=0;i<5;++i){need(s.thermal_array||!s.thermal[i].present,FUSION_CAPTURE_FORMAT,"table present without thermal handle array");
  need(!s.thermal[i].present||s.thermal[i].kind==1,FUSION_CAPTURE_FORMAT,"wrong thermal table kind");
  thermal[i]=static_cast<const fusion_birth_table_v1 *>(loaded.load(s.thermal[i]));}
 std::vector<fusion_beam_table_entry_v1> beam;
 beam.reserve(s.beam.size());
 for(const auto &reference:s.beam){need(!reference.table.present||reference.table.kind==2,FUSION_CAPTURE_FORMAT,"wrong beam table kind");
  beam.push_back({reference.channel,reference.slot,reference.cell,
   static_cast<const fusion_beam_birth_table_v1 *>(loaded.load(reference.table))});}
 fusion_beam_table_entry_v1 zero_beam{};fusion_inert_ion_v1 zero_inert{};
 const auto *beam_pointer=s.beam_array?(beam.empty()?&zero_beam:beam.data()):nullptr;
 const auto *inert_pointer=s.inert_array?(s.inert.empty()?&zero_inert:s.inert.data()):nullptr;
 const auto *thermal_pointer=s.thermal_array?thermal.data():nullptr;
 const auto *limits=s.floor_present?&s.floor:nullptr;
 Replay out;out.s.resize(count);out.t.resize(count);
 if(s.entry==2)out.mapped.resize(size_t(s.cells)*70);
 auto *floor=s.floor_present?&out.floor:nullptr;
 if(s.entry==1)out.status=fusion_c_coupled_sources_covered_trial(s.dt,&s.options,&s.fast,
  thermal_pointer,int(beam.size()),beam_pointer,s.effective,s.cells,s.edges.data(),s.n.data(),s.ue,s.ui,s.ne,
  s.z2.data(),s.inert_count,inert_pointer,s.logs.data(),s.s.data(),s.t.data(),s.birth.data(),s.escape.data(),
  out.n.data(),out.s.data(),out.t.data(),&out.result,&out.handoff,&out.usage,limits,floor,s.domain,&out.outside);
 else out.status=fusion_c_coupled_sources_packets_trial(s.dt,&s.options,&s.fast,
  thermal_pointer,int(beam.size()),beam_pointer,s.effective,s.cells,s.edges.data(),s.n.data(),s.ue,s.ui,s.ne,
  s.z2.data(),s.inert_count,inert_pointer,s.logs.data(),s.s.data(),s.t.data(),s.birth.data(),s.escape.data(),
  out.n.data(),out.s.data(),out.t.data(),&out.result,&out.handoff,&out.usage,limits,floor,s.domain,&out.outside,
  out.mapped.data(),&out.packets);
 return out;
}
} // namespace
extern "C" int fusion_capture_replay_file_v1(const char *path,uint64_t maximum,
 fusion_capture_replay_report_v1 *report) {
 if(report)std::memset(report,0,sizeof(*report));
 if(!report)return FUSION_CAPTURE_INVALID;
 try {
  auto s=pb11_diagnostics::read_snapshot(path,maximum);
  auto out=pb11_diagnostics::run_snapshot(s);
  auto digest=pb11_diagnostics::output_digest(out,s.entry,s.cells,s.floor_present);
  report->original_status=s.status;report->replay_status=out.status;
  report->status_matches=s.status==out.status;report->entry_point=s.entry;report->host_zone=s.zone;report->cells=s.cells;
  report->host_time_s=s.time;report->dt_s=s.dt;
  std::memcpy(report->kernel_identity,s.kernel.data(),64);
  std::memcpy(report->output_sha256,digest.data(),64);
  return FUSION_CAPTURE_OK;
 }catch(const pb11_diagnostics::detail::Error &error){return error.status;}
 catch(const std::bad_alloc &){return FUSION_CAPTURE_ALLOCATION;}
 catch(...){return FUSION_CAPTURE_IO;}
}
