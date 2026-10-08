#pragma once
#include "fusion_capture_trial.h"
#include <array>
#include <string>
#include <vector>
namespace pb11_diagnostics {
struct TableReference {
 bool present=false;
 int kind=0;
 uint64_t bytes=0;
 std::string digest,kernel,path;
};
struct BeamReference {
 int channel=0,slot=0,cell=0;
 TableReference table;
};
struct Snapshot {
 int entry=0,status=0,zone=0,cells=0,effective=0,domain=0,inert_count=0;
 double time=0,dt=0,ue=0,ui=0,ne=0;
 std::string kernel;
 bool thermal_array=false,beam_array=false,inert_array=false,floor_present=false;
 fusion_coupled_thermal_options_v1 options{};
 fusion_fast_target_options_v1 fast{};
 fusion_coupled_floor_limits_v1 floor{};
 std::array<TableReference,5> thermal;
 std::vector<BeamReference> beam;
 std::vector<fusion_inert_ion_v1> inert;
 std::vector<double> edges,n,z2,logs,s,t,birth,escape;
};
struct Replay {
 int status=0;
 std::array<double,6> n{};
 std::vector<double> s,t,mapped;
 fusion_coupled_thermal_v1 result{};
 fusion_handoff_diagnostics_v1 handoff{};
 fusion_beam_table_usage_v1 usage{};
 fusion_birth_floor_ledger_v1 floor{};
 fusion_birth_packets_v1 packets{};
 uint64_t outside=0;
};
// Internal owned decoder/replay surfaces support complete bitwise tests.
Snapshot read_snapshot(const char *path,uint64_t max_bytes);
Replay run_snapshot(const Snapshot &snapshot);
std::string output_digest(const Replay &out,int entry,int cells,bool floor);
}
