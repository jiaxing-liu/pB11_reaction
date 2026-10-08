#include "detail/trial_snapshot.hpp"
#include "detail/wire.hpp"

namespace pb11_diagnostics {
std::string output_digest(const Replay &out,int entry,int cells,bool floor) {
 detail::need(entry==1||entry==2,FUSION_CAPTURE_FORMAT,"invalid output entry");
 detail::need(cells>=1&&cells<=100000,FUSION_CAPTURE_FORMAT,"invalid output cells");
 const size_t kinetic=6*static_cast<size_t>(cells);
 detail::need(out.s.size()==kinetic&&out.t.size()==kinetic,
              FUSION_CAPTURE_FORMAT,"invalid output kinetic shape");
 if(entry==2)detail::need(out.mapped.size()==70*static_cast<size_t>(cells),
                         FUSION_CAPTURE_FORMAT,"invalid output packet shape");
 detail::HashSink h;
 // Canonical little-endian words: metadata, populations, result declaration
 // order, handoff, usage, outside, optional floor, optional packet rows/mapped.
 h.i32(entry);h.i32(cells);h.u32(floor?1:0);h.i32(out.status);
 h.doubles(out.n.data(),out.n.size());
 h.doubles(out.s.data(),out.s.size());h.doubles(out.t.data(),out.t.size());
 const auto &r=out.result;
 const auto &l=r.ledger;
 h.doubles(l.events_m3,5);
 h.doubles(l.nuclear_born_number_m3,6);h.doubles(l.nuclear_born_energy_J_m3,6);
 h.doubles(l.external_born_number_m3,6);h.doubles(l.external_born_energy_J_m3,6);
 h.doubles(l.thermal_consumed_number_m3,6);h.doubles(l.thermal_consumed_energy_J_m3,6);
 h.doubles(l.fast_consumed_number_m3,6);h.doubles(l.fast_consumed_energy_J_m3,6);
 h.doubles(l.escaped_number_m3,6);h.doubles(l.escaped_energy_J_m3,6);
 h.doubles(l.handed_off_number_m3,6);h.doubles(l.handed_off_energy_J_m3,6);
 h.doubles(l.heat_to_bath_J_m3,42);
 h.f64(l.neutron_number_m3);h.f64(l.neutron_energy_J_m3);
 h.doubles(r.inert_ion_heat_J_m3,6);
 h.f64(r.electron_energy_J_m3);h.f64(r.ion_energy_J_m3);
 h.doubles(r.particle_residual_m3,6);h.f64(r.energy_residual_J_m3);
 h.f64(r.max_source_rate_discrepancy);h.f64(r.max_source_debit_discrepancy);
 h.doubles(r.handoff_L1,6);h.doubles(r.handoff_mean_error,6);
 for(int v:r.handoff_projected)h.i32(v);
 const auto &d=out.handoff;
 h.doubles(d.candidate_number_m3,6);h.doubles(d.candidate_energy_J_m3,6);
 h.doubles(d.transferred_number_m3,6);h.doubles(d.transferred_energy_J_m3,6);
 h.doubles(d.target_kT_J,6);for(int v:d.tested)h.i32(v);
 const auto &u=out.usage;
 h.i32(u.direct_evaluations);h.i32(u.table_evaluations);
 h.f64(u.max_validated_rate_error);h.f64(u.max_validated_debit_error);
 h.f64(u.max_validated_number_L1);h.f64(u.max_validated_energy_L1);
 h.u64(out.outside);
 if(floor) {
  const auto &f=out.floor;
  h.doubles(f.born_number_m3,6);h.doubles(f.born_energy_J_m3,6);
  h.doubles(f.mapped_number_m3,6);h.doubles(f.mapped_energy_J_m3,6);
  h.doubles(f.ion_energy_correction_J_m3,6);h.doubles(f.energy_residual_J_m3,6);
  h.f64(f.remaining_ion_energy_J_m3);
 }
 if(entry==2) {
  const auto &p=out.packets;
  for(const auto &row:p.events_m3)h.doubles(row,5);
  for(const auto &source:p.below_number_m3)for(const auto &row:source)h.doubles(row,7);
  for(const auto &source:p.below_energy_J_m3)for(const auto &row:source)h.doubles(row,7);
  for(const auto &source:p.above_number_m3)for(const auto &row:source)h.doubles(row,7);
  for(const auto &source:p.above_energy_J_m3)for(const auto &row:source)h.doubles(row,7);
  h.doubles(out.mapped.data(),out.mapped.size());
 }
 return detail::hex(h.hash.digest());
}
}
