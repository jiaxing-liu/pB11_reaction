#ifndef FUSION_BEAM_BIRTH_TABLE_INTERNAL_H
#define FUSION_BEAM_BIRTH_TABLE_INTERNAL_H
#include "fusion_beam_birth_table.h"
namespace fusion_detail {
bool beam_birth_table_matches(const fusion_beam_birth_table_v1*, int channel,
 int slot, double projectile_energy_J, const fusion_beam_birth_options_v1&,
 int cells, const double* edges) noexcept;
}
#endif
