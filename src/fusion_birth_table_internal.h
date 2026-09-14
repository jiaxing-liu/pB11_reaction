#ifndef FUSION_BIRTH_TABLE_INTERNAL_H
#define FUSION_BIRTH_TABLE_INTERNAL_H
#include "fusion_birth_table.h"
namespace fusion_detail {
bool birth_table_matches(const fusion_birth_table_v1*,int channel,
 const fusion_thermal_birth_options_v1&,int cells,const double* edges) noexcept;
}
#endif
