#ifndef FUSION_BIRTH_TABLE_INTERNAL_H
#define FUSION_BIRTH_TABLE_INTERNAL_H
#include "fusion_birth_table.h"
#include "fusion_thermal_birth_budget.h"
namespace fusion_detail {
bool birth_table_matches(const fusion_birth_table_v1*,int channel,
 const fusion_thermal_birth_options_v1&,int cells,const double* edges) noexcept;
bool birth_table_budget_matches(const fusion_birth_table_v1*,int channel,
 const fusion_thermal_birth_options_v1&,const fusion_thermal_birth_budget_policy_v1&,
 int cells,const double*edges) noexcept;
}
#endif
