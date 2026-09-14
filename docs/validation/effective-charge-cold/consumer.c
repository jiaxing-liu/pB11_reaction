#include "fusion_coupled_thermal.h"
int main(void) {
 int status=fusion_c_coupled_thermal_table_trial_effective_charge(
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
 return status==PB11_STATUS_NULL_OUTPUT?0:1;
}
