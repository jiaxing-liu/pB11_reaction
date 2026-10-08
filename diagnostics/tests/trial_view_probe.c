/* ABI-only probe: never invokes numerical or capture/replay routines. */
#include "fusion_capture_trial.h"
#include <stddef.h>
#include <string.h>

int trial_view_probe(const fusion_capture_trial_input_v1 *view,
                     const void *const *pointers,
                     size_t *view_size,
                     fusion_capture_replay_report_v1 *report,
                     size_t *report_size)
{
    *view_size = sizeof(*view);
    *report_size = sizeof(*report);
    if (view->entry_point != 101) return 1;
    if (view->original_status != 102) return 2;
    if (view->host_zone != 103) return 3;
    if (view->host_time_s != 4.25) return 4;
    if (view->dt_s != 5.25) return 5;
    if ((const void *)view->options != pointers[0]) return 6;
    if ((const void *)view->fast_options != pointers[1]) return 7;
    if ((const void *)view->thermal_tables != pointers[2]) return 8;
    if (view->beam_table_count != 109) return 9;
    if ((const void *)view->beam_tables != pointers[3]) return 10;
    if (view->effective_charge != 111) return 11;
    if (view->cells != 112) return 12;
    if ((const void *)view->edges_J != pointers[4]) return 13;
    if ((const void *)view->thermal_number_m3 != pointers[5]) return 14;
    if (view->electron_energy_J_m3 != 15.25) return 15;
    if (view->ion_energy_J_m3 != 16.25) return 16;
    if (view->electron_density_m3 != 17.25) return 17;
    if ((const void *)view->thermal_charge_squared != pointers[6]) return 18;
    if (view->inert_count != 119) return 19;
    if ((const void *)view->inert != pointers[7]) return 20;
    if ((const void *)view->coulomb_logs != pointers[8]) return 21;
    if ((const void *)view->old_s_m3 != pointers[9]) return 22;
    if ((const void *)view->old_t_m3 != pointers[10]) return 23;
    if ((const void *)view->external_birth_m3_s != pointers[11]) return 24;
    if ((const void *)view->escape_s_inv != pointers[12]) return 25;
    if ((const void *)view->floor_limits != pointers[13]) return 26;
    if (view->table_domain_policy != 127) return 27;
    memset(report, 0, sizeof(*report));
    report->original_status = 201;
    report->replay_status = 202;
    report->status_matches = 203;
    report->entry_point = 204;
    report->host_zone = 205;
    report->cells = 206;
    report->host_time_s = 7.75;
    report->dt_s = 8.75;
    memset(report->kernel_identity, 'K', 64);
    memset(report->output_sha256, 'a', 64);
    return 0;
}
