#include "fusion_target_burn.h"
#include <math.h>
#include <stdio.h>

_Static_assert(sizeof(int) == 4, "C int ABI must be four bytes");
_Static_assert(sizeof(fusion_target_burn_v1) == 14 * sizeof(double),
               "target-burn v1 ABI must contain fourteen doubles");

int main(void) {
    const double kev = 1.602176634e-16;
    const double edges[2] = {0.0, 2.0};
    const double old_fast[1] = {1.0};
    const double K[1] = {1.0};
    const double M[1] = {0.25};
    double trial[1] = {-1.0}, loss[1] = {-1.0};
    fusion_target_burn_v1 out;
    const int status = fusion_c_target_burn_trial(
        1, 1.0, edges, old_fast, 1.0, 1.0, K, M, trial, loss, &out);
    const double exact = (sqrt(5.0) - 1.0) / 2.0;
    if (status != PB11_STATUS_OK || fabs(trial[0] - exact) > 5e-14 ||
        fabs(loss[0] - (1.0 - exact)) > 5e-14 ||
        out.reactions_m3 < 0.0) {
        fprintf(stderr, "target-burn trial smoke failed: status=%d trial=%.17g\n",
                status, trial[0]);
        return 1;
    }

    const double beam_edges[2] = {400.0 * kev, 600.0 * kev};
    const double beam_old[1] = {1.0e15};
    double beam_trial[1], beam_loss[1];
    fusion_beam_window_v1 window[1];
    fusion_target_burn_v1 beam_out;
    const int beam_status = fusion_c_beam_target_burn_window_trial(
        0, 1, 1.0, 1.67262192595e-27, 11.0 * 1.66053906892e-27,
        beam_edges, beam_old, 1.0e20, 0.0, beam_trial, beam_loss,
        window, &beam_out);
    if (beam_status != PB11_STATUS_OK ||
        window[0].resolved_reactivity_m3_s <= 0.0 ||
        window[0].domain_incomplete != 0 || beam_out.reactions_m3 <= 0.0) {
        fprintf(stderr, "window smoke failed: status=%d rate=%.17g\n",
                beam_status, window[0].resolved_reactivity_m3_s);
        return 1;
    }

    printf("PASS trial=%.17g window_rate=%.17g reactions=%.17g\n",
           trial[0], window[0].resolved_reactivity_m3_s,
           beam_out.reactions_m3);
    return 0;
}
