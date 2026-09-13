#include "fusion_alpha_spectrum.h"
#include <math.h>
#include <stdio.h>

_Static_assert(sizeof(int) == 4, "C int ABI must be four bytes");
_Static_assert(sizeof(fusion_alpha_spectrum_v1) == 96,
               "alpha spectrum v1 ABI must be 96 bytes");

int main(void) {
    const double mev = 1.602176634e-13;
    const double available = 9.3 * mev;
    const double edges[3] = {2.0 * mev, 3.0 * mev, 4.0 * mev};
    double birth[2] = {-1.0, -1.0};
    fusion_alpha_spectrum_v1 out;
    const int status = fusion_c_alpha_spectrum_grid(
        2, available, 0.001 * mev, 0.5, 0.0, 16, 16, 2,
        edges, birth, &out);
    const double number = out.mapped_number + out.below_number + out.above_number;
    const double energy = out.mapped_energy_J + out.below_energy_J +
                          out.above_energy_J;
    if (status != PB11_STATUS_OK || !isfinite(number) || !isfinite(energy) ||
        fabs(number - 3.0) > 4e-12 || fabs(energy - available) > 3e-11 * available ||
        out.below_number <= 0.0 || out.above_number <= 0.0 ||
        birth[0] < 0.0 || birth[1] < 0.0) {
        fprintf(stderr, "alpha-grid smoke failed: status=%d N=%.17g E=%.17g\n",
                status, number, energy);
        return 1;
    }
    printf("PASS status=%d N=%.17g E=%.17g below=%.17g above=%.17g\n",
           status, number, energy, out.below_number, out.above_number);
    return 0;
}
