#ifndef FUSION_CROSS_SECTIONS_H
#define FUSION_CROSS_SECTIONS_H
#include "fusion_network.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Centre-of-mass relative energy in joules; output cross section in m^2.
 * Bosch-Hale 1992 Tables IV/VI: DD Tp 0.5..5000 keV, DD He3n
 * 0.5..4900 keV, DT 0.5..4700 keV, DHe3 0.3..4800 keV.
 * DT switches at 530 keV following p622 text; DHe3 at 900 keV,
 * where the published piecewise fit has a small discontinuity.
 * pB uses existing Tentori-Belloni cross section, 0..9.76 MeV.
 * Zero energy returns the exact limiting value zero for all channels.
 * Positive energy below a channel's lower data limit is rejected, not clipped.
 * Status uses PB11_STATUS_*; outputs are zero on errors.
 */
int fusion_c_cross_section(int channel, double relative_energy_J,
                            double *cross_section_m2);
#ifdef __cplusplus
}
#endif
#endif
