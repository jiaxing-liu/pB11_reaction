#ifndef PB11_SOURCE_H
#define PB11_SOURCE_H

#include "pb11_c.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Outputs for the simplified instantaneous-thermalization source model.
 *
 * This model supplies thermal p-11B reaction and power sources only.  It is
 * deliberately not a fast-alpha transport or fast-alpha energy-closure model.
 */
typedef struct pb11_instant_source_v1 {
    double reaction_rate_m3_s;
    double proton_source_m3_s;
    double boron_source_m3_s;
    double helium4_source_m3_s;
    double fusion_power_W_m3;
    double electron_power_W_m3;
    double ion_power_W_m3;
} pb11_instant_source_v1;

/*
 * Evaluate the simplified instantaneous thermal p-11B source in SI units.
 *
 * KT_J is kT in joules, NP_M3 and NB_M3 are number densities in m^-3, and
 * ELECTRON_FRACTION is the fraction of fusion power deposited to electrons.
 * METHOD is one of PB11_REACTIVITY_INTEGRAL or PB11_REACTIVITY_FAST.
 *
 * All fields in OUT are set to zero before validation and remain zero on any
 * failure.  The function returns one of the PB11_STATUS_* values from
 * pb11_c.h.  This is an instantaneous thermalization approximation, not a
 * fast-alpha closure.
 */
int pb11_c_instant_thermal_source(double kT_J, double np_m3, double nB_m3,
                                  double electron_fraction, int method,
                                  pb11_instant_source_v1* out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* PB11_SOURCE_H */
