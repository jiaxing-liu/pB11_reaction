#ifndef FUSION_NETWORK_H
#define FUSION_NETWORK_H

#include "pb11_c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Stable IDs; no relationship to a host's species-array ordering. */
enum fusion_species_v1 {
    FUSION_PROTON = 0, FUSION_DEUTERON = 1, FUSION_TRITON = 2,
    FUSION_HELIUM3 = 3, FUSION_HELIUM4 = 4, FUSION_BORON11 = 5,
    FUSION_SPECIES_COUNT = 6
};
enum fusion_channel_v1 {
    FUSION_PB11_3ALPHA = 0, FUSION_DD_TP = 1, FUSION_DD_HE3N = 2,
    FUSION_DT_ALPHAN = 3, FUSION_DHE3_ALPHAP = 4,
    FUSION_CHANNEL_COUNT = 5
};

typedef struct fusion_particle_sources_v1 {
    /* All entries have units m^-3 s^-1. Loss and birth are nonnegative. */
    double reactant_loss[FUSION_SPECIES_COUNT];
    double product_birth[FUSION_SPECIES_COUNT];
    double net_source[FUSION_SPECIES_COUNT];
    double neutron_birth;
} fusion_particle_sources_v1;

/*
 * Assemble nuclear stoichiometry from reaction EVENT rates, m^-3 s^-1.
 * DD inputs are branch event rates and must already include the identical
 * pair factor 1/2. This function never inserts a second pair factor.
 *
 * Births are nuclear products, NOT necessarily thermalized host particles.
 * No slowing, deposition, transport, density or temperature evolution occurs.
 * Neutrons are recorded separately and must enter the caller's escape ledger.
 *
 * All input rates must be finite and nonnegative. Arrays have the fixed
 * lengths defined above. Returns PB11_STATUS_*; output is zero on failure.
 * Input and output memory must not overlap. NULL input is INVALID_ARGUMENT;
 * NULL output is NULL_OUTPUT. No C++ exception crosses this ABI.
 */
int fusion_c_particle_sources(const double event_rates[FUSION_CHANNEL_COUNT],
                              fusion_particle_sources_v1 *out);

#ifdef __cplusplus
}
#endif
#endif
