#include "fusion_network.h"

#include <cmath>
#include <limits>

namespace {
// Rows: pB, DD->T+p, DD->He3+n, DT, DHe3. Columns: p,D,T,He3,He4,B11.
constexpr int consumed[5][6] = {
    {1,0,0,0,0,1}, {0,2,0,0,0,0}, {0,2,0,0,0,0},
    {0,1,1,0,0,0}, {0,1,0,1,0,0}
};
constexpr int born[5][6] = {
    {0,0,0,0,3,0}, {1,0,1,0,0,0}, {0,0,0,1,0,0},
    {0,0,0,0,1,0}, {1,0,0,0,1,0}
};
constexpr int neutrons[5] = {0,0,1,1,0};
bool representable(long double value) {
    return std::isfinite(value) &&
        std::fabs(value) <= std::numeric_limits<double>::max();
}
}

extern "C" int fusion_c_particle_sources(const double *rates,
                                           fusion_particle_sources_v1 *out) {
    if (!out) return PB11_STATUS_NULL_OUTPUT;
    *out = {};
    if (!rates) return PB11_STATUS_INVALID_ARGUMENT;
    try {
        for (int c=0; c<FUSION_CHANNEL_COUNT; ++c) {
            if (!std::isfinite(rates[c])) return PB11_STATUS_INVALID_ARGUMENT;
            if (rates[c]<0) return PB11_STATUS_OUT_OF_RANGE;
        }
        fusion_particle_sources_v1 result{};
        for (int s=0; s<FUSION_SPECIES_COUNT; ++s) {
            long double loss=0, birth=0;
            for (int c=0; c<FUSION_CHANNEL_COUNT; ++c) {
                loss += static_cast<long double>(rates[c])*consumed[c][s];
                birth += static_cast<long double>(rates[c])*born[c][s];
            }
            if (!representable(loss) || !representable(birth) ||
                !representable(birth-loss)) return PB11_STATUS_NUMERICAL_FAILURE;
            result.reactant_loss[s]=static_cast<double>(loss);
            result.product_birth[s]=static_cast<double>(birth);
            result.net_source[s]=static_cast<double>(birth-loss);
        }
        long double nbirth=0;
        for (int c=0; c<FUSION_CHANNEL_COUNT; ++c)
            nbirth += static_cast<long double>(rates[c])*neutrons[c];
        if (!representable(nbirth)) return PB11_STATUS_NUMERICAL_FAILURE;
        result.neutron_birth=static_cast<double>(nbirth);
        *out=result;
        return PB11_STATUS_OK;
    } catch (...) {
        return PB11_STATUS_EXCEPTION;
    }
}
