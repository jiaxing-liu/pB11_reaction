#include "pb11_source.h"
#include "fusion_constants.hpp"

#include "pb11_c.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace {

constexpr double kJoulesPerKeV = 1.602176634e-16;
constexpr long double kJoulesPerKeVLong = 1.602176634e-16L;
constexpr long double kFusionEnergyPerReactionJ =
    fusion_constants::q_MeV[0] * fusion_constants::joules_per_MeV;

void clear_output(pb11_instant_source_v1* out) {
    out->reaction_rate_m3_s = 0.0;
    out->proton_source_m3_s = 0.0;
    out->boron_source_m3_s = 0.0;
    out->helium4_source_m3_s = 0.0;
    out->fusion_power_W_m3 = 0.0;
    out->electron_power_W_m3 = 0.0;
    out->ion_power_W_m3 = 0.0;
}

/*
 * Multiplying the smallest factor first avoids overflowing an intermediate
 * density product when the final three-factor product is representable.
 * long double also provides the extra exponent range needed for densities
 * whose pairwise product is above the double range.
 */
long double multiply_nonnegative(long double first, long double second,
                                long double third) {
    if (first == 0.0L || second == 0.0L || third == 0.0L) {
        return 0.0L;
    }

    std::array<long double, 3> factors{first, second, third};
    std::sort(factors.begin(), factors.end());

    long double product = 1.0L;
    const long double max_long_double =
        std::numeric_limits<long double>::max();
    for (const long double factor : factors) {
        if (product > max_long_double / factor) {
            return std::numeric_limits<long double>::infinity();
        }
        product *= factor;
    }
    return product;
}

/* Convert a computed value while rejecting both overflow and nonzero underflow. */
bool store_finite_double(long double value, double* destination) {
    if (!std::isfinite(value)) {
        return false;
    }
    if (value == 0.0L) {
        *destination = 0.0;
        return true;
    }

    const long double max_double =
        static_cast<long double>(std::numeric_limits<double>::max());
    if (std::fabs(value) > max_double) {
        return false;
    }

    const double converted = static_cast<double>(value);
    if (!std::isfinite(converted) || converted == 0.0) {
        return false;
    }
    *destination = converted;
    return true;
}

}  // namespace

extern "C" int pb11_c_instant_thermal_source(
    double kT_J, double np_m3, double nB_m3, double electron_fraction,
    int method, pb11_instant_source_v1* out) {
    if (out == nullptr) {
        return PB11_STATUS_NULL_OUTPUT;
    }

    clear_output(out);

    try {
        if (!std::isfinite(kT_J) || !std::isfinite(np_m3) ||
            !std::isfinite(nB_m3) || !std::isfinite(electron_fraction)) {
            return PB11_STATUS_INVALID_ARGUMENT;
        }
        if (kT_J <= 0.0 || np_m3 < 0.0 || nB_m3 < 0.0 ||
            electron_fraction < 0.0 || electron_fraction > 1.0) {
            return PB11_STATUS_OUT_OF_RANGE;
        }

        /* Check the SI-to-keV conversion before passing a double to the ABI. */
        const long double temperature_keV_long =
            static_cast<long double>(kT_J) / kJoulesPerKeVLong;
        const long double max_double =
            static_cast<long double>(std::numeric_limits<double>::max());
        if (!std::isfinite(temperature_keV_long) ||
            temperature_keV_long > max_double) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }

        double temperature_keV = kT_J / kJoulesPerKeV;
        if (!std::isfinite(temperature_keV) || temperature_keV <= 0.0) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }

        /*
         * A decimal SI literal for exactly 10 or 500 keV can round one ULP
         * above the endpoint when divided in double precision.  Canonicalize
         * that conversion-only ULP before the existing fast-domain check.
         */
        if (method == PB11_REACTIVITY_FAST) {
            const double fast_lower = 10.0;
            const double fast_upper = 500.0;
            if (temperature_keV > fast_upper &&
                temperature_keV <= std::nextafter(fast_upper,
                                                  std::numeric_limits<double>::infinity())) {
                temperature_keV = fast_upper;
            } else if (temperature_keV < fast_lower &&
                       temperature_keV >=
                           std::nextafter(fast_lower, 0.0)) {
                temperature_keV = fast_lower;
            }
        }

        double reactivity_m3_s = 0.0;
        const int reactivity_status =
            pb11_c_reactivity(temperature_keV, method, &reactivity_m3_s);
        if (reactivity_status != PB11_STATUS_OK) {
            return reactivity_status;
        }
        if (!std::isfinite(reactivity_m3_s) || reactivity_m3_s < 0.0) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }

        const long double reaction_rate = multiply_nonnegative(
            static_cast<long double>(np_m3), static_cast<long double>(nB_m3),
            static_cast<long double>(reactivity_m3_s));
        if (!std::isfinite(reaction_rate)) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }

        const long double helium_source = 3.0L * reaction_rate;
        const long double fusion_power =
            reaction_rate * kFusionEnergyPerReactionJ;
        const long double electron_power =
            fusion_power * static_cast<long double>(electron_fraction);
        const long double ion_power =
            fusion_power *
            (1.0L - static_cast<long double>(electron_fraction));
        if (!std::isfinite(helium_source) || !std::isfinite(fusion_power) ||
            !std::isfinite(electron_power) || !std::isfinite(ion_power)) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }

        pb11_instant_source_v1 result{};
        if (!store_finite_double(reaction_rate, &result.reaction_rate_m3_s) ||
            !store_finite_double(-reaction_rate,
                                 &result.proton_source_m3_s) ||
            !store_finite_double(-reaction_rate,
                                 &result.boron_source_m3_s) ||
            !store_finite_double(helium_source,
                                 &result.helium4_source_m3_s) ||
            !store_finite_double(fusion_power, &result.fusion_power_W_m3) ||
            !store_finite_double(electron_power,
                                 &result.electron_power_W_m3) ||
            !store_finite_double(ion_power, &result.ion_power_W_m3)) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }

        *out = result;
        return PB11_STATUS_OK;
    } catch (...) {
        clear_output(out);
        return PB11_STATUS_EXCEPTION;
    }
}
