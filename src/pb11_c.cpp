#include "pb11_c.h"

#include "pb11.hpp"

#include <cmath>
#include <exception>

namespace {

template <typename Validator, typename Evaluation>
int checked_call(double input, double* result, Validator validator,
                 Evaluation evaluation) noexcept {
    if (result == nullptr) {
        return PB11_STATUS_NULL_OUTPUT;
    }

    // Keeping the error value finite prevents an unchecked foreign caller
    // from accidentally injecting NaN into a solver state.
    *result = 0.0;

    if (!std::isfinite(input)) {
        return PB11_STATUS_INVALID_ARGUMENT;
    }
    if (!validator(input)) {
        return PB11_STATUS_OUT_OF_RANGE;
    }

    try {
        const double value = evaluation(input);
        if (!std::isfinite(value)) {
            return PB11_STATUS_NUMERICAL_FAILURE;
        }
        *result = value;
        return PB11_STATUS_OK;
    } catch (const std::exception&) {
        return PB11_STATUS_EXCEPTION;
    } catch (...) {
        return PB11_STATUS_EXCEPTION;
    }
}

bool fitted_energy(double energy_MeV) {
    return energy_MeV >= 0.0 && energy_MeV <= 9.76;
}

bool positive_temperature(double temperature_keV) {
    return temperature_keV > 0.0;
}

bool fast_temperature(double temperature_keV) {
    return temperature_keV >= 10.0 && temperature_keV <= 500.0;
}

}  // namespace

extern "C" int pb11_c_abi_version(void) {
    return 1;
}

extern "C" const char* pb11_c_status_message(int status) {
    switch (status) {
    case PB11_STATUS_OK:
        return "ok";
    case PB11_STATUS_NULL_OUTPUT:
        return "output pointer is null";
    case PB11_STATUS_INVALID_ARGUMENT:
        return "input is not finite";
    case PB11_STATUS_OUT_OF_RANGE:
        return "input is outside the published range";
    case PB11_STATUS_NUMERICAL_FAILURE:
        return "underlying calculation returned a non-finite value";
    case PB11_STATUS_EXCEPTION:
        return "underlying C++ calculation threw an exception";
    case PB11_STATUS_UNKNOWN_METHOD:
        return "unknown reactivity method";
    default:
        return "unknown status";
    }
}

extern "C" int pb11_c_sfactor(double E_MeV, double* result) {
    return checked_call(E_MeV, result, fitted_energy,
                        [](double energy) { return pb11_sfactor(energy); });
}

extern "C" int pb11_c_cross_section(double E_MeV, double* result) {
    return checked_call(
        E_MeV, result, fitted_energy,
        [](double energy) { return pb11_cross_section(energy); });
}

extern "C" int pb11_c_reactivity_integral(double T_keV, double* result) {
    return checked_call(
        T_keV, result, positive_temperature,
        [](double temperature) { return pb11_reactivity_integral(temperature); });
}

extern "C" int pb11_c_reactivity_fast(double T_keV, double* result) {
    return checked_call(
        T_keV, result, fast_temperature,
        [](double temperature) { return pb11_reactivity_fast(temperature); });
}

extern "C" int pb11_c_reactivity(double T_keV, int method, double* result) {
    if (result == nullptr) {
        return PB11_STATUS_NULL_OUTPUT;
    }
    *result = 0.0;

    switch (method) {
    case PB11_REACTIVITY_INTEGRAL:
        return pb11_c_reactivity_integral(T_keV, result);
    case PB11_REACTIVITY_FAST:
        return pb11_c_reactivity_fast(T_keV, result);
    default:
        return PB11_STATUS_UNKNOWN_METHOD;
    }
}
