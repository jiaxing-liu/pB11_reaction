#ifndef PB11_C_H
#define PB11_C_H

/*
 * Stable C ABI for the p-11B reaction library.
 *
 * The C++ API in pb11.hpp remains available for native C++ callers.  These
 * entry points are intended for Fortran and other foreign-function callers:
 * they never let a C++ exception cross the ABI boundary and never return a
 * non-finite value on an error path.  The result is set to zero when the
 * function returns a non-zero status.
 */

#ifdef __cplusplus
extern "C" {
#endif

enum pb11_status_code {
    PB11_STATUS_OK = 0,
    PB11_STATUS_NULL_OUTPUT = 1,
    PB11_STATUS_INVALID_ARGUMENT = 2,
    PB11_STATUS_OUT_OF_RANGE = 3,
    PB11_STATUS_NUMERICAL_FAILURE = 4,
    PB11_STATUS_EXCEPTION = 5,
    PB11_STATUS_UNKNOWN_METHOD = 6
};

enum pb11_reactivity_method {
    PB11_REACTIVITY_INTEGRAL = 0,
    PB11_REACTIVITY_FAST = 1
};

/* Returns the ABI revision.  It is currently 1. */
int pb11_c_abi_version(void);

/* Return a static human-readable string for a status code. */
const char* pb11_c_status_message(int status);

/*
 * Inputs and outputs use the same units as the native C++ API:
 *   energy E: MeV, S-factor: MeV barn, cross section: barn,
 *   temperature kT: keV, reactivity: m^3/s.
 */
int pb11_c_sfactor(double E_MeV, double* result);
int pb11_c_cross_section(double E_MeV, double* result);
int pb11_c_reactivity_integral(double T_keV, double* result);
int pb11_c_reactivity_fast(double T_keV, double* result);
int pb11_c_reactivity(double T_keV, int method, double* result);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif  /* PB11_C_H */
