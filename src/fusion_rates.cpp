#include "fusion_rates.h"
#include "fusion_constants.hpp"
#include <cmath>
#include <limits>

namespace {
constexpr double joules_per_keV=1.602176634e-16;
struct Fit {
    double bg, reduced_mass_keV;
    double c[7];
    double minimum_keV, maximum_keV;
};
// Bosch & Hale, Nuclear Fusion 32 (1992) 611, Table VII (p625).
// Internal order follows the public channels excluding pB: DD(Tp),DD(He3n),DT,DHe3.
constexpr Fit fits[] = {
    {31.3970,937814.,{5.65718e-12,3.41267e-3,1.99167e-3,0.,1.05060e-5,0.,0.},.2,100.},
    {31.3970,937814.,{5.43360e-12,5.85778e-3,7.68222e-3,0.,-2.96400e-6,0.,0.},.2,100.},
    {34.3827,1124656.,{1.17302e-9,1.51361e-2,7.51886e-2,4.60643e-3,1.35000e-2,-1.06750e-4,1.36600e-5},.2,100.},
    {68.7508,1124572.,{5.51036e-10,6.41918e-3,-2.02896e-3,-1.91080e-5,1.35776e-4,0.,0.},.5,190.}
};
constexpr int reactants[5][2]={{0,5},{1,1},{1,1},{1,2},{1,3}};
bool method_valid(int method) {
    return method==PB11_REACTIVITY_INTEGRAL || method==PB11_REACTIVITY_FAST;
}
bool domain(double &temperature, double lower, double upper) {
    if (temperature<lower && temperature>=std::nextafter(lower,0.)) temperature=lower;
    if (temperature>upper && temperature<=std::nextafter(upper,INFINITY)) temperature=upper;
    return temperature>=lower && temperature<=upper;
}
}

extern "C" int fusion_c_channel_q(int channel,double *out) {
    if (!out) return PB11_STATUS_NULL_OUTPUT;
    *out=0;
    if (channel<0 || channel>=FUSION_CHANNEL_COUNT) return PB11_STATUS_INVALID_ARGUMENT;
    *out=static_cast<double>(fusion_constants::q_MeV[channel]*fusion_constants::joules_per_MeV);
    return PB11_STATUS_OK;
}

extern "C" int fusion_c_thermal_reactivity(int channel, double kT_J,
                                             int pb_method, double *out) {
    if (!out) return PB11_STATUS_NULL_OUTPUT;
    *out=0;
    if (!std::isfinite(kT_J)) return PB11_STATUS_INVALID_ARGUMENT;
    if (kT_J<=0) return PB11_STATUS_OUT_OF_RANGE;
    if (channel<0 || channel>=FUSION_CHANNEL_COUNT) return PB11_STATUS_INVALID_ARGUMENT;
    if (!method_valid(pb_method)) return PB11_STATUS_UNKNOWN_METHOD;
    const double tk=kT_J/joules_per_keV;
    if (!std::isfinite(tk) || tk<=0) return PB11_STATUS_NUMERICAL_FAILURE;
    try {
        double temperature=tk;
        if (channel==FUSION_PB11_3ALPHA) {
            if (pb_method==PB11_REACTIVITY_FAST && !domain(temperature,10.,500.))
                return PB11_STATUS_OUT_OF_RANGE;
            return pb11_c_reactivity(temperature,pb_method,out);
        }
        const auto &fit=fits[channel-1];
        if (!domain(temperature,fit.minimum_keV,fit.maximum_keV)) return PB11_STATUS_OUT_OF_RANGE;
        const auto &c=fit.c;
        // Equations (12)-(14). Published result is cm^3/s; convert once to SI.
        const double theta=temperature/(1-temperature*(c[1]+temperature*(c[3]+temperature*c[5]))/
            (1+temperature*(c[2]+temperature*(c[4]+temperature*c[6]))));
        const double xi=std::cbrt(fit.bg*fit.bg/(4*theta));
        const double result=c[0]*theta*std::sqrt(xi/(fit.reduced_mass_keV*temperature*temperature*temperature))*
            std::exp(-3*xi)*1e-6;
        if (!std::isfinite(result) || result<0) return PB11_STATUS_NUMERICAL_FAILURE;
        *out=result;
        return PB11_STATUS_OK;
    } catch (...) {
        *out=0;
        return PB11_STATUS_EXCEPTION;
    }
}

extern "C" int fusion_c_thermal_rates(double kT_J, const double *density,
                                        int mask, int pb_method,
                                        fusion_thermal_rates_v1 *out) {
    if (!out) return PB11_STATUS_NULL_OUTPUT;
    *out={};
    if (!density || !std::isfinite(kT_J)) return PB11_STATUS_INVALID_ARGUMENT;
    if (kT_J<=0) return PB11_STATUS_OUT_OF_RANGE;
    if (mask<0 || (static_cast<unsigned int>(mask) & ~((1u<<FUSION_CHANNEL_COUNT)-1u))) return PB11_STATUS_INVALID_ARGUMENT;
    if (!method_valid(pb_method)) return PB11_STATUS_UNKNOWN_METHOD;
    for (int s=0;s<FUSION_SPECIES_COUNT;++s) {
        if (!std::isfinite(density[s])) return PB11_STATUS_INVALID_ARGUMENT;
        if (density[s]<0) return PB11_STATUS_OUT_OF_RANGE;
    }
    try {
        fusion_thermal_rates_v1 result{};
        long double total_power=0;
        for (int ch=0;ch<FUSION_CHANNEL_COUNT;++ch) {
            if (!(mask & (1u<<ch))) continue;
            const int status=fusion_c_thermal_reactivity(ch,kT_J,pb_method,&result.reactivity_m3_s[ch]);
            if (status!=PB11_STATUS_OK) return status;
            const int a=reactants[ch][0], b=reactants[ch][1];
            // Multiply by the small rate first; extended range avoids naive n*n overflow.
            const long double rate=static_cast<long double>(result.reactivity_m3_s[ch]) *
                density[a]*density[b]*(a==b ? .5L : 1.L);
            if (!std::isfinite(rate) || rate>std::numeric_limits<double>::max())
                return PB11_STATUS_NUMERICAL_FAILURE;
            result.event_rate_m3_s[ch]=static_cast<double>(rate);
            if (rate!=0 && result.event_rate_m3_s[ch]==0) return PB11_STATUS_NUMERICAL_FAILURE;
            const long double power=rate*fusion_constants::q_MeV[ch]*fusion_constants::joules_per_MeV;
            if (!std::isfinite(power) || power>std::numeric_limits<double>::max())
                return PB11_STATUS_NUMERICAL_FAILURE;
            result.nuclear_power_W_m3[ch]=static_cast<double>(power);
            if (power!=0 && result.nuclear_power_W_m3[ch]==0) return PB11_STATUS_NUMERICAL_FAILURE;
            total_power+=power;
        }
        if (!std::isfinite(total_power) || total_power>std::numeric_limits<double>::max())
            return PB11_STATUS_NUMERICAL_FAILURE;
        result.total_nuclear_power_W_m3=static_cast<double>(total_power);
        const int status=fusion_c_particle_sources(result.event_rate_m3_s,&result.particles);
        if (status!=PB11_STATUS_OK) return status;
        *out=result;
        return PB11_STATUS_OK;
    } catch (...) {
        return PB11_STATUS_EXCEPTION;
    }
}
