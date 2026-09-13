#include "fusion_two_component.h"
#include <cmath>
#include <limits>

extern "C" int fusion_c_coulomb_transfer_rate(double energy,double mass,
    double charge,const fusion_maxwellian_bath_v1 *bath,double *out) {
    if (!out) return PB11_STATUS_NULL_OUTPUT;
    *out=0;
    if (!bath) return PB11_STATUS_INVALID_ARGUMENT;
    if (!std::isfinite(energy) || !std::isfinite(mass) || !std::isfinite(charge) ||
        !std::isfinite(bath->density_m3) || !std::isfinite(bath->mass_kg) ||
        !std::isfinite(bath->mean_charge_squared) || !std::isfinite(bath->kT_J) ||
        !std::isfinite(bath->coulomb_log)) return PB11_STATUS_INVALID_ARGUMENT;
    if (energy<0 || mass<=0 || bath->density_m3<0 || bath->mass_kg<=0 ||
        bath->mean_charge_squared<0 || bath->kT_J<=0 || bath->coulomb_log<=0)
        return PB11_STATUS_OUT_OF_RANGE;
    if (charge==0 || bath->mean_charge_squared==0 || bath->density_m3==0)
        return PB11_STATUS_OK;
    // Logarithmic evaluation avoids spurious intermediate overflow at an
    // exponentially small high-energy tail. Long double retains the exponent
    // range of ratios of finite positive double inputs on supported platforms.
    using R=long double;
    const R pi=std::acos(-1.L),e=1.602176634e-19L,eps=8.8541878188e-12L;
    const R log_vti=(std::log(2.L)+std::log(R(bath->kT_J))-
                      std::log(R(bath->mass_kg)))/2;
    const R log_nu=std::log(4*pi)+std::log(R(bath->density_m3))+
        2*std::log(std::abs(R(charge)))+std::log(R(bath->mean_charge_squared))+
        2*(2*std::log(e)-std::log(4*pi*eps))+std::log(R(bath->coulomb_log))-
        std::log(R(mass))-std::log(R(bath->mass_kg));
    const R y2=(R(energy)/bath->kT_J)*(R(bath->mass_kg)/mass);
    const R log_rate=std::log(4/std::sqrt(pi))+log_nu-3*log_vti-y2;
    if (std::isnan(log_rate) || log_rate>std::log(R(std::numeric_limits<double>::max())))
        return PB11_STATUS_NUMERICAL_FAILURE;
    // -infinity can describe a vanishing Maxwellian tail on platforms with
    // no extended long-double exponent range. This is a physical underflow.
    const R rate=std::exp(log_rate);
    if (!std::isfinite(rate) || rate<0) return PB11_STATUS_NUMERICAL_FAILURE;
    *out=static_cast<double>(rate);
    return PB11_STATUS_OK;
}
