#include "fusion_coulomb.h"
#include <cmath>
#include <limits>

namespace {
using Real=long double;
constexpr Real electron_charge=1.602176634e-19L;
constexpr Real epsilon0=8.8541878188e-12L; // CODATA2022, F/m.
bool put(Real value,double &destination) {
    if (!std::isfinite(value) || std::abs(value)>std::numeric_limits<double>::max()) return false;
    destination=static_cast<double>(value);
    return value==0 || destination!=0;
}
}
extern "C" int fusion_c_coulomb_energy(double energy,double mass,double charge,
    const fusion_maxwellian_bath_v1 *bath,fusion_coulomb_energy_v1 *out) {
    if (!out) return PB11_STATUS_NULL_OUTPUT;
    *out={};
    if (!bath) return PB11_STATUS_INVALID_ARGUMENT;
    if (!std::isfinite(energy) || !std::isfinite(mass) || !std::isfinite(charge) ||
        !std::isfinite(bath->density_m3) || !std::isfinite(bath->mass_kg) ||
        !std::isfinite(bath->mean_charge_squared) || !std::isfinite(bath->kT_J) ||
        !std::isfinite(bath->coulomb_log)) return PB11_STATUS_INVALID_ARGUMENT;
    if (energy<0 || mass<=0 || bath->density_m3<0 || bath->mass_kg<=0 ||
        bath->mean_charge_squared<0 || bath->kT_J<=0 || bath->coulomb_log<=0)
        return PB11_STATUS_OUT_OF_RANGE;
    if (charge==0 || bath->mean_charge_squared==0 || bath->density_m3==0) return PB11_STATUS_OK;
    try {
        const Real pi=std::acos(-1.L),sqrtpi=std::sqrt(pi);
        const Real ma=mass,mb=bath->mass_kg,temp=bath->kT_J;
        const Real thermal_speed=std::sqrt(2*temp/mb);
        const Real y2=(static_cast<Real>(energy)/temp)*(mb/ma),y=std::sqrt(y2);
        Real h_over_y;
        if (y<.01L) {
            // H(y)/y, H=erf(y)-2y exp(-y^2)/sqrt(pi), stable at zero.
            h_over_y=4/sqrtpi*y2*(1.L/3-y2/5+y2*y2/14-y2*y2*y2/54+y2*y2*y2*y2/264);
        } else {
            h_over_y=(std::erf(y)-2*y*std::exp(-y2)/sqrtpi)/y;
        }
        const Real coulomb=electron_charge*electron_charge/(4*pi*epsilon0);
        const Real common=4*pi*static_cast<Real>(bath->density_m3)*charge*charge*
            bath->mean_charge_squared*coulomb*coulomb*bath->coulomb_log;
        const Real diffusion=common*temp/(mb*thermal_speed)*h_over_y;
        const Real drift=common/(ma*thermal_speed)*(2*std::exp(-y2)/sqrtpi-(ma/mb)*h_over_y);
        fusion_coulomb_energy_v1 result{};
        if (!put(diffusion,result.diffusion_J2_s) || !put(drift,result.mean_energy_rate_J_s))
            return PB11_STATUS_NUMERICAL_FAILURE;
        if (result.diffusion_J2_s<0) return PB11_STATUS_NUMERICAL_FAILURE;
        *out=result;
        return PB11_STATUS_OK;
    } catch (...) {return PB11_STATUS_EXCEPTION;}
}
