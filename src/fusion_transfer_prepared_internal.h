#ifndef FUSION_TRANSFER_PREPARED_INTERNAL_H
#define FUSION_TRANSFER_PREPARED_INTERNAL_H
#include "fusion_two_component.h"
#include <cmath>
#include <limits>
namespace fusion_detail {
struct prepared_transfer {
    long double mass, temperature, bath_mass;
    long double normalization;
    bool zero;
};
// Caller has validated finite/domain inputs, exactly as the scalar ABI does.
inline void prepare_transfer_validated(double mass,double charge,
    const fusion_maxwellian_bath_v1& bath,prepared_transfer& p) {
    p.mass=mass;p.temperature=bath.kT_J;p.bath_mass=bath.mass_kg;
    p.zero=charge==0 || bath.mean_charge_squared==0 || bath.density_m3==0;
    if(p.zero){p.normalization=0;return;}
    using R=long double;
    const R pi=std::acos(-1.L),e=1.602176634e-19L,eps=8.8541878188e-12L;
    const R log_vti=(std::log(2.L)+std::log(R(bath.kT_J))-
                      std::log(R(bath.mass_kg)))/2;
    const R log_nu=std::log(4*pi)+std::log(R(bath.density_m3))+
        2*std::log(std::abs(R(charge)))+std::log(R(bath.mean_charge_squared))+
        2*(2*std::log(e)-std::log(4*pi*eps))+std::log(R(bath.coulomb_log))-
        std::log(R(mass))-std::log(R(bath.mass_kg));
    p.normalization=std::log(4/std::sqrt(pi))+log_nu-3*log_vti;
}
inline int evaluate_transfer_prepared(double energy,const prepared_transfer& p,double* out) {
    *out=0;
    if(p.zero)return PB11_STATUS_OK;
    using R=long double;
    const R y2=(R(energy)/p.temperature)*(p.bath_mass/p.mass);
    const R log_rate=p.normalization-y2;
    if (std::isnan(log_rate) || log_rate>std::log(R(std::numeric_limits<double>::max())))
        return PB11_STATUS_NUMERICAL_FAILURE;
    const R rate=std::exp(log_rate);
    if (!std::isfinite(rate) || rate<0) return PB11_STATUS_NUMERICAL_FAILURE;
    *out=static_cast<double>(rate);
    return PB11_STATUS_OK;
}
}
#endif
