#include "fusion_transfer_prepared_internal.h"
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
    fusion_detail::prepared_transfer prepared;
    fusion_detail::prepare_transfer_validated(mass,charge,*bath,prepared);
    return fusion_detail::evaluate_transfer_prepared(energy,prepared,out);
}
