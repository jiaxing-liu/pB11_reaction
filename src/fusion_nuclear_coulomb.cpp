#include "fusion_nuclear_coulomb.h"
#include "fusion_coulomb_tables.inc"
#include <cmath>
namespace {
constexpr double mev=1.602176634e-13;
double cheb(const double *c,double x) {
    double b1=0,b2=0;
    for(int k=16;k>=1;--k){double b=2*x*b1-b2+c[k];b2=b1;b1=b;}
    return x*b1-b2+c[0];
}
}
extern "C" int fusion_c_nuclear_coulomb(int channel,double input,
    fusion_nuclear_coulomb_v1 *out) {
    if(!out)return PB11_STATUS_NULL_OUTPUT;
    *out={};
    if(channel<0 || channel>4 || !std::isfinite(input))return PB11_STATUS_INVALID_ARGUMENT;
    if(input<.001*mev || input>12*mev)return PB11_STATUS_OUT_OF_RANGE;
    double E=input/mev,x=std::log(E);
    // Rounding the conversion/log at inclusive limits is within one ulp.
    const double low=std::log(.001),high=std::log(12.);
    if(x<low)x=low;
    if(x>high)x=high;
    for(const auto &segment:fusion_table_data::segments) {
        if(segment.channel!=channel || x<segment.lo || x>segment.hi)continue;
        const double t=(2*x-segment.lo-segment.hi)/(segment.hi-segment.lo);
        fusion_nuclear_coulomb_v1 result{};
        result.log_penetrability=cheb(segment.coeff[0],t);
        result.shift=cheb(segment.coeff[1],t);
        result.phase_real=cheb(segment.coeff[2],t);
        result.phase_imag=cheb(segment.coeff[3],t);
        const double mu=3727.3794118*(channel<3?2./3:.5),a=channel<3?5.1:4.5;
        result.rho=a*std::sqrt(2*mu*E)/197.3269804;
        const double norm=std::hypot(result.phase_real,result.phase_imag);
        if(!std::isfinite(result.log_penetrability) || !std::isfinite(result.shift) ||
           !std::isfinite(norm) || std::abs(norm-1)>1e-6)return PB11_STATUS_NUMERICAL_FAILURE;
        // Preserve the evaluated interpolation error; do not renormalize phase.
        *out=result;return PB11_STATUS_OK;
    }
    return PB11_STATUS_NUMERICAL_FAILURE;
}
