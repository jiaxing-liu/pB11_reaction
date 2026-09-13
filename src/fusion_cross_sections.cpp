#include "fusion_cross_sections.h"
#include <cmath>
#include <limits>

namespace {
constexpr double keV_J=1.602176634e-16;
struct CrossFit {
    double bg;
    double a[5],b[4];
    double minimum,maximum;
};
// Bosch-Hale 1992 Table IV. Output from Eqs8/9 is millibarn.
constexpr CrossFit low[] = {
    {31.3970,{5.5576e4,2.1054e2,-3.2638e-2,1.4987e-6,1.8181e-10},{0,0,0,0},.5,5000},
    {31.3970,{5.3701e4,3.3027e2,-1.2706e-1,2.9327e-5,-2.5151e-9},{0,0,0,0},.5,4900},
    {34.3827,{6.927e4,7.454e8,2.050e6,5.2002e4,0},{6.38e1,-9.95e-1,6.981e-5,1.728e-4},.5,4700},
    {68.7508,{5.7501e6,2.5226e3,4.5566e1,0,0},{-3.1995e-3,-8.5530e-6,5.9014e-8,0},.3,4800}
};
// Table VI contains only a constant numerator, not inherited low-energy A2..A5.
constexpr CrossFit high[] = {
    {34.3827,{-1.4714e6,0,0,0,0},{-8.4127e-3,4.7983e-6,-1.0748e-9,8.5184e-14},530,4700},
    {68.7508,{-8.3993e5,0,0,0,0},{-2.6830e-3,1.1633e-6,-2.1332e-10,1.4250e-14},900,4800}
};
}
extern "C" int fusion_c_cross_section(int ch,double energy_J,double *out) {
    if (!out) return PB11_STATUS_NULL_OUTPUT;
    *out=0;
    if (!std::isfinite(energy_J) || ch<0 || ch>=FUSION_CHANNEL_COUNT)
        return PB11_STATUS_INVALID_ARGUMENT;
    if (energy_J<0) return PB11_STATUS_OUT_OF_RANGE;
    if (energy_J==0) return PB11_STATUS_OK;
    double energy=energy_J/keV_J;
    if (!std::isfinite(energy)) return PB11_STATUS_NUMERICAL_FAILURE;
    try {
        if (ch==0) {
            double barn=0;
            double mev=energy/1000.;
            if (mev>9.76 && mev<=std::nextafter(9.76,INFINITY)) mev=9.76;
            const int status=pb11_c_cross_section(mev,&barn);
            if (status) return status;
            *out=barn*1e-28;
            return PB11_STATUS_OK;
        }
        const auto &limits=low[ch-1];
        if (energy<limits.minimum && energy>=std::nextafter(limits.minimum,0.)) energy=limits.minimum;
        if (energy>limits.maximum && energy<=std::nextafter(limits.maximum,INFINITY)) energy=limits.maximum;
        if (energy<limits.minimum || energy>limits.maximum) return PB11_STATUS_OUT_OF_RANGE;
        const CrossFit *fit=&limits;
        if (ch==FUSION_DT_ALPHAN && energy>=530) fit=&high[0];
        if (ch==FUSION_DHE3_ALPHAP && energy>=900) fit=&high[1];
        const auto &a=fit->a; const auto &b=fit->b;
        const double s=(a[0]+energy*(a[1]+energy*(a[2]+energy*(a[3]+energy*a[4]))))/
            (1+energy*(b[0]+energy*(b[1]+energy*(b[2]+energy*b[3]))));
        const double sigma=s/energy*std::exp(-fit->bg/std::sqrt(energy))*1e-31;
        if (!std::isfinite(sigma) || sigma<0) return PB11_STATUS_NUMERICAL_FAILURE;
        *out=sigma;
        return PB11_STATUS_OK;
    } catch (...) {
        *out=0;
        return PB11_STATUS_EXCEPTION;
    }
}
