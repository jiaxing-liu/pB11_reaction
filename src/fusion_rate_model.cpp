#include "fusion_rate_model.h"
#include <cmath>
#include <limits>

extern "C" int fusion_c_cross_section_model(int ch,int policy,int low,
    double e,double *out) {
    if(!out)return PB11_STATUS_NULL_OUTPUT;
    *out=0;
    if(!std::isfinite(e))return PB11_STATUS_INVALID_ARGUMENT;
    if(e<0)return PB11_STATUS_OUT_OF_RANGE;
    double emin=0,emax=0;
    int status=fusion_c_cross_section_domain(ch,&emin,&emax);
    if(status)return status;
    if((policy!=FUSION_ENDPOINT_S && policy!=FUSION_HIGH_FLAT) ||
       low<FUSION_PB_LOW_TB || low>FUSION_PB_LOW_C0_PLUS12 ||
       (ch!=FUSION_PB11_3ALPHA && low!=FUSION_PB_LOW_TB))
        return PB11_STATUS_INVALID_ARGUMENT;
    if(e==0)return PB11_STATUS_OK;
    constexpr long double kev=1.602176634e-16L,barn=1e-28L;
    const long double E=e/kev;
    long double sigma=0;
    if(ch==FUSION_PB11_3ALPHA && low!=FUSION_PB_LOW_TB && e<=400*double(kev)) {
        const long double c0=197+(low==FUSION_PB_LOW_C0_MINUS12?-12:
                                  low==FUSION_PB_LOW_C0_PLUS12?12:0);
        const long double c1=low==FUSION_PB_LOW_NS?.240L:.269L;
        const long double c2=low==FUSION_PB_LOW_NS?2.31e-4L:2.54e-4L;
        const long double d=E-148;
        const long double S=c0+c1*E+c2*E*E+1.82e4L/(d*d+2.35L*2.35L);
        sigma=barn*std::exp(std::log(1000*S)-std::log(E)-std::sqrt(22589.L/E));
    } else if(e>=emin && e<=emax) {
        return fusion_c_cross_section(ch,e,out);
    } else {
        const double boundary=e<emin?emin:emax;
        double endpoint=0;
        status=fusion_c_cross_section(ch,boundary,&endpoint);
        if(status)return status;
        if(e>emax && policy==FUSION_HIGH_FLAT) sigma=endpoint;
        else {
            const long double BG[]={std::sqrt(22589.L),31.3970L,31.3970L,34.3827L,68.7508L};
            const long double Eb=boundary/kev;
            sigma=std::exp(std::log(static_cast<long double>(endpoint))+
                std::log(Eb)-std::log(E)+BG[ch]*(1/std::sqrt(Eb)-1/std::sqrt(E)));
        }
    }
    if(!std::isfinite(sigma) || sigma<0 || sigma>std::numeric_limits<double>::max())
        return PB11_STATUS_NUMERICAL_FAILURE;
    *out=static_cast<double>(sigma);
    return PB11_STATUS_OK;
}
