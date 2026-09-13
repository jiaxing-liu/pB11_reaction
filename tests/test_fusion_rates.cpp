#include "fusion_rates.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {
int failures=0;
void check(bool value,const char *name) {
    if (!value) { ++failures; std::cerr<<"FAIL: "<<name<<'\n'; }
}
bool close(double x,double y,double relative=1e-12) {
    return std::isfinite(x) && std::isfinite(y) &&
        std::abs(x-y)<=relative*std::max(std::abs(x),std::abs(y));
}
constexpr double kev=1.602176634e-16;
}
int main() {
    // Independent published values: Bosch-Hale Table VIII p625, units cm^3/s.
    // Columns: DT,DHe3,DD(Tp),DD(He3n). These are not computed by this test.
    constexpr double temperatures[]={1,5,10,20,50};
    constexpr double reference[][4]={
        {6.857e-21,3.057e-26,1.017e-22,9.933e-23},
        {1.366e-17,6.377e-21,9.024e-20,9.128e-20},
        {1.136e-16,2.126e-19,5.781e-19,6.023e-19},
        {4.330e-16,3.482e-18,2.399e-18,2.603e-18},
        {8.649e-16,5.554e-17,9.838e-18,1.133e-17}
    };
    constexpr int channels[]={FUSION_DT_ALPHAN,FUSION_DHE3_ALPHAP,FUSION_DD_TP,FUSION_DD_HE3N};
    constexpr double tolerance[]={.0032,.026,.004,.0035};
    for (int row=0;row<5;++row) for (int col=0;col<4;++col) {
        double rate=-1;
        check(fusion_c_thermal_reactivity(channels[col],temperatures[row]*kev,0,&rate)==0,"published rate status");
        check(close(rate,reference[row][col]*1e-6,tolerance[col]),"Table VIII reactivity including fit/rounding error");
    }
    for (int ch=1;ch<5;++ch) {
        const double low=ch==4?.5:.2, high=ch==4?190.:100.;
        double rate;
        check(fusion_c_thermal_reactivity(ch,low*kev,0,&rate)==0 && rate>0,"lower temperature endpoint");
        check(fusion_c_thermal_reactivity(ch,high*kev,0,&rate)==0 && rate>0,"upper temperature endpoint");
        check(fusion_c_thermal_reactivity(ch,(low-.001)*kev,0,&rate)==PB11_STATUS_OUT_OF_RANGE && rate==0,"below fit domain");
        check(fusion_c_thermal_reactivity(ch,(high+.001)*kev,0,&rate)==PB11_STATUS_OUT_OF_RANGE && rate==0,"above fit domain");
    }
    // Independent decimal subtraction of CODATA nuclear rest energies (MeV).
    constexpr double q_mev[]={8.68,4.03266389,3.26885694,17.58924794,18.35305489};
    for (int ch=0;ch<5;++ch) {
        double q;
        check(fusion_c_channel_q(ch,&q)==0 && close(q,q_mev[ch]*1.602176634e-13),"nuclear Q values and joule conversion");
    }
    double rate, old;
    for (int method=0;method<2;++method) {
        check(fusion_c_thermal_reactivity(0,100*kev,method,&rate)==0,"pB scalar status");
        check(pb11_c_reactivity(100,method,&old)==0 && close(rate,old),"pB scalar compatibility");
    }
    check(fusion_c_thermal_reactivity(99,10*kev,0,&rate)==PB11_STATUS_INVALID_ARGUMENT && rate==0,"bad channel");
    check(fusion_c_thermal_reactivity(3,10*kev,99,&rate)==PB11_STATUS_UNKNOWN_METHOD && rate==0,"bad method");
    check(fusion_c_thermal_reactivity(3,10*kev,0,nullptr)==PB11_STATUS_NULL_OUTPUT,"null rate pointer");
    double densities[6]={0,1e20,0,0,0,0};
    fusion_thermal_rates_v1 result{};
    const unsigned dd=(1u<<FUSION_DD_TP)|(1u<<FUSION_DD_HE3N);
    check(fusion_c_thermal_rates(10*kev,densities,dd,0,&result)==0,"DD network status");
    const double rt=result.event_rate_m3_s[1], rn=result.event_rate_m3_s[2];
    check(close(rt,.5e40*result.reactivity_m3_s[1]),"DD Tp factor one half");
    check(close(rn,.5e40*result.reactivity_m3_s[2]),"DD He3n factor one half");
    check(close(result.particles.reactant_loss[1],2*(rt+rn)),"DD consumes two deuterons");
    check(close(result.particles.product_birth[2],rt) && close(result.particles.product_birth[3],rn),"DD secondary reactants born");
    check(result.event_rate_m3_s[3]==0 && result.event_rate_m3_s[4]==0,"secondary channels off initially");
    densities[2]=1e18; densities[3]=2e18;
    check(fusion_c_thermal_rates(10*kev,densities,30u,0,&result)==0,"secondary network status");
    check(result.event_rate_m3_s[3]>0 && result.event_rate_m3_s[4]>0,"secondary burn responds to T and He3");
    check(result.particles.product_birth[4]>0 && result.particles.reactant_loss[2]>0 && result.particles.reactant_loss[3]>0,"secondary alpha birth and fuel depletion");
    check(fusion_c_thermal_rates(10*kev,densities,0u,0,&result)==0 && result.particles.neutron_birth==0,"all channels disabled");
    check(fusion_c_thermal_rates(10*kev,densities,32u,0,&result)==PB11_STATUS_INVALID_ARGUMENT && result.particles.product_birth[4]==0,"invalid mask clears output");
    check(fusion_c_thermal_rates(10*kev,densities,-1,0,&result)==PB11_STATUS_INVALID_ARGUMENT,"negative mask");
    densities[1]=0;
    check(fusion_c_thermal_rates(.1*kev,densities,dd,0,&result)==PB11_STATUS_OUT_OF_RANGE,"zero fuel does not hide invalid active domain");
    densities[1]=-1;
    check(fusion_c_thermal_rates(10*kev,densities,dd,0,&result)==PB11_STATUS_OUT_OF_RANGE,"negative density");
    densities[1]=std::numeric_limits<double>::quiet_NaN();
    check(fusion_c_thermal_rates(10*kev,densities,0u,0,&result)==PB11_STATUS_INVALID_ARGUMENT,"disabled channels still validate density");
    densities[1]=std::numeric_limits<double>::max();
    check(fusion_c_thermal_rates(10*kev,densities,dd,0,&result)==PB11_STATUS_NUMERICAL_FAILURE && result.event_rate_m3_s[1]==0,"event output overflow");
    check(fusion_c_thermal_rates(10*kev,nullptr,dd,0,&result)==PB11_STATUS_INVALID_ARGUMENT,"null density");
    check(fusion_c_thermal_rates(10*kev,densities,dd,0,nullptr)==PB11_STATUS_NULL_OUTPUT,"null network output");
    if (failures) return 1;
    std::cout<<"PASS: published thermal rates, domains, SI, DD symmetry, secondary channels and ABI errors\n";
}
