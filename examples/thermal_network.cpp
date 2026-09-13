#include "fusion_rates.h"
#include <iomanip>
#include <iostream>

int main() {
    // SI example: thermal DD fuel with explicitly enabled secondary DT/DHe3.
    // Secondary rates start at zero because no T or He3 is initially supplied.
    double densities[FUSION_SPECIES_COUNT]={0,1e20,0,0,0,0};
    fusion_thermal_rates_v1 result{};
    const int status=fusion_c_thermal_rates(10*1.602176634e-16,densities,30,
                                            PB11_REACTIVITY_INTEGRAL,&result);
    if (status) {std::cerr<<pb11_c_status_message(status)<<'\n'; return 1;}
    std::cout<<std::scientific<<std::setprecision(8);
    for (int c=0;c<FUSION_CHANNEL_COUNT;++c)
        std::cout<<"channel "<<c<<" events [m^-3 s^-1] "<<result.event_rate_m3_s[c]
                 <<" nuclear release [W m^-3] "<<result.nuclear_power_W_m3[c]<<'\n';
    std::cout<<"neutron births [m^-3 s^-1] "<<result.particles.neutron_birth<<'\n';
    std::cout<<"Product births above are not yet thermalized or deposited.\n";
}
