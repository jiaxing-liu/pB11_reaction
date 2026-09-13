#ifndef FUSION_CONSTANTS_HPP
#define FUSION_CONSTANTS_HPP
namespace fusion_constants {
inline constexpr long double joules_per_MeV=1.602176634e-13L;
// CODATA 2022 nuclear mass-energy equivalents, NIST allascii.txt (MeV).
inline constexpr long double p=938.27208943L, d=1875.61294500L,
    t=2808.92113668L, he3=2808.39161112L, he4=3727.3794118L,
    neutron=939.56542194L;
// pB preserves the explicitly rounded 8.68 MeV baseline. Other Q values
// follow nuclear (not neutral-atom) masses; no electron binding correction.
inline constexpr long double q_MeV[5]={8.68L,2*d-t-p,2*d-he3-neutron,
                                      d+t-he4-neutron,d+he3-he4-p};
}
#endif
