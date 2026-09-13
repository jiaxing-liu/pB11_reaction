#include "fusion_nuclear_data.h"
#include "fusion_rates.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using R = long double;
constexpr R c = 299792458.L, c2 = c * c;
constexpr R e = 1.602176634e-19L, mev = 1.602176634e-13L;
constexpr R u = 1.66053906892e-27L;

void check(bool ok, const char *what) { if (!ok) throw std::runtime_error(what); }
bool close(R a, R b, R tol = 3e-12L, R floor = 0.L) {
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    const R scale = std::max({std::abs(a), std::abs(b), floor});
    return std::abs(a - b) <= tol * scale;
}
bool zero_mass(const fusion_nuclear_mass_v1 &x) {
    return x.mass_kg == 0 && x.rest_energy_J == 0 && x.reference_mass_u == 0 &&
           x.known_uncertainty_scale_kg == 0 && x.nuclear_charge == 0 &&
           x.mass_number == 0;
}
void clear_pair(fusion_particle_four_vector_v1 p[2], double value) {
    for (int i = 0; i < 2; ++i) {
        p[i].mass_kg = value; p[i].kinetic_energy_J = value;
        for (double &x : p[i].momentum_kg_m_s) x = value;
    }
}
bool zero_pair(const fusion_particle_four_vector_v1 p[2]) {
    for (int i = 0; i < 2; ++i) {
        if (p[i].mass_kg != 0 || p[i].kinetic_energy_J != 0) return false;
        for (double x : p[i].momentum_kg_m_s) if (x != 0) return false;
    }
    return true;
}

void mass_and_channels(std::array<fusion_nuclear_mass_v1, 8> &m) {
    constexpr int charge[8] = {1, 1, 1, 2, 2, 5, 0, -1};
    constexpr int number[8] = {1, 2, 3, 3, 4, 11, 1, 0};
    for (int id = 0; id < 8; ++id) {
        check(fusion_c_nuclear_mass(id, &m[id]) == PB11_STATUS_OK, "valid mass ID");
        check(std::isfinite(m[id].mass_kg) && m[id].mass_kg > 0 &&
              std::isfinite(m[id].rest_energy_J) && m[id].rest_energy_J > 0 &&
              std::isfinite(m[id].reference_mass_u) && m[id].reference_mass_u > 0 &&
              std::isfinite(m[id].known_uncertainty_scale_kg) &&
              m[id].known_uncertainty_scale_kg >= 0, "finite mass record");
        check(m[id].nuclear_charge == charge[id] && m[id].mass_number == number[id],
              "mass charge/number bookkeeping");
        check(close(m[id].rest_energy_J, R(m[id].mass_kg) * c2), "mass rest energy");
    }
    const R me = 5.485799090441e-4L * u;
    check(close(m[7].mass_kg, me, 5e-14L), "CODATA electron mass");
    const R expected_b = 11.009305166L * u - 5.L * me +
                         670.9838405L * e / c2;
    check(close(m[5].mass_kg, expected_b, 5e-14L), "AME/NIST B11 nuclear mass");
    for (int id : {-1, 8, 99}) {
        fusion_nuclear_mass_v1 bad{1, 2, 3, 4, 5, 6};
        check(fusion_c_nuclear_mass(id, &bad) == PB11_STATUS_OUT_OF_RANGE &&
              zero_mass(bad), "invalid mass ID clears output");
    }
}

void channels_and_q(const std::array<fusion_nuclear_mass_v1, 8> &m) {
    constexpr int react[5][2] = {{0,5}, {1,1}, {1,1}, {1,2}, {1,3}};
    constexpr int prod[5][3] = {{4,4,4}, {2,0,-1}, {3,6,-1}, {4,6,-1}, {4,0,-1}};
    const R golden_q[5]={8.682378274225927L,4.032663898739255L,3.268856950827066L,
                         17.589247903737996L,18.353054851650186L};
    const R golden_unc_eV[5]={12.293309118248754L,.130088468564021L,.470496675140556L,
                              .542980194030242L,.154134294060342L};
    for (int ch = 0; ch < 5; ++ch) {
        fusion_nuclear_channel_v1 r{};
        check(fusion_c_nuclear_channel(ch, &r) == PB11_STATUS_OK, "valid channel");
        check(r.product_count == (ch == 0 ? 3 : 2), "channel product count");
        R q = 0;
        for (int j = 0; j < 2; ++j) {
            check(r.reactant_ids[j] == react[ch][j], "channel reactant IDs");
            q += m[react[ch][j]].mass_kg;
        }
        for (int j = 0; j < 3; ++j) {
            check(r.product_ids[j] == prod[ch][j], "channel product IDs");
            if (j < r.product_count) q -= m[prod[ch][j]].mass_kg;
        }
        check(close(r.q_J, q * c2, 2e-11L), "Q equals returned mass difference");
        check(close(r.q_J/mev,golden_q[ch],2e-11L),"independent decimal Q reference");
        check(close(r.known_uncertainty_scale_J/e,golden_unc_eV[ch],2e-10L),
              "independent linear uncertainty reference");
        int charge=0,baryon=0;
        for(int id:r.reactant_ids){charge+=m[id].nuclear_charge;baryon+=m[id].mass_number;}
        for(int j=0;j<r.product_count;++j){int id=r.product_ids[j];
            charge-=m[id].nuclear_charge;baryon-=m[id].mass_number;}
        check(charge==0&&baryon==0,"reaction charge and baryon conservation");
        check(std::isfinite(r.known_uncertainty_scale_J) &&
              r.known_uncertainty_scale_J >= 0, "channel uncertainty record");
    }
    fusion_nuclear_channel_v1 bad{1, 2, {1,1}, {1,1,1}, 3};
    check(fusion_c_nuclear_channel(-1, &bad) == PB11_STATUS_OUT_OF_RANGE &&
          bad.q_J == 0 && bad.known_uncertainty_scale_J == 0 &&
          bad.product_count == 0, "invalid channel clears output");
    double old_q = -1;
    check(fusion_c_channel_q(FUSION_PB11_3ALPHA, &old_q) == PB11_STATUS_OK &&
          close(old_q, 8.68L * mev, 2e-14L), "legacy pB Q remains 8.68 MeV");
    fusion_nuclear_channel_v1 pb{};
    check(fusion_c_nuclear_channel(FUSION_PB11_3ALPHA, &pb) == PB11_STATUS_OK &&
          close(pb.q_J, 8.6823782742L * mev, 3e-10L), "new pB Q value");
}

void two_body(const std::array<fusion_nuclear_mass_v1, 8> &m) {
    constexpr int product[5][2] = {{-1,-1}, {2,0}, {3,6}, {4,6}, {4,0}};
    const double d[3] = {0.6, -0.8, 0.0};
    for (int ch = 1; ch < 5; ++ch) {
        fusion_nuclear_channel_v1 r{};
        check(fusion_c_nuclear_channel(ch, &r) == PB11_STATUS_OK, "two-body channel");
        for (double E_MeV : {0.0, 0.01, 1.0}) {
            fusion_particle_four_vector_v1 p[2];
            check(fusion_c_nuclear_two_body_cm(ch, E_MeV * double(mev), d, p) == PB11_STATUS_OK,
                  "two-body status");
            const R available = E_MeV * mev + r.q_J;
            check(p[0].mass_kg == m[product[ch][0]].mass_kg &&
                  p[1].mass_kg == m[product[ch][1]].mass_kg, "two-body product masses");
            check(std::isfinite(p[0].kinetic_energy_J) && p[0].kinetic_energy_J >= 0 &&
                  std::isfinite(p[1].kinetic_energy_J) && p[1].kinetic_energy_J >= 0 &&
                  close(R(p[0].kinetic_energy_J) + p[1].kinetic_energy_J, available),
                  "product kinetic energy is Ecm plus Q");
            R p2[3]{};
            for (int j = 0; j < 3; ++j) {
                p2[j] = R(p[0].momentum_kg_m_s[j]) + p[1].momentum_kg_m_s[j];
                check(close(p2[j], 0, 3e-12L, 1e-40L), "CM momenta are opposite");
            }
            const R initial = (R(m[r.reactant_ids[0]].mass_kg) +
                               R(m[r.reactant_ids[1]].mass_kg)) * c2 + E_MeV * mev;
            const R final = (R(p[0].mass_kg) + R(p[1].mass_kg)) * c2 +
                            p[0].kinetic_energy_J + p[1].kinetic_energy_J;
            check(close(initial, final, 5e-12L), "total rest plus kinetic energy");
            if (E_MeV > 0) {
                const R pm = std::sqrt(R(p[0].momentum_kg_m_s[0]) * p[0].momentum_kg_m_s[0] +
                                       R(p[0].momentum_kg_m_s[1]) * p[0].momentum_kg_m_s[1]);
                check(pm > 0 && close(R(p[0].momentum_kg_m_s[0]), .6L * pm, 5e-12L) &&
                      close(R(p[0].momentum_kg_m_s[1]), -.8L * pm, 5e-12L),
                      "product-0 momentum follows direction");
            }
        }
    }
    fusion_particle_four_vector_v1 p[2];
    clear_pair(p, -1);
    check(fusion_c_nuclear_two_body_cm(FUSION_PB11_3ALPHA, 0, d, p) ==
          PB11_STATUS_OUT_OF_RANGE && zero_pair(p), "pB is not a two-body channel");
    double nan = std::numeric_limits<double>::quiet_NaN();
    clear_pair(p, -1);
    check(fusion_c_nuclear_two_body_cm(FUSION_DT_ALPHAN, nan, d, p) ==
          PB11_STATUS_INVALID_ARGUMENT && zero_pair(p), "NaN energy is rejected");
    const double bad[3] = {1, 1, 0};
    clear_pair(p, -1);
    check(fusion_c_nuclear_two_body_cm(FUSION_DT_ALPHAN, 0, bad, p) ==
          PB11_STATUS_OUT_OF_RANGE && zero_pair(p), "bad direction is rejected");
    const double nan_dir[3] = {nan, 0, 1};
    clear_pair(p, -1);
    check(fusion_c_nuclear_two_body_cm(FUSION_DT_ALPHAN, 0, nan_dir, p) ==
          PB11_STATUS_INVALID_ARGUMENT && zero_pair(p), "NaN direction is rejected");
}
}

int main() {
    try {
        std::array<fusion_nuclear_mass_v1, 8> masses{};
        mass_and_channels(masses);
        channels_and_q(masses);
        two_body(masses);
        std::cout << "All nuclear-data tests passed\n";
        return 0;
    } catch (const std::exception &ex) {
        std::cerr << "FAIL: " << ex.what() << '\n';
        return 1;
    }
}
