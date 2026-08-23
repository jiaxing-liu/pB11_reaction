#include "pb11.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <array>
#include <cmath>
#include <limits>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

constexpr double kEnergyBoundary1MeV = 0.400;
constexpr double kEnergyBoundary2MeV = 0.668;
constexpr double kEnergyMaxMeV = 9.76;
constexpr double kGamowEnergyMeV = 22.589;
constexpr double kReducedMassEnergyMeV = 859.526;

// Exact SI definitions: 1 eV = 1.602176634e-19 J, 1 barn = 1e-28 m^2,
// and c = 299792458 m/s.  The reduced mass is obtained from mu = (mu c^2)/c^2.
constexpr double kMeVToJoule = 1.602176634e-13;
constexpr double kKeVToJoule = 1.602176634e-16;
constexpr double kBarnToSquareMetre = 1.0e-28;
constexpr double kSpeedOfLightMetresPerSecond = 299792458.0;
constexpr double kReducedMassKg =
    kReducedMassEnergyMeV * kMeVToJoule /
    (kSpeedOfLightMetresPerSecond * kSpeedOfLightMetresPerSecond);

// Table 1, "This work", plus the narrow-resonance parameters explicitly
// retained from Nevins & Swain as stated below equation (5).
constexpr double kC0 = 197.0;
constexpr double kC1 = 0.269;
constexpr double kC2 = 2.54e-4;
constexpr double kNarrowAmplitudeMeVBarn = 1.82e4;
constexpr double kNarrowEnergyKeV = 148.0;
constexpr double kNarrowWidthKeV = 2.35;

constexpr double kD0 = 346.0;
constexpr double kD1 = 150.0;
constexpr double kD2 = -59.9;
constexpr double kD5 = -0.460;

constexpr double kBackgroundMeVBarn = 0.381;
constexpr std::array<double, 4> kAmplitudesMeVBarn{
    1.98e6, 3.89e6, 1.36e6, 3.71e6};
constexpr std::array<double, 4> kResonanceEnergiesKeV{
    640.9, 1211.0, 2340.0, 3294.0};
constexpr std::array<double, 4> kResonanceWidthsKeV{
    85.5, 414.0, 221.0, 351.0};

double quiet_nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool valid_fitted_energy(double energy_MeV) {
    return std::isfinite(energy_MeV) && energy_MeV >= 0.0 &&
           energy_MeV <= kEnergyMaxMeV;
}

double sfactor_low(double energy_keV) {
    const double narrow_offset = energy_keV - kNarrowEnergyKeV;
    const double narrow =
        kNarrowAmplitudeMeVBarn /
        (narrow_offset * narrow_offset +
         kNarrowWidthKeV * kNarrowWidthKeV);
    return kC0 + kC1 * energy_keV +
           kC2 * energy_keV * energy_keV + narrow;
}

double sfactor_middle(double energy_MeV) {
    const double x = (energy_MeV - kEnergyBoundary1MeV) / 0.100;
    const double x2 = x * x;
    return kD0 + kD1 * x + kD2 * x2 + kD5 * x2 * x2 * x;
}

double sfactor_high(double energy_keV) {
    double result = kBackgroundMeVBarn;
    for (std::size_t i = 0; i < kAmplitudesMeVBarn.size(); ++i) {
        const double offset = energy_keV - kResonanceEnergiesKeV[i];
        const double width = kResonanceWidthsKeV[i];
        result += kAmplitudesMeVBarn[i] /
                  (offset * offset + width * width);
    }
    return result;
}

double integrate_cross_section_kernel(double temperature_MeV) {
    // Splitting out a symmetric interval around the 148 keV resonance makes
    // its peak a Gauss-Kronrod node.  The remaining breakpoints cover the
    // published piece boundaries and the centres of the broad resonances.
    constexpr double narrow_half_interval_MeV =
        20.0 * kNarrowWidthKeV / 1000.0;
    constexpr std::array<double, 10> breakpoints{
        0.0,
        kNarrowEnergyKeV / 1000.0 - narrow_half_interval_MeV,
        kNarrowEnergyKeV / 1000.0 + narrow_half_interval_MeV,
        kEnergyBoundary1MeV,
        kEnergyBoundary2MeV,
        1.211,
        2.340,
        3.294,
        5.700,
        kEnergyMaxMeV};

    const auto kernel = [temperature_MeV](double energy_MeV) {
        if (energy_MeV == 0.0) {
            return 0.0;
        }
        return energy_MeV * pb11_cross_section(energy_MeV) *
               std::exp(-energy_MeV / temperature_MeV);
    };

    using Integrator = boost::math::quadrature::gauss_kronrod<double, 61>;
    constexpr unsigned max_depth = 15;
    constexpr double relative_tolerance = 1.0e-10;

    double result = 0.0;
    for (std::size_t i = 1; i < breakpoints.size(); ++i) {
        result += Integrator::integrate(
            kernel, breakpoints[i - 1], breakpoints[i], max_depth,
            relative_tolerance);
    }
    return result;
}

double low_temperature_reactivity(double temperature_keV) {
    // Equations (10)-(14).  Energies used in dimensionless ratios remain in
    // keV; the dimensional factors are converted to SI before equation (10).
    constexpr double gamow_energy_keV = 1000.0 * kGamowEnergyMeV;
    const double effective_energy_keV =
        std::cbrt(gamow_energy_keV * temperature_keV * temperature_keV / 4.0);
    const double effective_width_keV =
        4.0 * std::sqrt(temperature_keV * effective_energy_keV / 3.0);
    const double tau = 3.0 * effective_energy_keV / temperature_keV;

    const double effective_sfactor_MeVBarn =
        kC0 * (1.0 + 5.0 / (12.0 * tau)) +
        kC1 * (effective_energy_keV + (35.0 / 36.0) * temperature_keV) +
        kC2 * (effective_energy_keV * effective_energy_keV +
               (89.0 / 36.0) * effective_energy_keV * temperature_keV);

    const double temperature_J = temperature_keV * kKeVToJoule;
    const double effective_width_J = effective_width_keV * kKeVToJoule;
    const double effective_sfactor_Jm2 =
        effective_sfactor_MeVBarn * kMeVToJoule * kBarnToSquareMetre;
    const double nonresonant =
        std::sqrt(2.0 / kReducedMassKg) * effective_width_J *
        effective_sfactor_Jm2 /
        std::pow(temperature_J, 1.5) * std::exp(-tau);

    // Equation (15).  A_L/delta_E_L has units of area: MeV barn / keV.
    const double resonance_area_m2 =
        (1000.0 * kNarrowAmplitudeMeVBarn / kNarrowWidthKeV) *
        kBarnToSquareMetre;
    const double thermal_speed =
        std::sqrt(8.0 * kPi * temperature_J / kReducedMassKg);
    const double resonant =
        thermal_speed * (1.0 / (temperature_keV * temperature_keV)) *
        resonance_area_m2 *
        std::exp(-std::sqrt(gamow_energy_keV / kNarrowEnergyKeV) -
                 kNarrowEnergyKeV / temperature_keV);

    return nonresonant + resonant;
}

double high_temperature_reactivity(double temperature_keV) {
    // Table 2, "This work".  With kT, theta and mu*c^2 expressed in keV,
    // equation (7) and the tabulated unit of P1 produce m^3/s directly.
    constexpr double p1_keVm3PerSecond = 9.7827e-14;
    constexpr double p2_perKeV = -5.1610e-2;
    constexpr double p3_perKeV = 1.3240e-1;
    constexpr double p4_perKeV2 = 3.8446e-4;
    constexpr double p5_perKeV2 = 1.2499e-3;
    constexpr double p6_perKeV3 = -5.6715e-6;
    constexpr double p7_perKeV3 = 1.1615e-6;
    constexpr double gamow_energy_keV = 1000.0 * kGamowEnergyMeV;
    constexpr double reduced_mass_energy_keV =
        1000.0 * kReducedMassEnergyMeV;

    const double numerator =
        p2_perKeV +
        temperature_keV *
            (p4_perKeV2 + temperature_keV * p6_perKeV3);
    const double denominator =
        1.0 + temperature_keV *
                  (p3_perKeV +
                   temperature_keV *
                       (p5_perKeV2 + temperature_keV * p7_perKeV3));
    const double theta_keV =
        temperature_keV /
        (1.0 - temperature_keV * numerator / denominator);
    const double xi = std::cbrt(gamow_energy_keV / (4.0 * theta_keV));

    return p1_keVm3PerSecond * theta_keV *
           std::sqrt(xi /
                     (reduced_mass_energy_keV *
                      temperature_keV * temperature_keV * temperature_keV)) *
           std::exp(-3.0 * xi);
}

}  // namespace

double pb11_sfactor(double E_MeV) {
    if (!valid_fitted_energy(E_MeV)) {
        return quiet_nan();
    }
    if (E_MeV <= kEnergyBoundary1MeV) {
        return sfactor_low(1000.0 * E_MeV);
    }
    if (E_MeV <= kEnergyBoundary2MeV) {
        return sfactor_middle(E_MeV);
    }
    return sfactor_high(1000.0 * E_MeV);
}

double pb11_cross_section(double E_MeV) {
    if (!valid_fitted_energy(E_MeV)) {
        return quiet_nan();
    }
    if (E_MeV == 0.0) {
        return 0.0;
    }

    // Writing (S/E)*exp(-sqrt(E_G/E)) in log form avoids inf*0 for tiny E.
    const double exponent =
        -std::sqrt(kGamowEnergyMeV / E_MeV) - std::log(E_MeV);
    return pb11_sfactor(E_MeV) * std::exp(exponent);
}

double pb11_reactivity_integral(double T_keV) {
    if (!std::isfinite(T_keV) || T_keV <= 0.0) {
        return quiet_nan();
    }

    const double temperature_MeV = T_keV / 1000.0;
    const double integral_MeV2Barn =
        integrate_cross_section_kernel(temperature_MeV);

    // In equation (6), E*sigma(E)*dE is converted from MeV^2 barn to
    // J^2 m^2.  Combined with mu in kg and kT in J, the result is m^3/s.
    const double integral_J2m2 =
        integral_MeV2Barn * kMeVToJoule * kMeVToJoule *
        kBarnToSquareMetre;
    const double temperature_J = T_keV * kKeVToJoule;
    const double prefactor =
        std::sqrt(8.0 / (kPi * kReducedMassKg)) /
        std::pow(temperature_J, 1.5);
    return prefactor * integral_J2m2;
}

double pb11_reactivity_fast(double T_keV) {
    if (!std::isfinite(T_keV) || T_keV < 10.0 || T_keV > 500.0) {
        return quiet_nan();
    }
    if (T_keV < 70.0) {
        return low_temperature_reactivity(T_keV);
    }
    return high_temperature_reactivity(T_keV);
}
