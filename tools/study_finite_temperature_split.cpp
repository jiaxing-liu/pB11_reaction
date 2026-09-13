#include "fusion_two_component.h"
// Finite-temperature two-component trace-alpha diagnostic.
//
// This is a numerical prototype only.  ST and TH are two populations on the
// same energy grid.  The ST escape and TH birth below are an internal linear
// transfer used to factor one full backward-Euler step; they are not a
// physical escape boundary or an instantaneous fluid ash model.
#include "fusion_coulomb.h"
#include "fusion_kinetics.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kJoulesPerKeV = 1.602176634e-16;
constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kEpsilon0 = 8.8541878188e-12;
constexpr double kAlphaMass = 6.644657345e-27;
constexpr double kElectronMass = 9.1093837139e-31;
constexpr double kDeuteronMass = 3.3435837768e-27;
constexpr double kTritonMass = 5.0073567512e-27;
constexpr double kPopulation = 1.0e12;
constexpr double kInitialAlphaEnergy = 3000.0 * kJoulesPerKeV;
constexpr double kBathTemperature = 10.0 * kJoulesPerKeV;
constexpr double kCoulombLog = 15.0;
constexpr double kDuration = 1.0;

struct Metrics {
    double st_number_fraction = 0.0;
    double th_number_fraction = 0.0;
    double full_number_fraction = 0.0;
    double st_energy_fraction = 0.0;
    double th_energy_fraction = 0.0;
    double full_energy_fraction = 0.0;
    std::array<double, 3> heat_delta_fraction = {0.0, 0.0, 0.0};
    double internal_transfer_number_fraction = 0.0;
    double internal_transfer_energy_fraction = 0.0;
    double internal_transfer_number_residual = 0.0;
    double internal_transfer_energy_residual = 0.0;
    double split_number_residual = 0.0;
    double split_energy_residual = 0.0;
    double full_number_residual = 0.0;
    double full_energy_residual = 0.0;
    double normalized_half_l1_vs_full = 0.0;
    double th_maxwellian_half_l1 = 0.0;
    double th_maxwellian_relative_entropy = 0.0;
};

struct MaxErrors {
    double normalized_half_l1_vs_full = 0.0;
    std::array<double, 3> heat_delta_fraction = {0.0, 0.0, 0.0};
    double internal_transfer_number_residual = 0.0;
    double internal_transfer_energy_residual = 0.0;
    double split_number_residual = 0.0;
    double split_energy_residual = 0.0;
    double full_number_residual = 0.0;
    double full_energy_residual = 0.0;
    double th_maxwellian_half_l1 = 0.0;
    double th_maxwellian_relative_entropy = 0.0;
};

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

double sum(const std::vector<double> &values) {
    long double result = 0.0L;
    for (double value : values) result += static_cast<long double>(value);
    require(std::isfinite(result), "non-finite population sum");
    require(result >= 0.0L && result <= std::numeric_limits<double>::max(),
            "population sum is outside double range");
    return static_cast<double>(result);
}

double weighted_energy(const std::vector<double> &population,
                       const std::vector<double> &centers) {
    require(population.size() == centers.size(), "energy vector extent mismatch");
    long double result = 0.0L;
    for (std::size_t i = 0; i < population.size(); ++i)
        result += static_cast<long double>(population[i]) *
                  static_cast<long double>(centers[i]);
    require(std::isfinite(result), "non-finite kinetic energy");
    require(result >= 0.0L && result <= std::numeric_limits<double>::max(),
            "kinetic energy is outside double range");
    return static_cast<double>(result);
}

std::vector<double> make_hybrid_edges(int cells) {
    require(cells >= 8, "cells must be >=8");
    const int split = cells / 4;
    require(split >= 2 && split < cells, "invalid hybrid-grid split");
    constexpr double low_keV = 0.001;
    constexpr double join_keV = 20.0;
    constexpr double high_keV = 6000.0;
    std::vector<double> edges(static_cast<std::size_t>(cells) + 1U);
    edges[0] = 0.0;
    for (int i = 1; i <= split; ++i) {
        const double fraction = static_cast<double>(i - 1) /
                                static_cast<double>(split - 1);
        edges[static_cast<std::size_t>(i)] =
            low_keV * std::pow(join_keV / low_keV, fraction) * kJoulesPerKeV;
    }
    for (int i = split + 1; i <= cells; ++i) {
        const double fraction = static_cast<double>(i - split) /
                                static_cast<double>(cells - split);
        edges[static_cast<std::size_t>(i)] =
            (join_keV + (high_keV - join_keV) * fraction) * kJoulesPerKeV;
    }
    require(edges[0] == 0.0 && edges[1] > 0.0,
            "hybrid grid zero endpoint failed");
    require(std::abs(edges[static_cast<std::size_t>(split)] / kJoulesPerKeV -
                     join_keV) < 1.0e-12,
            "hybrid grid join failed");
    for (int i = 1; i <= cells; ++i)
        require(std::isfinite(edges[static_cast<std::size_t>(i)]) &&
                    edges[static_cast<std::size_t>(i)] >
                        edges[static_cast<std::size_t>(i - 1)],
                "hybrid grid is not strictly increasing");
    require(std::abs(edges.back() / kJoulesPerKeV - high_keV) < 1.0e-10,
            "hybrid grid upper endpoint failed");
    return edges;
}

std::vector<double> centers_of(const std::vector<double> &edges) {
    std::vector<double> centers(edges.size() - 1U);
    for (std::size_t i = 0; i < centers.size(); ++i)
        centers[i] = 0.5 * (edges[i] + edges[i + 1]);
    return centers;
}

void initialize_alpha_birth(std::vector<double> &population,
                            const std::vector<double> &centers) {
    std::fill(population.begin(), population.end(), 0.0);
    const auto right_it = std::upper_bound(centers.begin(), centers.end(),
                                           kInitialAlphaEnergy);
    require(right_it != centers.begin() && right_it != centers.end(),
            "3 MeV birth is outside cell-center range");
    const std::size_t right = static_cast<std::size_t>(right_it - centers.begin());
    const double weight_right =
        (kInitialAlphaEnergy - centers[right - 1]) /
        (centers[right] - centers[right - 1]);
    population[right] = kPopulation * weight_right;
    population[right - 1] = kPopulation * (1.0 - weight_right);
    require(std::abs(sum(population) - kPopulation) < 1.0e-3,
            "birth number is not conserved");
    require(std::abs(weighted_energy(population, centers) -
                     kPopulation * kInitialAlphaEnergy) <
                1.0e-15 * kPopulation * kInitialAlphaEnergy,
            "birth energy is not conserved");
}

std::array<fusion_maxwellian_bath_v1, 3> make_baths() {
    return {{{1.0e20, kElectronMass, 1.0, kBathTemperature, kCoulombLog},
             {5.0e19, kDeuteronMass, 1.0, kBathTemperature, kCoulombLog},
             {5.0e19, kTritonMass, 1.0, kBathTemperature, kCoulombLog}}};
}

std::vector<double> make_diffusion(
    const std::vector<double> &edges,
    const std::array<fusion_maxwellian_bath_v1, 3> &baths) {
    const int cells = static_cast<int>(edges.size() - 1U);
    std::vector<double> diffusion(3U * static_cast<std::size_t>(cells - 1));
    for (int b = 0; b < 3; ++b) {
        for (int face = 0; face < cells - 1; ++face) {
            fusion_coulomb_energy_v1 coefficient{};
            const int status = fusion_c_coulomb_energy(
                edges[static_cast<std::size_t>(face + 1)], kAlphaMass, 2.0,
                &baths[static_cast<std::size_t>(b)], &coefficient);
            require(status == PB11_STATUS_OK,
                    "Coulomb diffusion coefficient failed");
            require(std::isfinite(coefficient.diffusion_J2_s) &&
                        coefficient.diffusion_J2_s >= 0.0,
                    "Coulomb diffusion coefficient is invalid");
            diffusion[static_cast<std::size_t>(b) *
                          static_cast<std::size_t>(cells - 1) +
                      static_cast<std::size_t>(face)] =
                coefficient.diffusion_J2_s;
        }
    }
    return diffusion;
}

std::array<long double, 2> ion_collision_rates() {
    const long double pi = std::acos(-1.0L);
    const long double charge = static_cast<long double>(kElementaryCharge);
    const long double coulomb =
        charge * charge /
        (4.0L * pi * static_cast<long double>(kEpsilon0));
    const long double za2 = 4.0L;
    const long double zi2 = 1.0L;
    const long double density = 5.0e19L;
    const long double ln_lambda = static_cast<long double>(kCoulombLog);
    const long double c = 4.0L * pi * density * za2 * zi2 *
                          coulomb * coulomb * ln_lambda;
    return {{c / (static_cast<long double>(kAlphaMass) *
                    static_cast<long double>(kDeuteronMass)),
             c / (static_cast<long double>(kAlphaMass) *
                    static_cast<long double>(kTritonMass))}};
}

std::vector<double> make_finite_temperature_lambda(
    const std::vector<double> &centers,
    const std::array<long double, 2> &nu, double &lambda_min,
    double &lambda_max, double &lambda_at_birth) {
    require(centers.size() >= 1U, "empty energy grid");
    const long double pi = std::acos(-1.0L);
    const long double sqrt_pi = std::sqrt(pi);
    const long double alpha_mass = static_cast<long double>(kAlphaMass);
    const long double alpha_temperature =
        static_cast<long double>(kBathTemperature);
    std::vector<double> lambda(centers.size());
    long double min_value = std::numeric_limits<long double>::max();
    long double max_value = 0.0L;
    std::size_t birth_index = 0U;
    double birth_distance = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < centers.size(); ++i) {
        const long double energy = static_cast<long double>(centers[i]);
        const long double v2 = 2.0L * energy / alpha_mass;
        long double total = 0.0L;
        for (std::size_t ion = 0; ion < 2U; ++ion) {
            const long double mass = ion == 0U
                                         ? static_cast<long double>(kDeuteronMass)
                                         : static_cast<long double>(kTritonMass);
            const long double ti = alpha_temperature;
            const long double vti2 = 2.0L * ti / mass;
            const long double vti = std::sqrt(vti2);
            const long double prefactor =
                4.0L * nu[ion] / (sqrt_pi * vti * vti * vti);
            const long double exponent = -v2 / vti2;
            const long double term = prefactor * std::exp(exponent);
            require(std::isfinite(term) && term >= 0.0L,
                    "invalid finite-temperature lambda term");
            total += term;
        }
        require(std::isfinite(total) && total >= 0.0L,
                "invalid finite-temperature lambda");
        require(total <= static_cast<long double>(
                             std::numeric_limits<double>::max()),
                "finite-temperature lambda exceeds double range");
        lambda[i] = static_cast<double>(total);
        require(std::isfinite(lambda[i]) && lambda[i] >= 0.0,
                "finite-temperature lambda conversion failed");
        min_value = std::min(min_value, total);
        max_value = std::max(max_value, total);
        const double distance =
            std::abs(centers[i] - kInitialAlphaEnergy);
        if (distance < birth_distance) {
            birth_distance = distance;
            birth_index = i;
        }
    }
    lambda_min = static_cast<double>(min_value);
    lambda_max = static_cast<double>(max_value);
    lambda_at_birth = lambda[birth_index];
    require(std::isfinite(lambda_min) && std::isfinite(lambda_max) &&
                std::isfinite(lambda_at_birth),
            "finite-temperature lambda diagnostics failed");
    return lambda;
}

void call_fp(int cells, double dt, const std::vector<double> &edges,
             const std::vector<double> &old,
             const double *temperatures,
             const std::vector<double> &diffusion,
             const std::vector<double> &birth,
             const std::vector<double> &escape,
             std::vector<double> &trial,
             std::vector<double> &heat,
             fusion_kinetic_ledger_v1 &ledger,
             const std::string &label) {
    const int status = fusion_c_energy_fp_trial(
        cells, 3, dt, edges.data(), old.data(), temperatures, diffusion.data(),
        birth.data(), escape.data(), 0.0, trial.data(), heat.data(), &ledger);
    if (status != PB11_STATUS_OK) {
        std::ostringstream message;
        message << label << " FP trial failed with status " << status;
        throw std::runtime_error(message.str());
    }
}

std::vector<long double> make_log_maxwellian_weights(
    const std::vector<double> &edges, const std::vector<double> &centers) {
    require(edges.size() == centers.size() + 1U,
            "Maxwellian grid extent mismatch");
    std::vector<long double> log_weight(centers.size());
    long double largest = -std::numeric_limits<long double>::infinity();
    for (std::size_t i = 0; i < centers.size(); ++i) {
        const long double width = static_cast<long double>(edges[i + 1]) -
                                  static_cast<long double>(edges[i]);
        const long double center = static_cast<long double>(centers[i]);
        require(width > 0.0L && center > 0.0L,
                "invalid Maxwellian cell geometry");
        const long double value = std::log(width) + 0.5L * std::log(center) -
                                  center / static_cast<long double>(kBathTemperature);
        require(std::isfinite(value), "invalid Maxwellian log weight");
        log_weight[i] = value;
        largest = std::max(largest, value);
    }
    long double scaled_sum = 0.0L;
    for (long double value : log_weight)
        scaled_sum += std::exp(value - largest);
    require(std::isfinite(scaled_sum) && scaled_sum > 0.0L,
            "invalid Maxwellian log normalization");
    const long double log_normalization = largest + std::log(scaled_sum);
    for (long double &value : log_weight) {
        value -= log_normalization;
        require(std::isfinite(value), "invalid normalized Maxwellian log weight");
    }
    return log_weight;
}

Metrics measure(const std::vector<double> &st,
                const std::vector<double> &th,
                const std::vector<double> &full,
                const std::vector<double> &centers,
                const std::array<double, 3> &split_heat,
                const std::array<double, 3> &full_heat,
                double escaped_number,
                double born_number,
                double escaped_energy,
                double born_energy,
                const std::vector<long double> &log_maxwellian) {
    require(st.size() == th.size() && st.size() == full.size() &&
                st.size() == centers.size() && st.size() == log_maxwellian.size(),
            "state extent mismatch");
    const double st_n = sum(st);
    const double th_n = sum(th);
    const double full_n = sum(full);
    const double st_e = weighted_energy(st, centers);
    const double th_e = weighted_energy(th, centers);
    const double full_e = weighted_energy(full, centers);
    std::vector<double> total(st.size());
    for (std::size_t i = 0; i < total.size(); ++i) {
        require(std::isfinite(st[i] + th[i]) && st[i] + th[i] >= 0.0,
                "invalid combined population");
        total[i] = st[i] + th[i];
    }
    const double total_n = sum(total);
    const double total_e = weighted_energy(total, centers);
    double l1 = 0.0;
    for (std::size_t i = 0; i < total.size(); ++i)
        l1 += std::abs(total[i] / kPopulation - full[i] / kPopulation);
    require(th_n > 0.0, "thermal component remained empty");
    long double th_l1 = 0.0L;
    long double th_entropy = 0.0L;
    for (std::size_t i = 0; i < th.size(); ++i) {
        const long double probability =
            static_cast<long double>(th[i]) / static_cast<long double>(th_n);
        require(std::isfinite(probability) && probability >= 0.0L,
                "invalid normalized TH probability");
        const long double reference = std::exp(log_maxwellian[i]);
        require(std::isfinite(reference) && reference > 0.0L,
                "invalid sampled Maxwellian probability");
        th_l1 += std::abs(probability - reference);
        if (probability > 0.0L)
            th_entropy += probability *
                          (std::log(probability) - log_maxwellian[i]);
    }
    const double number_scale = kPopulation;
    const double energy_scale = kPopulation * kInitialAlphaEnergy;
    Metrics result;
    result.st_number_fraction = st_n / number_scale;
    result.th_number_fraction = th_n / number_scale;
    result.full_number_fraction = full_n / number_scale;
    result.st_energy_fraction = st_e / energy_scale;
    result.th_energy_fraction = th_e / energy_scale;
    result.full_energy_fraction = full_e / energy_scale;
    for (int b = 0; b < 3; ++b)
        result.heat_delta_fraction[static_cast<std::size_t>(b)] =
            (split_heat[static_cast<std::size_t>(b)] -
             full_heat[static_cast<std::size_t>(b)]) /
            energy_scale;
    result.internal_transfer_number_fraction = escaped_number / number_scale;
    result.internal_transfer_energy_fraction = escaped_energy / energy_scale;
    result.internal_transfer_number_residual =
        (born_number - escaped_number) / number_scale;
    result.internal_transfer_energy_residual =
        (born_energy - escaped_energy) / energy_scale;
    result.split_number_residual = total_n / number_scale - 1.0;
    result.split_energy_residual =
        (total_e + split_heat[0] + split_heat[1] + split_heat[2]) /
            energy_scale -
        1.0;
    result.full_number_residual = full_n / number_scale - 1.0;
    result.full_energy_residual =
        (full_e + full_heat[0] + full_heat[1] + full_heat[2]) /
            energy_scale -
        1.0;
    result.normalized_half_l1_vs_full = 0.5 * l1;
    result.th_maxwellian_half_l1 = static_cast<double>(0.5L * th_l1);
    result.th_maxwellian_relative_entropy = static_cast<double>(th_entropy);
    require(std::isfinite(result.st_number_fraction) &&
                std::isfinite(result.th_number_fraction) &&
                std::isfinite(result.full_number_fraction) &&
                std::isfinite(result.st_energy_fraction) &&
                std::isfinite(result.th_energy_fraction) &&
                std::isfinite(result.full_energy_fraction) &&
                std::isfinite(result.internal_transfer_number_fraction) &&
                std::isfinite(result.internal_transfer_energy_fraction) &&
                std::isfinite(result.internal_transfer_number_residual) &&
                std::isfinite(result.internal_transfer_energy_residual) &&
                std::isfinite(result.split_number_residual) &&
                std::isfinite(result.split_energy_residual) &&
                std::isfinite(result.full_number_residual) &&
                std::isfinite(result.full_energy_residual) &&
                std::isfinite(result.normalized_half_l1_vs_full) &&
                std::isfinite(result.th_maxwellian_half_l1) &&
                std::isfinite(result.th_maxwellian_relative_entropy),
            "non-finite split diagnostic");
    require(result.th_maxwellian_half_l1 >= 0.0 &&
                result.th_maxwellian_relative_entropy > -1.0e-12,
            "invalid TH Maxwellian diagnostic");
    for (double value : result.heat_delta_fraction)
        require(std::isfinite(value), "non-finite bath heat difference");
    return result;
}

void update_maxima(MaxErrors &maxima, const Metrics &metrics) {
    maxima.normalized_half_l1_vs_full =
        std::max(maxima.normalized_half_l1_vs_full,
                 metrics.normalized_half_l1_vs_full);
    for (int b = 0; b < 3; ++b)
        maxima.heat_delta_fraction[static_cast<std::size_t>(b)] = std::max(
            maxima.heat_delta_fraction[static_cast<std::size_t>(b)],
            std::abs(metrics.heat_delta_fraction[static_cast<std::size_t>(b)]));
    maxima.internal_transfer_number_residual = std::max(
        maxima.internal_transfer_number_residual,
        std::abs(metrics.internal_transfer_number_residual));
    maxima.internal_transfer_energy_residual = std::max(
        maxima.internal_transfer_energy_residual,
        std::abs(metrics.internal_transfer_energy_residual));
    maxima.split_number_residual = std::max(
        maxima.split_number_residual, std::abs(metrics.split_number_residual));
    maxima.split_energy_residual = std::max(
        maxima.split_energy_residual, std::abs(metrics.split_energy_residual));
    maxima.full_number_residual = std::max(
        maxima.full_number_residual, std::abs(metrics.full_number_residual));
    maxima.full_energy_residual = std::max(
        maxima.full_energy_residual, std::abs(metrics.full_energy_residual));
    maxima.th_maxwellian_half_l1 =
        std::max(maxima.th_maxwellian_half_l1, metrics.th_maxwellian_half_l1);
    maxima.th_maxwellian_relative_entropy = std::max(
        maxima.th_maxwellian_relative_entropy,
        std::abs(metrics.th_maxwellian_relative_entropy));
}

void write_header(std::ofstream &csv) {
    csv << "time_s,st_number_fraction,th_number_fraction,full_number_fraction,"
           "st_energy_fraction,th_energy_fraction,full_energy_fraction,"
           "split_e_heat_fraction,split_D_heat_fraction,split_T_heat_fraction,"
           "full_e_heat_fraction,full_D_heat_fraction,full_T_heat_fraction,"
           "heat_delta_e_fraction,heat_delta_D_fraction,heat_delta_T_fraction,"
           "internal_transfer_number_fraction,internal_transfer_energy_fraction,"
           "internal_transfer_number_residual,internal_transfer_energy_residual,"
           "split_number_residual,split_energy_residual,full_number_residual,"
           "full_energy_residual,normalized_half_L1_vs_FULL,th_maxwellian_half_L1,"
           "th_maxwellian_relative_entropy\n";
}

int run(int cells, int steps, const std::string &csv_path,
        const std::string &log_path) {
    std::ofstream log(log_path);
    if (!log) throw std::runtime_error("cannot open log " + log_path);
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open CSV " + csv_path);
    log << std::setprecision(17);
    csv << std::setprecision(16);
    const auto started = std::chrono::steady_clock::now();
    try {
        require(cells >= 8 && steps >= 1, "cells must be >=8 and steps >=1");
        const double dt = kDuration / static_cast<double>(steps);
        const std::vector<double> edges = make_hybrid_edges(cells);
        const std::vector<double> centers = centers_of(edges);
        const auto baths = make_baths();
        const std::vector<double> diffusion = make_diffusion(edges, baths);
        const std::vector<long double> log_maxwellian =
            make_log_maxwellian_weights(edges, centers);
        const auto nu = ion_collision_rates();
        require(std::isfinite(nu[0]) && std::isfinite(nu[1]) &&
                    nu[0] > 0.0L && nu[1] > 0.0L,
                "invalid ion collision rate");
        double lambda_min = 0.0, lambda_max = 0.0, lambda_at_birth = 0.0;
        const std::vector<double> lambda = make_finite_temperature_lambda(
            centers, nu, lambda_min, lambda_max, lambda_at_birth);
        const double temperatures[3] = {kBathTemperature, kBathTemperature,
                                        kBathTemperature};
        std::vector<double> st(centers.size());
        std::vector<double> st_next(centers.size());
        std::vector<double> th(centers.size(), 0.0);
        std::vector<double> th_next(centers.size());
        std::vector<double> full(centers.size());
        std::vector<double> full_next(centers.size());
        std::vector<double> zero_birth(centers.size(), 0.0);
        std::vector<double> zero_escape(centers.size(), 0.0);
        std::vector<double> th_birth(centers.size(), 0.0);
        initialize_alpha_birth(st, centers);
        full = st;
        std::array<double, 3> split_heat = {0.0, 0.0, 0.0};
        std::array<double, 3> full_heat = {0.0, 0.0, 0.0};
        double escaped_number = 0.0;
        double born_number = 0.0;
        double escaped_energy = 0.0;
        double born_energy = 0.0;
        MaxErrors maxima;
        const int sample_every = std::max(1, steps / 200);

        write_header(csv);
        log << "cells=" << cells << " steps=" << steps << " dt_s=" << dt
            << " split=" << cells / 4
            << " e1_keV=" << edges[1] / kJoulesPerKeV
            << " join_keV="
            << edges[static_cast<std::size_t>(cells / 4)] / kJoulesPerKeV
            << " lambda_min_s_inv=" << lambda_min
            << " lambda_max_s_inv=" << lambda_max
            << " lambda_at_3MeV_s_inv=" << lambda_at_birth
            << " nu_D_v3_s=" << static_cast<double>(nu[0])
            << " nu_T_v3_s=" << static_cast<double>(nu[1]) << '\n';
        log << "lambda(E)=sum_i 4*nu_i/(sqrt(pi)*vti_i^3)*exp(-v(E)^2/vti_i^2); "
               "ST escape and TH birth are internal transfer only\n";

        for (int step = 1; step <= steps; ++step) {
#ifdef FUSION_STUDY_PUBLIC_API
            // Public mode only verifies total N/U, bath heat and FULL parity.
            // Its duplicated internal-transfer fields are bookkeeping, not
            // independent S-loss/T-birth subledger measurements.
            fusion_two_component_ledger_v1 combined{};
            double combined_heat[3]={};
            const int combined_status=fusion_c_two_component_trial(cells,3,dt,
                edges.data(),st.data(),th.data(),temperatures,diffusion.data(),
                zero_birth.data(),zero_birth.data(),zero_escape.data(),lambda.data(),
                st_next.data(),th_next.data(),combined_heat,&combined);
            require(combined_status==0,"public trial rejected step "+std::to_string(step)+
                " status "+std::to_string(combined_status));
            for(int b=0;b<3;b++) split_heat[b]+=combined_heat[b];
            fusion_kinetic_ledger_v1 st_ledger{},th_ledger{};
            st_ledger.escaped_number_m3=th_ledger.born_number_m3=combined.transferred_number_m3;
            st_ledger.escaped_energy_J_m3=th_ledger.born_energy_J_m3=combined.transferred_energy_J_m3;
#else
            fusion_kinetic_ledger_v1 st_ledger{};
            std::vector<double> st_heat(3, 0.0);
            call_fp(cells, dt, edges, st, temperatures, diffusion, zero_birth,
                    lambda, st_next, st_heat, st_ledger,
                    "ST/full step " + std::to_string(step));
            for (int b = 0; b < 3; ++b)
                split_heat[static_cast<std::size_t>(b)] +=
                    st_heat[static_cast<std::size_t>(b)];

            for (std::size_t i = 0; i < th_birth.size(); ++i) {
                const long double value =
                    static_cast<long double>(lambda[i]) *
                    static_cast<long double>(st_next[i]);
                require(std::isfinite(value) && value >= 0.0L &&
                            value <= static_cast<long double>(
                                         std::numeric_limits<double>::max()),
                        "invalid finite-temperature TH birth");
                th_birth[i] = static_cast<double>(value);
            }

            fusion_kinetic_ledger_v1 th_ledger{};
            std::vector<double> th_heat(3, 0.0);
            call_fp(cells, dt, edges, th, temperatures, diffusion, th_birth,
                    zero_escape, th_next, th_heat, th_ledger,
                    "TH/full step " + std::to_string(step));
            for (int b = 0; b < 3; ++b)
                split_heat[static_cast<std::size_t>(b)] +=
                    th_heat[static_cast<std::size_t>(b)];

#endif
            fusion_kinetic_ledger_v1 full_ledger{};
            std::vector<double> full_heat_step(3, 0.0);
            call_fp(cells, dt, edges, full, temperatures, diffusion, zero_birth,
                    zero_escape, full_next, full_heat_step, full_ledger,
                    "FULL step " + std::to_string(step));
            for (int b = 0; b < 3; ++b)
                full_heat[static_cast<std::size_t>(b)] +=
                    full_heat_step[static_cast<std::size_t>(b)];

            require(std::isfinite(st_ledger.escaped_number_m3) &&
                        std::isfinite(st_ledger.escaped_energy_J_m3) &&
                        std::isfinite(th_ledger.born_number_m3) &&
                        std::isfinite(th_ledger.born_energy_J_m3) &&
                        st_ledger.escaped_number_m3 >= 0.0 &&
                        st_ledger.escaped_energy_J_m3 >= 0.0 &&
                        th_ledger.born_number_m3 >= 0.0 &&
                        th_ledger.born_energy_J_m3 >= 0.0,
                    "invalid internal transfer ledger");
            escaped_number += st_ledger.escaped_number_m3;
            born_number += th_ledger.born_number_m3;
            escaped_energy += st_ledger.escaped_energy_J_m3;
            born_energy += th_ledger.born_energy_J_m3;
            require(std::isfinite(escaped_number) && std::isfinite(born_number) &&
                        std::isfinite(escaped_energy) && std::isfinite(born_energy),
                    "non-finite cumulative internal transfer");

            st.swap(st_next);
            th.swap(th_next);
            full.swap(full_next);
            const Metrics metrics = measure(
                st, th, full, centers, split_heat, full_heat, escaped_number,
                born_number, escaped_energy, born_energy, log_maxwellian);
            update_maxima(maxima, metrics);
            require(maxima.split_number_residual <= 1.0e-8 &&
                        maxima.split_energy_residual <= 1.0e-8 &&
                        maxima.full_number_residual <= 1.0e-8 &&
                        maxima.full_energy_residual <= 1.0e-8,
                    "cumulative N/U conservation check failed");

            require(maxima.normalized_half_l1_vs_full <= 1.e-11,
                    "full-distribution parity check failed");
            for (double delta : maxima.heat_delta_fraction)
                require(delta <= 1.e-11, "per-bath heat parity check failed");
            if (step % sample_every == 0 || step == steps) {
                csv << step * dt << ',' << metrics.st_number_fraction << ','
                    << metrics.th_number_fraction << ','
                    << metrics.full_number_fraction << ','
                    << metrics.st_energy_fraction << ','
                    << metrics.th_energy_fraction << ','
                    << metrics.full_energy_fraction << ','
                    << split_heat[0] / (kPopulation * kInitialAlphaEnergy) << ','
                    << split_heat[1] / (kPopulation * kInitialAlphaEnergy) << ','
                    << split_heat[2] / (kPopulation * kInitialAlphaEnergy) << ','
                    << full_heat[0] / (kPopulation * kInitialAlphaEnergy) << ','
                    << full_heat[1] / (kPopulation * kInitialAlphaEnergy) << ','
                    << full_heat[2] / (kPopulation * kInitialAlphaEnergy) << ','
                    << metrics.heat_delta_fraction[0] << ','
                    << metrics.heat_delta_fraction[1] << ','
                    << metrics.heat_delta_fraction[2] << ','
                    << metrics.internal_transfer_number_fraction << ','
                    << metrics.internal_transfer_energy_fraction << ','
                    << metrics.internal_transfer_number_residual << ','
                    << metrics.internal_transfer_energy_residual << ','
                    << metrics.split_number_residual << ','
                    << metrics.split_energy_residual << ','
                    << metrics.full_number_residual << ','
                    << metrics.full_energy_residual << ','
                    << metrics.normalized_half_l1_vs_full << ','
                    << metrics.th_maxwellian_half_l1 << ','
                    << metrics.th_maxwellian_relative_entropy << '\n';
            }
        }

        const Metrics final = measure(st, th, full, centers, split_heat,
                                      full_heat, escaped_number, born_number,
                                      escaped_energy, born_energy, log_maxwellian);
        const auto finished = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(finished - started).count();
        log << "final st_N=" << final.st_number_fraction
            << " st_U=" << final.st_energy_fraction
            << " th_N=" << final.th_number_fraction
            << " th_U=" << final.th_energy_fraction
            << " full_N=" << final.full_number_fraction
            << " full_U=" << final.full_energy_fraction
            << " normalized_half_L1_vs_FULL=" << final.normalized_half_l1_vs_full
            << " final_heat_delta_e=" << final.heat_delta_fraction[0]
            << " final_heat_delta_D=" << final.heat_delta_fraction[1]
            << " final_heat_delta_T=" << final.heat_delta_fraction[2]
            << " internal_transfer_N_fraction="
            << final.internal_transfer_number_fraction
            << " internal_transfer_E_fraction="
            << final.internal_transfer_energy_fraction
            << " internal_transfer_N_residual="
            << final.internal_transfer_number_residual
            << " internal_transfer_E_residual="
            << final.internal_transfer_energy_residual
            << " split_N_residual=" << final.split_number_residual
            << " split_E_residual=" << final.split_energy_residual
            << " full_N_residual=" << final.full_number_residual
            << " full_E_residual=" << final.full_energy_residual
            << " th_maxwellian_half_L1=" << final.th_maxwellian_half_l1
            << " th_maxwellian_relative_entropy="
            << final.th_maxwellian_relative_entropy << '\n';
        log << "max_over_all_steps normalized_half_L1_vs_FULL="
            << maxima.normalized_half_l1_vs_full
            << " abs_heat_delta_e=" << maxima.heat_delta_fraction[0]
            << " abs_heat_delta_D=" << maxima.heat_delta_fraction[1]
            << " abs_heat_delta_T=" << maxima.heat_delta_fraction[2]
            << " abs_internal_transfer_N_residual="
            << maxima.internal_transfer_number_residual
            << " abs_internal_transfer_E_residual="
            << maxima.internal_transfer_energy_residual
            << " abs_split_N_residual=" << maxima.split_number_residual
            << " abs_split_E_residual=" << maxima.split_energy_residual
            << " abs_full_N_residual=" << maxima.full_number_residual
            << " abs_full_E_residual=" << maxima.full_energy_residual
            << " max_th_maxwellian_half_L1=" << maxima.th_maxwellian_half_l1
            << " max_abs_th_maxwellian_relative_entropy="
            << maxima.th_maxwellian_relative_entropy
            << " elapsed_s=" << elapsed << '\n';
        log << "status=OK; finite-temperature ST/TH transfer is internal bookkeeping only; "
               "this is the supplied numerical prototype, not a physical escape or fluid ash model\n";
        return 0;
    } catch (const std::exception &error) {
        log << "status=ERROR message=" << error.what() << '\n';
        return 1;
    }
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc < 3 || argc > 5) {
            std::cerr << "usage: study_finite_temperature_split cells steps [csv] [log]\n";
            return 2;
        }
        const int cells = std::stoi(argv[1]);
        const int steps = std::stoi(argv[2]);
        const std::string base = "/tmp/finite-temperature-split-n" +
                                 std::to_string(cells) + "-s" +
                                 std::to_string(steps);
        const std::string csv = argc >= 4 ? argv[3] : base + ".csv";
        const std::string log = argc >= 5 ? argv[4] : base + ".log";
        return run(cells, steps, csv, log);
    } catch (const std::exception &error) {
        std::cerr << "status=ERROR message=" << error.what() << '\n';
        return 1;
    }
}
