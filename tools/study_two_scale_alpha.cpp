// Numerical two-scale trace-alpha prototype.  This is a transcription of
// the explicitly supplied recipe, not a production API or a claim of a full
// Peigney-spectrum/thermal-ash reproduction.
#include "fusion_coulomb.h"
#include "fusion_kinetics.h"

#include <algorithm>
#include <array>
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
// TH is linear in its cell-integrated population and source. This numerical
// scale keeps the first tiny source energy representable by the FP kernel;
// populations and heats are divided by it immediately after each TH call.
constexpr double kThermalPopulationScale = 1.0e20;

struct Metrics {
    double st_number_fraction = 0.0;
    double st_energy_fraction = 0.0;
    double th_number_fraction = 0.0;
    double th_energy_fraction = 0.0;
    double full_number_fraction = 0.0;
    double full_energy_fraction = 0.0;
    double two_scale_number_residual = 0.0;
    double two_scale_energy_residual = 0.0;
    double full_number_residual = 0.0;
    double full_energy_residual = 0.0;
    double normalized_half_l1_vs_full = 0.0;
    double maxwellian_overlap = 0.0;
};

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

double sum(const std::vector<double> &values) {
    long double result = 0.0L;
    for (double value : values) result += value;
    require(std::isfinite(result), "non-finite population sum");
    return static_cast<double>(result);
}

double weighted_energy(const std::vector<double> &population,
                       const std::vector<double> &centers) {
    require(population.size() == centers.size(), "energy vector extent mismatch");
    long double result = 0.0L;
    for (std::size_t i = 0; i < population.size(); ++i)
        result += static_cast<long double>(population[i]) * centers[i];
    require(std::isfinite(result), "non-finite kinetic energy");
    return static_cast<double>(result);
}

std::vector<double> make_hybrid_edges(int cells) {
    require(cells >= 8, "cells must leave at least two logarithmic intervals");
    const int split = cells / 4;
    require(split >= 2 && split < cells, "invalid hybrid-grid split");
    const double low_keV = 0.001;
    const double join_keV = 20.0;
    const double high_keV = 6000.0;
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
    require(edges[0] == 0.0 && edges[1] > 0.0, "hybrid grid zero endpoint failed");
    require(std::abs(edges[split] / kJoulesPerKeV - join_keV) < 1.0e-12,
            "hybrid grid logarithmic/linear join failed");
    for (int i = 1; i <= cells; ++i)
        require(std::isfinite(edges[static_cast<std::size_t>(i)]) &&
                    edges[static_cast<std::size_t>(i)] > edges[static_cast<std::size_t>(i - 1)],
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

std::vector<double> make_diffusion(const std::vector<double> &edges,
                                   const std::array<fusion_maxwellian_bath_v1, 3> &baths) {
    const int cells = static_cast<int>(edges.size() - 1U);
    std::vector<double> diffusion(3U * static_cast<std::size_t>(cells - 1));
    for (int b = 0; b < 3; ++b) {
        for (int face = 0; face < cells - 1; ++face) {
            fusion_coulomb_energy_v1 coefficient{};
            const int status = fusion_c_coulomb_energy(
                edges[static_cast<std::size_t>(face + 1)], kAlphaMass, 2.0,
                &baths[static_cast<std::size_t>(b)], &coefficient);
            require(status == PB11_STATUS_OK, "Coulomb diffusion coefficient failed");
            require(std::isfinite(coefficient.diffusion_J2_s) &&
                        coefficient.diffusion_J2_s >= 0.0,
                    "Coulomb diffusion coefficient is invalid");
            diffusion[static_cast<std::size_t>(b) * static_cast<std::size_t>(cells - 1) +
                      static_cast<std::size_t>(face)] = coefficient.diffusion_J2_s;
        }
    }
    return diffusion;
}

std::vector<double> make_electron_diffusion(const std::vector<double> &edges,
                                             const std::array<fusion_maxwellian_bath_v1, 3> &baths) {
    const int cells = static_cast<int>(edges.size() - 1U);
    std::vector<double> diffusion(static_cast<std::size_t>(cells - 1));
    for (int face = 0; face < cells - 1; ++face) {
        fusion_coulomb_energy_v1 coefficient{};
        const int status = fusion_c_coulomb_energy(
            edges[static_cast<std::size_t>(face + 1)], kAlphaMass, 2.0,
            &baths[0], &coefficient);
        require(status == PB11_STATUS_OK, "electron diffusion coefficient failed");
        require(std::isfinite(coefficient.diffusion_J2_s) &&
                    coefficient.diffusion_J2_s >= 0.0,
                "electron diffusion coefficient is invalid");
        diffusion[static_cast<std::size_t>(face)] = coefficient.diffusion_J2_s;
    }
    return diffusion;
}

std::array<long double, 2> ion_collision_rates() {
    const long double pi = std::acos(-1.0L);
    const long double charge = static_cast<long double>(kElementaryCharge);
    const long double coulomb = charge * charge /
                                (4.0L * pi * static_cast<long double>(kEpsilon0));
    const long double za2 = 4.0L;
    const long double density = 5.0e19L;
    const long double ln_lambda = static_cast<long double>(kCoulombLog);
    const long double common = 4.0L * pi * density * za2 * coulomb * coulomb * ln_lambda /
                               static_cast<long double>(kAlphaMass);
    return {{common / static_cast<long double>(kDeuteronMass),
             common / static_cast<long double>(kTritonMass)}};
}

std::vector<double> make_velocity_rates(const std::vector<double> &edges,
                                        long double nu) {
    const int cells = static_cast<int>(edges.size() - 1U);
    std::vector<double> rates(static_cast<std::size_t>(cells));
    for (int i = 0; i < cells; ++i) {
        const long double vlo = std::sqrt(2.0L * static_cast<long double>(edges[static_cast<std::size_t>(i)]) /
                                          static_cast<long double>(kAlphaMass));
        const long double vhi = std::sqrt(2.0L * static_cast<long double>(edges[static_cast<std::size_t>(i + 1)]) /
                                          static_cast<long double>(kAlphaMass));
        const long double volume = (vhi * vhi * vhi - vlo * vlo * vlo) / 3.0L;
        require(std::isfinite(volume) && volume > 0.0L, "invalid velocity-cell volume");
        const long double rate = nu / volume;
        require(std::isfinite(rate) && rate >= 0.0L &&
                    std::abs(rate) <= std::numeric_limits<double>::max(),
                "invalid ion velocity rate");
        rates[static_cast<std::size_t>(i)] = static_cast<double>(rate);
    }
    return rates;
}

std::vector<double> make_source_shape(const std::vector<double> &edges,
                                      const std::vector<double> &centers,
                                      double temperature,
                                      double &mean_energy,
                                      double &mean_relative_error) {
    require(temperature > 0.0 && std::isfinite(temperature),
            "invalid source temperature");
    std::vector<double> shape(centers.size());
    long double normalization = 0.0L;
    long double energy = 0.0L;
    for (std::size_t i = 0; i < shape.size(); ++i) {
        const long double width = static_cast<long double>(edges[i + 1] - edges[i]);
        const long double center = static_cast<long double>(centers[i]);
        const long double value = width * std::sqrt(center) *
                                  std::exp(-center / static_cast<long double>(temperature));
        require(std::isfinite(value) && value >= 0.0L, "invalid Maxwellian source bin");
        shape[i] = static_cast<double>(value);
        normalization += value;
        energy += value * center;
    }
    require(std::isfinite(normalization) && normalization > 0.0L,
            "empty Maxwellian source normalization");
    const long double mean = energy / normalization;
    mean_energy = static_cast<double>(mean);
    mean_relative_error = static_cast<double>(
        (mean - 1.5L * static_cast<long double>(temperature)) /
        (1.5L * static_cast<long double>(temperature)));
    for (double &value : shape) value /= static_cast<double>(normalization);
    require(std::isfinite(mean_energy) && std::isfinite(mean_relative_error),
            "invalid Maxwellian source mean");
    return shape;
}

void call_fp(int cells, int baths, double dt, const std::vector<double> &edges,
             const std::vector<double> &old, const double *temperatures,
             const std::vector<double> &diffusion, const std::vector<double> &birth,
             const std::vector<double> &escape, std::vector<double> &trial,
             std::vector<double> &heat, fusion_kinetic_ledger_v1 &ledger,
             const std::string &label) {
    const int status = fusion_c_energy_fp_trial(
        cells, baths, dt, edges.data(), old.data(), temperatures,
        diffusion.empty() ? nullptr : diffusion.data(), birth.data(), escape.data(),
        0.0, trial.data(), heat.empty() ? nullptr : heat.data(), &ledger);
    if (status != PB11_STATUS_OK) {
        std::ostringstream message;
        message << label << " FP trial failed with status " << status;
        throw std::runtime_error(message.str());
    }
}

Metrics measure(const std::vector<double> &st, const std::vector<double> &th,
                const std::vector<double> &full, const std::vector<double> &centers,
                const std::vector<double> &maxwellian, const std::array<double, 3> &two_heat,
                const std::array<double, 3> &full_heat, const std::array<double, 2> &cold_heat) {
    const double st_n = sum(st), th_n = sum(th), full_n = sum(full);
    const double st_e = weighted_energy(st, centers);
    const double th_e = weighted_energy(th, centers);
    const double full_e = weighted_energy(full, centers);
    const double energy_scale = kPopulation * kInitialAlphaEnergy;
    std::vector<double> total(st.size());
    for (std::size_t i = 0; i < total.size(); ++i) total[i] = st[i] + th[i];
    const double total_n = sum(total);
    const double total_e = weighted_energy(total, centers);
    double l1 = 0.0;
    double overlap = 0.0;
    for (std::size_t i = 0; i < total.size(); ++i) {
        const double two_density = total[i] / kPopulation;
        const double full_density = full[i] / kPopulation;
        l1 += std::abs(two_density - full_density);
        overlap += std::min(two_density, maxwellian[i]);
    }
    Metrics result;
    result.st_number_fraction = st_n / kPopulation;
    result.st_energy_fraction = st_e / energy_scale;
    result.th_number_fraction = th_n / kPopulation;
    result.th_energy_fraction = th_e / energy_scale;
    result.full_number_fraction = full_n / kPopulation;
    result.full_energy_fraction = full_e / energy_scale;
    result.two_scale_number_residual = total_n / kPopulation - 1.0;
    result.two_scale_energy_residual =
        (total_e + two_heat[0] + two_heat[1] + two_heat[2]) / energy_scale - 1.0;
    result.full_number_residual = full_n / kPopulation - 1.0;
    result.full_energy_residual =
        (full_e + full_heat[0] + full_heat[1] + full_heat[2]) / energy_scale - 1.0;
    result.normalized_half_l1_vs_full = 0.5 * l1;
    result.maxwellian_overlap = overlap;
    require(std::isfinite(result.two_scale_number_residual) &&
                std::isfinite(result.two_scale_energy_residual) &&
                std::isfinite(result.full_number_residual) &&
                std::isfinite(result.full_energy_residual) &&
                std::isfinite(result.normalized_half_l1_vs_full) &&
                std::isfinite(result.maxwellian_overlap),
            "non-finite output metric");
    (void)cold_heat;
    return result;
}

void write_header(std::ofstream &csv) {
    csv << "time_s,st_number_fraction,st_energy_fraction,th_number_fraction,"
           "th_energy_fraction,full_number_fraction,full_energy_fraction,"
           "two_scale_electron_heat_fraction,two_scale_D_heat_fraction,"
           "two_scale_T_heat_fraction,full_electron_heat_fraction,"
           "full_D_heat_fraction,full_T_heat_fraction,ion_cold_D_heat_fraction,"
           "ion_cold_T_heat_fraction,two_scale_number_residual,"
           "two_scale_energy_residual,full_number_residual,full_energy_residual,"
           "normalized_half_L1_vs_FULL,maxwellian_overlap,source_mean_relative_error_D,"
           "source_mean_relative_error_T\n";
}

int run(int cells, int steps, const std::string &csv_path,
        const std::string &log_path) {
    std::ofstream log(log_path);
    if (!log) throw std::runtime_error("cannot open log " + log_path);
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("cannot open CSV " + csv_path);
    log << std::setprecision(17);
    csv << std::setprecision(16);
    try {
        require(cells >= 8 && steps >= 1, "cells must be >=8 and steps >=1");
        const double dt = kDuration / static_cast<double>(steps);
        const std::vector<double> edges = make_hybrid_edges(cells);
        const std::vector<double> centers = centers_of(edges);
        const auto baths = make_baths();
        const auto diffusion_full = make_diffusion(edges, baths);
        const auto diffusion_electron = make_electron_diffusion(edges, baths);
        const auto nu_i = ion_collision_rates();
        const long double nu_total = nu_i[0] + nu_i[1];
        require(std::isfinite(nu_total) && nu_total > 0.0L, "invalid total ion collision rate");
        const std::array<double, 2> ion_weights = {
            static_cast<double>(nu_i[0] / nu_total),
            static_cast<double>(nu_i[1] / nu_total)};
        const std::vector<double> rate = make_velocity_rates(edges, nu_total);
        const double source_temperature_D =
            (kAlphaMass / kDeuteronMass) * kBathTemperature;
        const double source_temperature_T =
            (kAlphaMass / kTritonMass) * kBathTemperature;
        double source_mean_D = 0.0, source_error_D = 0.0;
        double source_mean_T = 0.0, source_error_T = 0.0;
        const std::vector<double> source_D = make_source_shape(
            edges, centers, source_temperature_D, source_mean_D, source_error_D);
        const std::vector<double> source_T = make_source_shape(
            edges, centers, source_temperature_T, source_mean_T, source_error_T);
        std::vector<double> maxwellian(centers.size());
        long double maxwellian_norm = 0.0L;
        for (std::size_t i = 0; i < centers.size(); ++i) {
            const long double width = static_cast<long double>(edges[i + 1] - edges[i]);
            const long double value = width * std::sqrt(static_cast<long double>(centers[i])) *
                                      std::exp(-static_cast<long double>(centers[i]) /
                                               static_cast<long double>(kBathTemperature));
            maxwellian[i] = static_cast<double>(value);
            maxwellian_norm += value;
        }
        require(maxwellian_norm > 0.0L && std::isfinite(maxwellian_norm),
                "invalid reference Maxwellian");
        for (double &value : maxwellian)
            value /= static_cast<double>(maxwellian_norm);

        std::vector<double> st(centers.size()), st_fp(centers.size()), st_adv(centers.size());
        std::vector<double> th(centers.size()), th_next(centers.size());
        std::vector<double> th_scaled(centers.size()), th_next_scaled(centers.size());
        std::vector<double> source_birth_scaled(centers.size());
        std::vector<double> full(centers.size()), full_next(centers.size());
        std::vector<double> zero_birth(centers.size(), 0.0);
        std::vector<double> zero_escape(centers.size(), 0.0);
        std::vector<double> source_birth(centers.size(), 0.0);
        initialize_alpha_birth(st, centers);
        th.assign(th.size(), 0.0);
        full = st;
        std::array<double, 3> two_heat = {0.0, 0.0, 0.0};
        std::array<double, 3> full_heat = {0.0, 0.0, 0.0};
        std::array<double, 2> cold_heat = {0.0, 0.0};
        const double temperatures[3] = {kBathTemperature, kBathTemperature,
                                        kBathTemperature};
        const double electron_temperature[1] = {kBathTemperature};
        const int sample_every = std::max(1, steps / 200);

        write_header(csv);
        log << "cells=" << cells << " steps=" << steps << " dt_s=" << dt
            << " split=" << cells / 4 << " e1_keV=" << edges[1] / kJoulesPerKeV
            << " join_keV=" << edges[static_cast<std::size_t>(cells / 4)] / kJoulesPerKeV
            << " nu_D_v3_s=" << static_cast<double>(nu_i[0])
            << " nu_T_v3_s=" << static_cast<double>(nu_i[1])
            << " weight_D=" << ion_weights[0] << " weight_T=" << ion_weights[1]
            << " source_mean_D_keV=" << source_mean_D / kJoulesPerKeV
            << " source_mean_T_keV=" << source_mean_T / kJoulesPerKeV
            << " source_error_D=" << source_error_D
            << " source_error_T=" << source_error_T
            << " thermal_population_scale=" << kThermalPopulationScale << '\n';

        for (int step = 1; step <= steps; ++step) {
            fusion_kinetic_ledger_v1 st_ledger{};
            std::vector<double> st_heat(1, 0.0);
            call_fp(cells, 1, dt, edges, st, electron_temperature,
                    diffusion_electron, zero_birth, zero_escape, st_fp, st_heat,
                    st_ledger, "ST/electron step " + std::to_string(step));
            two_heat[0] += st_heat[0];

            const double st_before_advection_energy = weighted_energy(st_fp, centers);
            std::fill(st_adv.begin(), st_adv.end(), 0.0);
            for (int j = cells - 1; j >= 0; --j) {
                const double incoming =
                    j + 1 < cells ? dt * rate[static_cast<std::size_t>(j + 1)] *
                                         st_adv[static_cast<std::size_t>(j + 1)]
                                   : 0.0;
                const double denominator =
                    1.0 + dt * rate[static_cast<std::size_t>(j)];
                st_adv[static_cast<std::size_t>(j)] =
                    (st_fp[static_cast<std::size_t>(j)] + incoming) / denominator;
                require(std::isfinite(st_adv[static_cast<std::size_t>(j)]) &&
                            st_adv[static_cast<std::size_t>(j)] >= 0.0,
                        "negative/non-finite ST advection result");
            }
            const double st_after_advection_energy = weighted_energy(st_adv, centers);
            const double dUcold = st_before_advection_energy - st_after_advection_energy;
            require(std::isfinite(dUcold), "non-finite cold-ion energy decrement");
            const double R = dt * rate[0] * st_adv[0];
            require(std::isfinite(R) && R >= 0.0, "invalid low-energy injection count");
            const double injected_D = R * ion_weights[0];
            const double injected_T = R * ion_weights[1];
            for (std::size_t i = 0; i < source_birth.size(); ++i)
                source_birth[i] =
                    (injected_D * source_D[i] + injected_T * source_T[i]) / dt;
            const double Uinject_D = injected_D * source_mean_D;
            const double Uinject_T = injected_T * source_mean_T;
            cold_heat[0] += ion_weights[0] * dUcold - Uinject_D;
            cold_heat[1] += ion_weights[1] * dUcold - Uinject_T;
            two_heat[1] += ion_weights[0] * dUcold - Uinject_D;
            two_heat[2] += ion_weights[1] * dUcold - Uinject_T;

            for (std::size_t i = 0; i < th.size(); ++i) {
                th_scaled[i] = kThermalPopulationScale * th[i];
                source_birth_scaled[i] = kThermalPopulationScale * source_birth[i];
            }
            fusion_kinetic_ledger_v1 th_ledger{};
            std::vector<double> th_heat_scaled(3, 0.0);
            call_fp(cells, 3, dt, edges, th_scaled, temperatures, diffusion_full,
                    source_birth_scaled, zero_escape, th_next_scaled,
                    th_heat_scaled, th_ledger,
                    "TH/full step " + std::to_string(step));
            for (std::size_t i = 0; i < th.size(); ++i) {
                th_next[i] = th_next_scaled[i] / kThermalPopulationScale;
                require(std::isfinite(th_next[i]) && th_next[i] >= 0.0,
                        "invalid unscaled TH population");
            }
            for (int b = 0; b < 3; ++b)
                two_heat[static_cast<std::size_t>(b)] +=
                    th_heat_scaled[static_cast<std::size_t>(b)] /
                    kThermalPopulationScale;

            fusion_kinetic_ledger_v1 full_ledger{};
            std::vector<double> full_heat_step(3, 0.0);
            call_fp(cells, 3, dt, edges, full, temperatures, diffusion_full,
                    zero_birth, zero_escape, full_next, full_heat_step,
                    full_ledger, "FULL step " + std::to_string(step));
            for (int b = 0; b < 3; ++b) full_heat[static_cast<std::size_t>(b)] += full_heat_step[static_cast<std::size_t>(b)];

            st.swap(st_adv);
            th.swap(th_next);
            full.swap(full_next);

            if (step % sample_every == 0 || step == steps) {
                const Metrics metrics = measure(st, th, full, centers, maxwellian,
                                                two_heat, full_heat, cold_heat);
                csv << step * dt << ',' << metrics.st_number_fraction << ','
                    << metrics.st_energy_fraction << ',' << metrics.th_number_fraction << ','
                    << metrics.th_energy_fraction << ',' << metrics.full_number_fraction << ','
                    << metrics.full_energy_fraction << ',' << two_heat[0] / (kPopulation * kInitialAlphaEnergy)
                    << ',' << two_heat[1] / (kPopulation * kInitialAlphaEnergy) << ','
                    << two_heat[2] / (kPopulation * kInitialAlphaEnergy) << ','
                    << full_heat[0] / (kPopulation * kInitialAlphaEnergy) << ','
                    << full_heat[1] / (kPopulation * kInitialAlphaEnergy) << ','
                    << full_heat[2] / (kPopulation * kInitialAlphaEnergy) << ','
                    << cold_heat[0] / (kPopulation * kInitialAlphaEnergy) << ','
                    << cold_heat[1] / (kPopulation * kInitialAlphaEnergy) << ','
                    << metrics.two_scale_number_residual << ','
                    << metrics.two_scale_energy_residual << ','
                    << metrics.full_number_residual << ',' << metrics.full_energy_residual << ','
                    << metrics.normalized_half_l1_vs_full << ',' << metrics.maxwellian_overlap << ','
                    << source_error_D << ',' << source_error_T << '\n';
            }
        }
        const Metrics final = measure(st, th, full, centers, maxwellian,
                                      two_heat, full_heat, cold_heat);
        log << "final st_N=" << final.st_number_fraction
            << " st_U=" << final.st_energy_fraction
            << " th_N=" << final.th_number_fraction
            << " th_U=" << final.th_energy_fraction
            << " full_N=" << final.full_number_fraction
            << " full_U=" << final.full_energy_fraction
            << " two_scale_N_residual=" << final.two_scale_number_residual
            << " two_scale_E_residual=" << final.two_scale_energy_residual
            << " full_N_residual=" << final.full_number_residual
            << " full_E_residual=" << final.full_energy_residual
            << " normalized_half_L1_vs_FULL=" << final.normalized_half_l1_vs_full
            << " maxwellian_overlap=" << final.maxwellian_overlap
            << " cumulative_two_heat_e=" << two_heat[0]
            << " cumulative_two_heat_D=" << two_heat[1]
            << " cumulative_two_heat_T=" << two_heat[2]
            << " cumulative_full_heat_e=" << full_heat[0]
            << " cumulative_full_heat_D=" << full_heat[1]
            << " cumulative_full_heat_T=" << full_heat[2]
            << " cumulative_cold_heat_D=" << cold_heat[0]
            << " cumulative_cold_heat_T=" << cold_heat[1] << '\n';
        log << "status=OK; this is the supplied two-scale numerical prototype only\n";
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
            std::cerr << "usage: study_two_scale_alpha cells steps [csv] [log]\n";
            return 2;
        }
        const int cells = std::stoi(argv[1]);
        const int steps = std::stoi(argv[2]);
        const std::string base = "/tmp/two-scale-n" + std::to_string(cells) +
                                 "-s" + std::to_string(steps);
        const std::string csv = argc >= 4 ? argv[3] : base + ".csv";
        const std::string log = argc >= 5 ? argv[4] : base + ".log";
        return run(cells, steps, csv, log);
    } catch (const std::exception &error) {
        std::cerr << "status=ERROR message=" << error.what() << '\n';
        return 1;
    }
}
