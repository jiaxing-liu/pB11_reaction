#include "fusion_pb_birth.h"
#include "fusion_nuclear_data.h"
#include "fusion_products.h"
#include "fusion_alpha_spectrum.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using R = long double;
// Keep MeV as the library's canonical literal.  Deriving it by multiplying
// the keV literal can put the 1-keV cutoff one ULP below the accepted bound.
constexpr double mev = 1.602176634e-13;
constexpr double kev = .001 * mev;

void require(bool condition, const char *what) {
    if (!condition) throw std::runtime_error(what);
}

bool close_value(R a, R b, R relative = 3e-11L, R absolute = 1e-30L) {
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    return std::abs(a - b) <= absolute + relative * std::max(std::abs(a), std::abs(b));
}

std::vector<double> uniform_edges(int cells, double lo, double hi) {
    std::vector<double> result(static_cast<std::size_t>(cells) + 1);
    for (int i = 0; i <= cells; ++i)
        result[static_cast<std::size_t>(i)] = lo + (hi - lo) * i / cells;
    return result;
}

struct AlphaResult {
    std::vector<double> birth;
    fusion_pb_birth_v1 ledger{};
};

AlphaResult evaluate_alpha(double available, double q, const std::vector<double> &edges) {
    AlphaResult result;
    result.birth.assign(edges.size() - 1, 0.0);
    require(fusion_c_pb_alpha0_grid(available, q, static_cast<int>(result.birth.size()),
                                    edges.data(), result.birth.data(), &result.ledger) ==
                PB11_STATUS_OK,
            "alpha0 grid status");
    return result;
}

void check_common_inventory(const std::vector<double> &edges, double available,
                            const fusion_pb_birth_v1 &ledger,
                            const char *label, bool require_spill = false) {
    require(std::isfinite(ledger.mapped_number) && std::isfinite(ledger.mapped_energy_J), label);
    require(std::isfinite(ledger.below_number) && std::isfinite(ledger.below_energy_J) &&
                std::isfinite(ledger.above_number) && std::isfinite(ledger.above_energy_J),
            label);
    require(ledger.below_number >= 0.0 && ledger.above_number >= 0.0 &&
                ledger.mapped_number >= 0.0,
            label);
    require(close_value(static_cast<R>(ledger.mapped_number) + ledger.below_number +
                            ledger.above_number,
                        3.0L, 3e-11L, 2e-13L),
            label);
    require(close_value(static_cast<R>(ledger.mapped_energy_J) + ledger.below_energy_J +
                            ledger.above_energy_J,
                        available, 3e-11L, 2e-30L),
            label);
    require(std::abs(ledger.number_residual) < 4e-12 &&
                std::abs(ledger.energy_residual_J) < 4e-12 * available,
            label);
    if (require_spill)
        require(ledger.below_number > 0.0 && ledger.above_number > 0.0, label);
}

void check_alpha_result(const AlphaResult &result, const std::vector<double> &edges,
                        double available, const char *label, bool require_spill = false) {
    R mapped_number = 0.0L, mapped_energy = 0.0L;
    for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
        require(std::isfinite(result.birth[i]) && result.birth[i] >= 0.0, label);
        const R center = (static_cast<R>(edges[i]) + edges[i + 1]) / 2.0L;
        mapped_number += result.birth[i];
        mapped_energy += center * result.birth[i];
    }
    require(close_value(mapped_number, result.ledger.mapped_number, 4e-12L, 2e-14L), label);
    require(close_value(mapped_energy, result.ledger.mapped_energy_J, 4e-12L, 2e-30L), label);
    check_common_inventory(edges, available, result.ledger, label, require_spill);
    require(result.ledger.alpha0_fraction == 1.0 && result.ledger.low_alpha1_fraction == 0.0 &&
                result.ledger.broad_alpha1_fraction == 0.0,
            label);
}

fusion_three_body_cm_v1 canonical_alpha_event(double available, double q, double cosine = 1.0) {
    fusion_nuclear_mass_v1 alpha{};
    require(fusion_c_nuclear_mass(FUSION_HELIUM4, &alpha) == PB11_STATUS_OK &&
                std::isfinite(alpha.mass_kg) && alpha.mass_kg > 0.0,
            "canonical alpha mass");
    fusion_three_body_cm_v1 event{};
    require(fusion_c_three_equal_sequential_cm(alpha.mass_kg, available, q, cosine, &event) ==
                PB11_STATUS_OK,
            "canonical alpha sequential event");
    return event;
}

void compare_alpha_to_packet_map(const AlphaResult &analytic, double available, double q,
                                 const std::vector<double> &edges) {
    const fusion_three_body_cm_v1 event = canonical_alpha_event(available, q);
    const double low = std::min(event.kinetic_energy_J[1], event.kinetic_energy_J[2]);
    const double high = std::max(event.kinetic_energy_J[1], event.kinetic_energy_J[2]);
    constexpr int secondary_packets = 8192;
    std::vector<double> energies(static_cast<std::size_t>(secondary_packets) + 1);
    std::vector<double> rates(energies.size());
    energies[0] = event.kinetic_energy_J[0];
    rates[0] = 1.0;
    require(high > low, "nondegenerate alpha0 box for midpoint reference");
    for (int i = 0; i < secondary_packets; ++i) {
        energies[static_cast<std::size_t>(i) + 1] =
            low + (i + 0.5) * (high - low) / secondary_packets;
        rates[static_cast<std::size_t>(i) + 1] = 2.0 / secondary_packets;
    }
    std::vector<double> expected(edges.size() - 1, 0.0);
    fusion_birth_mapping_v1 mapped{};
    require(fusion_c_map_birth_packets(static_cast<int>(expected.size()), edges.data(),
                                       static_cast<int>(energies.size()), energies.data(),
                                       rates.data(), expected.data(), &mapped) ==
                PB11_STATUS_OK,
            "dense midpoint packet map");

    R l1 = 0.0L;
    for (std::size_t i = 0; i < expected.size(); ++i)
        l1 += std::abs(static_cast<R>(expected[i]) - analytic.birth[i]);
    // The reference deliberately uses a finite midpoint mesh.  The tolerance
    // is loose compared with the 8192-point quadrature, while still detecting
    // a wrong box interval, particle count, or center projection.
    require(l1 / 3.0L < 2e-4L, "analytic alpha0 box agrees with dense midpoint map");
    require(close_value(mapped.mapped_number_m3_s, analytic.ledger.mapped_number, 4e-4L,
                        3e-7L),
            "analytic alpha0 mapped count agrees with midpoint map");
    require(close_value(mapped.mapped_energy_W_m3, analytic.ledger.mapped_energy_J, 4e-4L,
                        3e-30L),
            "analytic alpha0 mapped energy agrees with midpoint map");
    require(close_value(mapped.below_number_m3_s, analytic.ledger.below_number, 4e-4L,
                        3e-7L) &&
                close_value(mapped.above_number_m3_s, analytic.ledger.above_number, 4e-4L,
                            3e-7L),
            "analytic alpha0 spill count agrees with midpoint map");
    require(close_value(mapped.below_energy_W_m3, analytic.ledger.below_energy_J, 4e-4L,
                        3e-30L) &&
                close_value(mapped.above_energy_W_m3, analytic.ledger.above_energy_J, 4e-4L,
                            3e-30L),
            "analytic alpha0 spill energy agrees with midpoint map");
}

void alpha0_conservation_and_shapes() {
    constexpr double available = 8.68 * mev;
    const double q = 0.35 * available;
    const fusion_three_body_cm_v1 event = canonical_alpha_event(available, q);
    const double minimum = std::min({event.kinetic_energy_J[0], event.kinetic_energy_J[1],
                                     event.kinetic_energy_J[2]});
    const double maximum = std::max({event.kinetic_energy_J[0], event.kinetic_energy_J[1],
                                     event.kinetic_energy_J[2]});
    require(maximum > minimum && minimum > 0.0, "alpha0 event energy interval");

    const std::vector<std::vector<double>> grids = {
        uniform_edges(24, 0.0, available),
        {0.0, .07 * available, .16 * available, .29 * available, .45 * available,
         .64 * available, .81 * available, available},
        // Move both arithmetic-center boundaries inside the physical range;
        // all excluded source weight must remain in the explicit spill ledger.
        uniform_edges(12, minimum + .25 * (maximum - minimum),
                      maximum - .25 * (maximum - minimum))};
    for (std::size_t i = 0; i < grids.size(); ++i) {
        const AlphaResult result = evaluate_alpha(available, q, grids[i]);
        check_alpha_result(result, grids[i], available, "alpha0 conservation on varied grid",
                           i == 2);
    }

    // The box has a useful width for a direct independent midpoint check,
    // including a grid whose center hull clips both ends of the source.
    compare_alpha_to_packet_map(evaluate_alpha(available, q, grids[2]), available, q, grids[2]);

    // One cell is a valid grid.  Choose its sole center between the event's
    // minimum and maximum so that both sides exercise explicit clipping.
    const std::vector<double> one_cell = {minimum / 2.0, maximum + minimum / 2.0};
    const AlphaResult one = evaluate_alpha(available, q, one_cell);
    check_alpha_result(one, one_cell, available, "single-cell alpha0 conservation", true);
    require(one.birth.size() == 1 && one.birth[0] >= 0.0, "single-cell output exists");
}

void alpha0_degenerate_endpoints() {
    constexpr double available = 8.68 * mev;
    const std::vector<double> edges = uniform_edges(16, 0.0, available);
    for (double q : {0.0, available}) {
        const fusion_three_body_cm_v1 event = canonical_alpha_event(available, q);
        require(event.kinetic_energy_J[1] == event.kinetic_energy_J[2],
                "degenerate q gives equal secondary energies");
        const AlphaResult result = evaluate_alpha(available, q, edges);
        check_alpha_result(result, edges, available, "degenerate alpha0 conservation");
        require(result.ledger.secondary_min_J == result.ledger.secondary_max_J &&
                    result.ledger.primary_alpha0_energy_J == event.kinetic_energy_J[0] &&
                    result.ledger.secondary_min_J == event.kinetic_energy_J[1],
                "degenerate alpha0 diagnostics are exact");

        const double energies[3] = {event.kinetic_energy_J[0], event.kinetic_energy_J[1],
                                    event.kinetic_energy_J[2]};
        const double rates[3] = {1.0, 1.0, 1.0};
        std::vector<double> expected(edges.size() - 1, 0.0);
        fusion_birth_mapping_v1 mapped{};
        require(fusion_c_map_birth_packets(static_cast<int>(expected.size()), edges.data(), 3,
                                           energies, rates, expected.data(), &mapped) ==
                    PB11_STATUS_OK,
                "degenerate direct packet map");
        for (std::size_t i = 0; i < expected.size(); ++i)
            require(close_value(expected[i], result.birth[i], 2e-13L, 2e-14L),
                    "degenerate box is exactly two secondary deltas");
        require(close_value(mapped.below_number_m3_s, result.ledger.below_number, 2e-13L,
                            2e-14L) &&
                    close_value(mapped.above_number_m3_s, result.ledger.above_number, 2e-13L,
                                2e-14L),
                "degenerate spill is exact");
    }
}

void poison(fusion_pb_birth_v1 &out) {
    out.mapped_number = out.mapped_energy_J = out.below_number = out.below_energy_J = 9.0;
    out.above_number = out.above_energy_J = out.number_residual = out.energy_residual_J = 9.0;
    out.alpha0_fraction = out.low_alpha1_fraction = out.broad_alpha1_fraction = 9.0;
    out.primary_alpha0_energy_J = out.secondary_min_J = out.secondary_max_J = 9.0;
}

bool cleared(const fusion_pb_birth_v1 &out) {
    return out.mapped_number == 0.0 && out.mapped_energy_J == 0.0 && out.below_number == 0.0 &&
           out.below_energy_J == 0.0 && out.above_number == 0.0 && out.above_energy_J == 0.0 &&
           out.number_residual == 0.0 && out.energy_residual_J == 0.0 &&
           out.alpha0_fraction == 0.0 && out.low_alpha1_fraction == 0.0 &&
           out.broad_alpha1_fraction == 0.0 && out.primary_alpha0_energy_J == 0.0 &&
           out.secondary_min_J == 0.0 && out.secondary_max_J == 0.0;
}

void invalid_alpha0_clears_outputs() {
    constexpr double available = 8.68 * mev, q = 0.35 * available;
    const double edges[3] = {0.0, available / 2.0, available};
    const double decreasing[3] = {0.0, available / 2.0, available / 3.0};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    double birth[2];
    fusion_pb_birth_v1 out{};

    birth[0] = birth[1] = 7.0;
    poison(out);
    require(fusion_c_pb_alpha0_grid(nan, q, 2, edges, birth, &out) ==
                PB11_STATUS_INVALID_ARGUMENT &&
                birth[0] == 0.0 && birth[1] == 0.0 && cleared(out),
            "NaN alpha0 input clears outputs");
    birth[0] = birth[1] = 7.0;
    poison(out);
    require(fusion_c_pb_alpha0_grid(available, available + 1.0, 2, edges, birth, &out) ==
                PB11_STATUS_OUT_OF_RANGE && birth[0] == 0.0 && birth[1] == 0.0 && cleared(out),
            "out-of-range alpha0 q clears outputs");
    birth[0] = birth[1] = 7.0;
    poison(out);
    require(fusion_c_pb_alpha0_grid(available, q, 2, decreasing, birth, &out) ==
                PB11_STATUS_OUT_OF_RANGE && birth[0] == 0.0 && birth[1] == 0.0 && cleared(out),
            "non-increasing alpha0 grid clears outputs");
    double nan_edges[3] = {0.0, nan, available};
    birth[0] = birth[1] = 7.0;
    poison(out);
    require(fusion_c_pb_alpha0_grid(available, q, 2, nan_edges, birth, &out) ==
                PB11_STATUS_INVALID_ARGUMENT && birth[0] == 0.0 && birth[1] == 0.0 && cleared(out),
            "NaN alpha0 grid clears outputs");
    birth[0] = birth[1] = 7.0;
    poison(out);
    require(fusion_c_pb_alpha0_grid(available, q, 0, edges, birth, &out) ==
                PB11_STATUS_INVALID_ARGUMENT && birth[0] == 7.0 && birth[1] == 7.0 &&
                cleared(out),
            "zero-cell alpha0 input clears ledger");
    birth[0] = birth[1] = 7.0;
    poison(out);
    require(fusion_c_pb_alpha0_grid(available, q, 2, edges, birth, nullptr) ==
                PB11_STATUS_NULL_OUTPUT && birth[0] == 0.0 && birth[1] == 0.0,
            "null alpha0 ledger clears birth output");
    poison(out);
    require(fusion_c_pb_alpha0_grid(available, q, 2, edges, nullptr, &out) ==
                PB11_STATUS_NULL_OUTPUT && cleared(out),
            "null alpha0 birth output clears ledger");
}

struct SpectrumResult {
    std::vector<double> birth;
    fusion_alpha_spectrum_v1 ledger{};
};

SpectrumResult evaluate_spectrum(int mode, int policy, double available,
                                 const std::vector<double> &edges, double cutoff, double k,
                                 double phase, int nq = 16, int ncos = 16) {
    SpectrumResult result;
    result.birth.assign(edges.size() - 1, 0.0);
    require(fusion_c_alpha_spectrum_model_grid(
                mode, policy, available, cutoff, k, phase, nq, ncos,
                static_cast<int>(result.birth.size()), edges.data(), result.birth.data(),
                &result.ledger) == PB11_STATUS_OK,
            "alpha1 spectrum status");
    return result;
}

AlphaResult evaluate_mixture(double available, double q, double f0, double flow, int mode,
                             int policy, double cutoff, double k, double phase,
                             const std::vector<double> &edges, int nq = 16, int ncos = 16) {
    AlphaResult result;
    result.birth.assign(edges.size() - 1, 0.0);
    const int status = fusion_c_pb_cm_source_grid(
                available, q, f0, flow, mode, policy, cutoff, k, phase, nq, ncos,
                static_cast<int>(result.birth.size()), edges.data(), result.birth.data(),
                &result.ledger);
    require(status == PB11_STATUS_OK,
            "CM source mixture status");
    return result;
}

void compare_common(const fusion_pb_birth_v1 &a, const fusion_pb_birth_v1 &b,
                    const char *label) {
    require(close_value(a.mapped_number, b.mapped_number, 4e-11L, 2e-13L) &&
                close_value(a.mapped_energy_J, b.mapped_energy_J, 4e-11L, 2e-30L) &&
                close_value(a.below_number, b.below_number, 4e-11L, 2e-13L) &&
                close_value(a.below_energy_J, b.below_energy_J, 4e-11L, 2e-30L) &&
                close_value(a.above_number, b.above_number, 4e-11L, 2e-13L) &&
                close_value(a.above_energy_J, b.above_energy_J, 4e-11L, 2e-30L) &&
                close_value(a.number_residual, b.number_residual, 4e-9L, 2e-14L) &&
                close_value(a.energy_residual_J, b.energy_residual_J, 4e-9L, 2e-30L),
            label);
}

void compare_mixture_to_spectrum(const AlphaResult &mixture, const SpectrumResult &spectrum,
                                const char *label) {
    require(mixture.birth.size() == spectrum.birth.size(), label);
    for (std::size_t i = 0; i < mixture.birth.size(); ++i)
        require(close_value(mixture.birth[i], spectrum.birth[i], 4e-11L, 2e-13L), label);
    require(close_value(mixture.ledger.mapped_number, spectrum.ledger.mapped_number, 4e-11L,
                        2e-13L) &&
                close_value(mixture.ledger.mapped_energy_J, spectrum.ledger.mapped_energy_J,
                            4e-11L, 2e-30L) &&
                close_value(mixture.ledger.below_number, spectrum.ledger.below_number, 4e-11L,
                            2e-13L) &&
                close_value(mixture.ledger.above_number, spectrum.ledger.above_number, 4e-11L,
                            2e-13L),
            label);
}

void mixture_corners_and_composition() {
    constexpr double available = 9.3 * mev, q = 0.35 * available;
    constexpr double cutoff = .001 * mev, k = .71, phase = .43;
    const std::vector<double> edges = uniform_edges(32, 0.0, available);
    const AlphaResult alpha = evaluate_alpha(available, q, edges);

    const AlphaResult alpha_corner = evaluate_mixture(available, q, 1.0, 0.0, 13, 0,
                                                      cutoff, k, phase, edges);
    for (std::size_t i = 0; i < edges.size() - 1; ++i)
        require(close_value(alpha_corner.birth[i], alpha.birth[i], 4e-11L, 2e-13L),
                "f0=1 corner equals alpha0 API");
    compare_common(alpha_corner.ledger, alpha.ledger, "f0=1 ledger equals alpha0 API");
    require(alpha_corner.ledger.alpha0_fraction == 1.0 &&
                alpha_corner.ledger.low_alpha1_fraction == 0.0 &&
                alpha_corner.ledger.broad_alpha1_fraction == 0.0,
            "f0=1 fractions");

    const SpectrumResult low = evaluate_spectrum(2, 0, available, edges, cutoff, k, phase);
    const AlphaResult low_corner = evaluate_mixture(available, q, 0.0, 1.0, 13, 0,
                                                    cutoff, k, phase, edges);
    compare_mixture_to_spectrum(low_corner, low, "flow=1 corner equals low alpha1 API");
    require(low_corner.ledger.alpha0_fraction == 0.0 && low_corner.ledger.low_alpha1_fraction == 1.0 &&
                low_corner.ledger.broad_alpha1_fraction == 0.0,
            "flow=1 fractions");

    const SpectrumResult broad = evaluate_spectrum(13, 0, available, edges, cutoff, k, phase);
    const AlphaResult broad_corner = evaluate_mixture(available, q, 0.0, 0.0, 13, 0,
                                                      cutoff, k, phase, edges);
    compare_mixture_to_spectrum(broad_corner, broad, "broad=1 corner equals broad alpha1 API");
    require(broad_corner.ledger.alpha0_fraction == 0.0 &&
                broad_corner.ledger.low_alpha1_fraction == 0.0 &&
                broad_corner.ledger.broad_alpha1_fraction == 1.0,
            "broad=1 fractions");

    constexpr double f0 = .31, flow = .27, fbroad = 1.0 - f0 - flow;
    const AlphaResult interior = evaluate_mixture(available, q, f0, flow, 13, 0, cutoff, k,
                                                  phase, edges);
    for (std::size_t i = 0; i < edges.size() - 1; ++i) {
        const R expected = f0 * alpha.birth[i] + flow * low.birth[i] + fbroad * broad.birth[i];
        require(close_value(interior.birth[i], expected, 5e-10L, 3e-13L),
                "interior mixture weighted birth");
    }
    const fusion_pb_birth_v1 &ia = interior.ledger;
    const R expected_mapped_number = f0 * alpha.ledger.mapped_number +
                                     flow * low.ledger.mapped_number +
                                     fbroad * broad.ledger.mapped_number;
    const R expected_mapped_energy = f0 * alpha.ledger.mapped_energy_J +
                                     flow * low.ledger.mapped_energy_J +
                                     fbroad * broad.ledger.mapped_energy_J;
    require(close_value(ia.mapped_number, expected_mapped_number, 5e-10L, 3e-13L) &&
                close_value(ia.mapped_energy_J, expected_mapped_energy, 5e-10L, 3e-30L),
            "interior mixture weighted mapped ledger");
    require(close_value(ia.below_number,
                        f0 * alpha.ledger.below_number + flow * low.ledger.below_number +
                            fbroad * broad.ledger.below_number,
                        5e-10L, 3e-13L) &&
                close_value(ia.above_number,
                            f0 * alpha.ledger.above_number + flow * low.ledger.above_number +
                                fbroad * broad.ledger.above_number,
                            5e-10L, 3e-13L),
            "interior mixture weighted spill ledger");
    require(ia.alpha0_fraction == f0 && ia.low_alpha1_fraction == flow &&
                close_value(ia.broad_alpha1_fraction, fbroad, 2e-15L, 2e-15L),
            "interior mixture fractions");
    check_common_inventory(edges, available, interior.ledger, "interior mixture conservation");
}

void invalid_mixture_controls_even_when_inactive() {
    constexpr double available = 9.3 * mev, q = .35 * available;
    constexpr double cutoff = .001 * mev, k = .71, phase = .43;
    const std::vector<double> edges = uniform_edges(8, 0.0, available);
    double birth[8];
    fusion_pb_birth_v1 out{};
    const auto expect = [&](double f0, double flow, int mode, int policy, double local_cutoff,
                            double local_k, double local_phase, int nq, int ncos, int expected) {
        std::fill(birth, birth + 8, 7.0);
        poison(out);
        require(fusion_c_pb_cm_source_grid(available, q, f0, flow, mode, policy, local_cutoff,
                                           local_k, local_phase, nq, ncos, 8, edges.data(), birth,
                                           &out) == expected,
                "invalid mixture control status");
        for (double value : birth) require(value == 0.0, "invalid mixture clears birth");
        require(cleared(out), "invalid mixture clears ledger");
    };
    // The broad branch is inactive at f0=1, but its mode, FSCI policy, and
    // numerical controls are still part of the validated interface contract.
    expect(1.0, 0.0, 2, 0, cutoff, k, phase, 16, 16, PB11_STATUS_INVALID_ARGUMENT);
    expect(1.0, 0.0, 13, 7, cutoff, k, phase, 16, 16, PB11_STATUS_INVALID_ARGUMENT);
    expect(1.0, 0.0, 13, 0, cutoff, 2.0, phase, 16, 16, PB11_STATUS_OUT_OF_RANGE);
    expect(1.0, 0.0, 13, 0, cutoff, k, phase, 3, 16, PB11_STATUS_INVALID_ARGUMENT);
    expect(-.1, 1.0, 13, 0, cutoff, k, phase, 16, 16, PB11_STATUS_OUT_OF_RANGE);
    expect(.6, .5, 13, 0, cutoff, k, phase, 16, 16, PB11_STATUS_OUT_OF_RANGE);
    expect(std::numeric_limits<double>::quiet_NaN(), 0.0, 13, 0, cutoff, k, phase, 16, 16,
           PB11_STATUS_INVALID_ARGUMENT);
}
}

int main() {
    try {
        alpha0_conservation_and_shapes();
        alpha0_degenerate_endpoints();
        invalid_alpha0_clears_outputs();
        mixture_corners_and_composition();
        invalid_mixture_controls_even_when_inactive();
        std::cout << "PASS: alpha0 analytic source, dense box projection, mixture parity and validation\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
