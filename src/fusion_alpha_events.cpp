#include "fusion_alpha_events.h"
#include "fusion_alpha_amplitudes.h"

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace fusion_detail {
namespace {

constexpr double pi = 3.1415926535897932384626433832795;
constexpr double mev = 1.602176634e-13;
using C = std::complex<double>;

struct AmplitudeEvent {
    std::array<C, 5> a{}, b{};
    double energy[3];
    double measure;
};

double norm(const std::array<C, 5> &a) {
    double s = 0;
    for (C x : a) s += std::norm(x);
    return s;
}

} // namespace

std::vector<QuadNode> gauss_legendre(int n) {
    if (n < 1) throw std::invalid_argument("Gauss-Legendre order must be positive");

    std::vector<QuadNode> a(n);
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double x = std::cos(pi * (i + .75) / (n + .5)), derivative = 0;
        bool converged = false;
        for (int it = 0; it < 64; ++it) {
            double p = 1, prev = 0;
            for (int j = 1; j <= n; ++j) {
                double next = ((2 * j - 1) * x * p - (j - 1) * prev) / j;
                prev = p;
                p = next;
            }
            derivative = n * (x * p - prev) / (x * x - 1);
            double change = p / derivative;
            x -= change;
            if (std::abs(change) < 4e-16) {
                converged = true;
                break;
            }
        }
        if (!converged || !std::isfinite(x) || std::abs(x) >= 1)
            throw std::runtime_error("Gauss-Legendre root did not converge");

        // Recompute derivative at the final abscissa.
        double p = 1, prev = 0;
        for (int j = 1; j <= n; ++j) {
            double next = ((2 * j - 1) * x * p - (j - 1) * prev) / j;
            prev = p;
            p = next;
        }
        derivative = n * (x * p - prev) / (x * x - 1);
        double w = 2 / ((1 - x * x) * derivative * derivative);
        if (!std::isfinite(w) || w <= 0)
            throw std::runtime_error("Invalid quadrature weight");
        a[i] = {-x, w};
        a[n - 1 - i] = {x, w};
    }
    long double weight = 0;
    for (const auto &v : a) weight += v.w;
    if (std::abs(weight - 2) > 1e-12)
        throw std::runtime_error("Quadrature weight closure");
    return a;
}

int alpha_events(int mode, int policy, double A, double cutoff, double k,
                 double phase, int nq, int ncos, AlphaEvents &out) {
    out = {};
    if ((policy != 0 && policy != 1) ||
        (mode != 1 && mode != 2 && mode != 3 && mode != 13) ||
        nq < 4 || nq > 1024 || ncos < 4 || ncos > 1024 ||
        !std::isfinite(A) || !std::isfinite(cutoff) || !std::isfinite(k) ||
        !std::isfinite(phase))
        return PB11_STATUS_INVALID_ARGUMENT;
    if (A <= 0 || A > 12 * mev || cutoff < .001 * mev ||
        cutoff > .01 * mev || k < 0 || k > 1)
        return PB11_STATUS_OUT_OF_RANGE;

    try {
        const auto qnodes = gauss_legendre(nq);
        const auto cnodes = gauss_legendre(ncos);
        std::vector<AmplitudeEvent> events;
        events.reserve(nq * ncos);
        long double n1 = 0, n3 = 0;
        int pruned_events = 0;
        for (const auto &u : qnodes) {
            double theta = (u.x + 1) * pi / 4, s = std::sin(theta), q = A * s * s;
            double dq = A * std::sin(2 * theta) * u.w * pi / 4;
            for (const auto &v : cnodes) {
                AmplitudeEvent e{};
                fusion_alpha_amplitudes_v1 a{}, b{};
                int pa = 0, pb = 0;
                int status = fusion_c_alpha_amplitudes_fsci_cutoff(
                    mode == 13 ? 1 : mode, policy, A, q, v.x, cutoff, &a, &pa);
                if (status) return status;
                if (mode == 13) {
                    status = fusion_c_alpha_amplitudes_fsci_cutoff(
                        3, policy, A, q, v.x, cutoff, &b, &pb);
                    if (status) return status;
                }
                e.measure = a.phase_space_J * dq * v.w;
                for (int M = 0; M < 5; ++M) {
                    e.a[M] = {a.sym_real[M], a.sym_imag[M]};
                    e.b[M] = {b.sym_real[M], b.sym_imag[M]};
                }
                n1 += static_cast<long double>(e.measure) * norm(e.a);
                n3 += static_cast<long double>(e.measure) * norm(e.b);
                if (pa || pb) ++pruned_events;

                // Positive squared nonrelativistic momenta avoid endpoint cancellation.
                double p0 = std::sqrt(4 * (A - q) / 3), star = std::sqrt(q);
                double px = star * std::sqrt((1 - v.x) * (1 + v.x));
                e.energy[0] = p0 * p0 / 2;
                e.energy[1] = (px * px + std::pow(-p0 / 2 - star * v.x, 2)) / 2;
                e.energy[2] = (px * px + std::pow(-p0 / 2 + star * v.x, 2)) / 2;
                events.push_back(e);
            }
        }
        if (!(n1 > 0) || (mode == 13 && !(n3 > 0)))
            return PB11_STATUS_NUMERICAL_FAILURE;

        const C factor = mode == 13
            ? std::sqrt((1 - k) * static_cast<double>(n1 / n3)) *
                  std::polar(1., phase)
            : C(0, 0);
        std::vector<double> rates(events.size());
        long double total = 0;
        for (std::size_t i = 0; i < events.size(); ++i) {
            std::array<C, 5> combined = events[i].a;
            if (mode == 13)
                for (int M = 0; M < 5; ++M)
                    combined[M] = std::sqrt(k) * events[i].a[M] +
                        factor * events[i].b[M];
            rates[i] = events[i].measure * norm(combined);
            total += rates[i];
        }
        if (!(total > 0) || !std::isfinite(total))
            return PB11_STATUS_NUMERICAL_FAILURE;

        out.events.reserve(events.size());
        for (std::size_t i = 0; i < events.size(); ++i) {
            AlphaEvent event{};
            for (int j = 0; j < 3; ++j) event.energy_J[j] = events[i].energy[j];
            event.weight = static_cast<double>(rates[i] / total);
            out.events.push_back(event);
        }
        out.normalization_J2 = static_cast<double>(total);
        out.l1_normalization_J2 = static_cast<double>(n1);
        out.l3_normalization_J2 = mode == 13 ? static_cast<double>(n3) : 0;
        out.pruned_events = pruned_events;
        return PB11_STATUS_OK;
    } catch (...) {
        out = {};
        return PB11_STATUS_EXCEPTION;
    }
}

} // namespace fusion_detail
