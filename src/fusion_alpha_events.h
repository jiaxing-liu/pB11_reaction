#ifndef FUSION_ALPHA_EVENTS_H
#define FUSION_ALPHA_EVENTS_H

#include <vector>

namespace fusion_detail {

struct QuadNode {
    double x, w;
};

std::vector<QuadNode> gauss_legendre(int n);

struct AlphaEvent {
    double energy_J[3];
    double weight;
};

struct AlphaEvents {
    std::vector<AlphaEvent> events;
    double normalization_J2 = 0;
    double l1_normalization_J2 = 0;
    double l3_normalization_J2 = 0;
    int pruned_events = 0;
};

int alpha_events(int mode, int policy, double A, double cutoff, double k,
                 double phase, int nq, int ncos, AlphaEvents &out);

} // namespace fusion_detail

#endif
