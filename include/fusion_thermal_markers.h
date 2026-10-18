#ifndef FUSION_THERMAL_MARKERS_H
#define FUSION_THERMAL_MARKERS_H
#include "fusion_magnetic_push.h"
#ifdef __cplusplus
extern "C" {
#endif
enum { FUSION_THERMAL_ZERO_DRIFT_UNIFORM_VOLUME = 1 };
typedef struct fusion_birth_spatial_node_v1 { double x_m[3],volume_weight_m3; } fusion_birth_spatial_node_v1;
typedef struct fusion_birth_direction_v1 { double direction[3],probability_weight; } fusion_birth_direction_v1;
typedef struct fusion_thermal_marker_v1 {
 fusion_orbit_state_v1 initial;
 double number_weight,kinetic_energy_J;
} fusion_thermal_marker_v1;
/* Deterministic tensor quadrature for ONE accepted thermal charged mapped bin.
 * N is extensive particle amount; K[J] is the production energy-cell center.
 * Model1 explicitly assumes zero drift, one-particle isotropy and uniform birth
 * per supplied physical volume. This is not inferred from a zone-total packet.
 * mass[kg]>0, K>=0, N>=0, volume>0, tolerance in(0,1e-3]. Positive spatial and
 * angular weights must sum to volume and1 within relative tolerance. Unit
 * vectors, first angular moment0 and second moment I/3 are checked with this
 * absolute dimensionless tolerance; they do NOT prove quadrature convergence.
 * Weights are used as supplied (w_volume/volume)*w_direction: no normalization,
 * clipping or omission. Direction vectors are used as supplied, not renormalized.
 * Order is spatial-major/angular-minor. Proper velocity u=gamma*v is Cartesian
 * LAB SI; full solid-angle sampling already includes all gyrophases.
 * Counts1..1e6, product<=1e6; capacity must equal product. All arrays disjoint.
 * All output slots clear at entry when capacity0..1e6. Candidate markers are
 * staged internally then copied only after total N/E validation. Positive
 * values that cannot be represented as double cause failure, not silent loss.
 * No callbacks, host indexing, random numbers, persistent state, orbit tracing,
 * tail reconstruction or loss feedback. Beam marginals and neutron packets are
 * NOT eligible. Caller owns accepted/field epoch and canonical species mass.
 * Below/above physical N/E remain separate unsampled accounts; do not invent a
 * marker at their mean energy. Returned quadrature needs independent spatial,
 * angular and orbit convergence tests before any physical interpretation.
 */
int fusion_c_thermal_markers(int model,double number,double kinetic_J,double mass_kg,
 double source_volume_m3,double tolerance,int spatial_count,const fusion_birth_spatial_node_v1*spatial,
 int direction_count,const fusion_birth_direction_v1*directions,int capacity,fusion_thermal_marker_v1*out);
/* Stochastic compression of ONE validated monoenergy thermal quadrature bin.
 * Counts/draws 1..1e6; input, uniforms and output arrays must be disjoint.
 * All x/u/N/K finite, N>=0, equal K>=0 and extensive sum N>0. Caller owns
 * mass/species/epoch and the original spatial/angular validation; arbitrary
 * inputs do not establish isotropy, and beam marginals are not eligible.
 * External independent uniforms satisfy 0<=q<1; no RNG or hidden state.
 * Half-open inverse CDF uses long-double accumulation. Any positive weight
 * whose cumulative interval disappears is rejected; zero weights are skipped.
 * Finite precision uniforms may still fail to reach very small intervals:
 * this API does NOT guarantee ideal real-number categorical probabilities.
 * Copies sampled x/u/K with replacement and assigns sum(N)/draws to each.
 * N and N*K must close within 64*double epsilon relative to the input totals;
 * otherwise fail (including excessive subnormal rounding). Moments and
 * correlations are not exact.
 * Sampling error, original quadrature error and orbit classification uncertainty
 * are separate caller responsibilities. No tail inference or source feedback.
 * Positive unrepresentable output weights/energies or extensive totals fail.
 * Valid-sized output clears at entry; publication is staged/all-or-nothing.
 * Inputs remain unchanged. Invalid output capacity cannot be cleared safely.
 */
int fusion_c_thermal_marker_resample(int count,const fusion_thermal_marker_v1*input,
 int draws,const double*uniforms,fusion_thermal_marker_v1*output);
#ifdef __cplusplus
}
#endif
#endif
