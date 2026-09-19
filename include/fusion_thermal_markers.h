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
#ifdef __cplusplus
}
#endif
#endif
