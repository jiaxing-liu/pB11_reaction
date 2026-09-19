#ifndef FUSION_MAGNETIC_PUSH_H
#define FUSION_MAGNETIC_PUSH_H
#include "pb11_c.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct fusion_orbit_state_v1 { double x_m[3],u_m_s[3]; } fusion_orbit_state_v1;
/* u=gamma*v, not velocity or momentum; Cartesian SI, synchronized x/u.
 * Callback supplies static B[T] at the requested Cartesian midpoint. It must
 * reject invalid/out-of-domain queries BEFORE field evaluation; no exceptions.
 * Callback nonzero status propagates unchanged. It must fill all three B values.
 */
typedef int (*fusion_magnetic_field_callback_v1)(void*context,const double*x_m,double*B_T);
/* Exposes the exact midpoint request used by the wrapper, for preflight drift
 * segment/domain validation. Signed dt, output[3], errors clear output. */
int fusion_c_magnetic_midpoint(const fusion_orbit_state_v1*input,double dt_s,
 double*out_midpoint_m);
/* E=0 collisionless test particle: half drift, exact frozen-B magnetic rotation,
 * half drift. Second-order spatial trajectory, NOT exact full-orbit motion.
 * mass>0, signed charge[C], signed dt[s] (negative supports reversibility).
 * Callback is called once at x+dt*u/(2*gamma), except dt=0 or charge=0.
 * No boundary/segment/endpoint acceptance or adaptive error control is supplied.
 * Caller must validate both drift segments/endpoints and refine dt before using
 * a trial as an accepted orbit or interpreting a boundary event. No loss ledger.
 * out_midpoint is returned only on success. Nonoverlapping inputs/outputs.
 * No hidden state/allocation/I/O; errors clear all available outputs.
 */
int fusion_c_magnetic_push(double mass_kg,double charge_C,double dt_s,
 const fusion_orbit_state_v1*input,fusion_magnetic_field_callback_v1 field,
 void*context,fusion_orbit_state_v1*output,double*out_midpoint_m);
#ifdef __cplusplus
}
#endif
#endif
