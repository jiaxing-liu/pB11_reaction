# Fast-particle / shared-target burn trial

This increment connects nuclear-window beam moments to conservative depletion
of fast energy bins and a shared target pool. It is a reaction-only operator;
product births, nuclear Q, collisions, escape, ash and transport remain owned
by separate operators. It is not the full pB kinetic closure or a host run.

## Discrete equations and positivity

For fixed rate coefficients K_i (m^3/s), fast cell populations N_i (m^-3),
initial target density B and duration dt, backward Euler gives

```
N_i_new = N_i / (1 + dt*K_i*B_new)
B_new + sum_i (N_i - N_i_new) = B.
```

For B>0 set x=B_new/B and a_i=dt*K_i*B. Solve

```
f(x) = x - 1 + sum_i (N_i/B)*a_i*x/(1+a_i*x) = 0.
```

f is strictly increasing, f(0)=-1 and f(1)>=0. Thus the root is unique in
(0,1], and both pools remain nonnegative without clipping. All bins compete
for the same B_new; independently burning each against the full original B
would permit unphysical overconsumption. A positive lower bracket is
max[1/(1+dt*sum_i(K_i*N_i)), 1-sum_i(N_i)/B, 0]. Safeguarded Newton with
geometric bracket subdivision resolves targets surviving across widely
separated scales. Zero target density is an exact no-reaction case.

The implementation uses extended intermediates, validates representability,
and reports balances using the double-precision state actually returned.
Per-bin removed numbers are retained independently of subtraction cancellation
in large inventories. Double rounding/underflow remains subject to the
reported conservation tolerance; there is no artificial population cutoff.

## Energy and ownership

Cell energies E_i are arithmetic edge midpoints. Supplied target moment M_i
is integral sigma*v*E_target (J m^3/s). For R_i events per volume in the step:

```
removed_fast_energy   = sum_i R_i*E_i
removed_target_energy = sum_i R_i*M_i/K_i
U_target_new         = U_target - removed_target_energy.
```

K_i=0 requires M_i=0. M_i/K_i is the reaction-conditioned target energy,
not an automatic 3kT/2. The physical wrapper evaluates K_i and M_i with
fusion_c_beam_maxwellian_window, and initializes target U=3*n*kT/2.
Coefficients are frozen during this first-order trial. A large step can
exhaust target energy before target number; the trial then fails with
NUMERICAL_FAILURE and cleared outputs. The caller must reduce dt or iterate
its background approximation, not repair the energy by clipping it to zero.

The reaction birth owner must account for Q*sum(R_i) plus both removed
reactant energies in its product budget. The removed kinetic energies are
not deposited directly as electron/ion heat. This operator does not establish
a product energy distribution merely by closing that scalar budget.
Fast and thermal pools are distinct, so no identical-pair 1/2 is inserted.
Thermal-thermal identical-pair counting remains in the thermal network;
fast-fast reactions are not included here.

## Interface and trial state

- fusion_c_target_burn_trial accepts explicit coefficient/moment arrays.
- fusion_c_beam_target_burn_window_trial computes them for a monoenergetic
  representative of each fast cell against a Maxwellian target. Explicit
  masses must match the selected nuclear channel.
- Outputs are the trial fast populations, per-bin event counts and a fourteen
  double ledger of both pool inventories, removals and residuals.
- The physical wrapper also returns one full nuclear-window diagnostic per
  cell. A nonzero target temperature has a mathematically incomplete nuclear
  window. OK is a valid WINDOW burn result, not a claim that missing data
  are zero. Unknown probability is not a reaction-rate error bound (I008).

Inputs are never advanced in place. Retrying the same old state returns the
same trial; discarding outputs performs rollback. The caller accepts by
copying the returned state exactly once. No hidden counters, host grids,
BALDUR species slots or filesystem operations occur in this library operator.
The external host still needs a supported time-integrator and restart policy.

## Independent validation

A one-bin N=B=K=dt=1 problem has the backward-Euler quadratic root
(sqrt(5)-1)/2. The test uses this analytic value and energy-conditioned debits.
For the continuous unequal-pool problem N(0)=2, B(0)=1, K=1, the invariant
N-B=1 yields B(t)=exp(-t)/(2-exp(-t)). At t=1 the absolute errors for
40/80/160 accepted steps are 0.0067541/0.00339717/0.00170368. Ratios approach
2, demonstrating first-order convergence. The initial 80-step error exceeded
the frozen 0.003 test threshold; the grid was refined to 160 steps without
relaxing that threshold. This is numerical resolution evidence, not a solver
failure repaired by altering the model.

Additional tests cover multiple bins sharing a scarce target, K=1e200 with a
surplus target, unchanged inactive bins, both particle balances and energy,
energy-overdraw rejection, repeated trial/no input mutation, and cold/warm
pB nuclear-window composition with reaction-selected target energy. Compiler
and installed-interface results are recorded separately in validation/.
