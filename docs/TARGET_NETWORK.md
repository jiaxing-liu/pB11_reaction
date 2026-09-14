# Simultaneous fast-bin / thermal-target competition

`fusion_c_target_network_trial` is the depletion part of the full fast-fuel
composition. It generalizes the existing single-target trial to multiple
thermal targets and parallel reaction channels. It does not construct products,
release Q as heat, integrate collisions, or accept a host state. All those
amounts still need one complete atomic owner.

## Discrete equations and uniqueness

Let edge e connect fast bin i(e) to target j(e), with frozen K_e and selected
target-energy moment M_e. Fast and thermal pools are distinct even for D+D;
there is no identical-pair half factor. Parallel edges share BOTH pools.
For a=dt*K, aggregate a_ij over parallel edges. Backward Euler gives

    F_i' = F_i / (1 + sum_j a_ij B_j')
    B_j' - B_j + sum_i F_i*a_ij*B_j'/(1+sum_k a_ik*B_k') = 0.

Zero old target remains zero. For positive targets write B_j'=exp(y_j).
The residual is the gradient of

    Phi(y) = sum_j [exp(y_j)-B_j*y_j]
             + sum_i F_i*log(1+sum_j a_ij*exp(y_j)).

Its Hessian is diag(B') plus sum_i F_i*[diag(p_i)-p_i*p_i^T], with
p_ij=a_ij*B_j'/(1+sum_k a_ik*B_k') and p_i0=1/(1+sum_k a_ik*B_k').
The bracket is a categorical covariance including the unreacted category0;
the positive diagonal makes the Hessian positive definite. Phi grows at either
infinite end of any positive-target coordinate, so there is one finite minimum.
The resulting equations also give B_j'=B_j/(1+sum_i a_ij F_i') <= B_j.
No inventory clipping or sequential channel priority is used.

Newton solves at most six target log variables, with row scaling, a bounded
log-step and Armijo line search. A diagonal descent direction is available
when the Newton system cannot supply a descent direction. Iteration/allocation/
representability failures return an explicit error; there is no partial result.
The initial positive lower bound is B_j/(1+sum_i a_ij F_i).

## Stiff arithmetic and acceptance

The near-exhaustion residual must not subtract two almost equal unit-sized
reaction fractions and mistake zero for an accurate surviving population.
For p_ij>1/2 the contribution is accumulated as
F_i - F_i'*(1+sum_{k!=j} a_ik B_k'), using compensated sums. The Hessian diagonal
uses p_ij*(p_i0+sum_{k!=j}p_ik), not p_ij*(1-p_ij). Objective differences use
a dominant-category log shift and compensated sums for the same reason.

Convergence requires target residual normalized by old target <=8 double eps
AND gradient/diagonal-Hessian <=2e-14. The latter guards surviving-inventory
accuracy in stiff balanced pools. An initial prototype used only old-inventory
residual and incorrectly returned tiny but unequal survivors at K=1e100;
this failure is retained in validation. Another initial 2e-16 iteration norm
stopped one constructed case at a 2.29e-16 floating residual. The final explicit
8-eps numerical stopping rule addresses stored-double input precision; the
separate final particle/energy conservation gate remains2e-12.

Per-edge amounts dt*K_e*F_i'*B_j' are converted to output double first; all
returned debit ledgers are summed from those actual amounts. Thus weak events
remain measurable when old-minus-new inventory subtraction rounds to zero.
Target energy debit is sum_e loss_e*M_e/K_e. Negative remaining energy fails;
no target population or energy is repaired. Coefficients/backgrounds are frozen,
so callers still need timestep convergence and thermal background consistency.

## Interface ownership

Explicit fast representative energies allow flattened species/energy bins
without BALDUR indices. At most6 targets,600000 fast bins,2000000 edges are
accepted. K=0 requires M=0. Zero-edge trials preserve inputs exactly. All valid
outputs clear on failures, inputs remain unchanged, and repeated trials have
no accumulated state. The C ABI is additive; the previous single-target API
and nuclear-window wrapper retain their behavior.

Birth composition must use K and M from the SAME beam-source quadrature,
then scale each edge spectrum and spills by loss_e/K_e. No spectrum can be
reconstructed from these energy moments alone. Products must not be returned
to thermal fuel immediately; nuclear_born/fast_consumed/thermal_consumed,
neutrons, collisions and handoff remain separate terms in the complete ledger.

## Bounded validation

Independent tests cover analytic one-bin roots, parallel channels, competing
targets, manufactured multibin solutions, permutations, weak events, rejection,
output clearing, deterministic retry, and extreme balanced stiffness.
Root500 manufactured networks (7 fast bins,3 targets) complete with no failures;
max solution relative error3.493e-11, at most16 iterations. Their rounded initial
inventories are constructed from prescribed final values, so this also contains
input-rounding/conditioning effects.

The root stiff analytic sweep K=1e-200,1,1e20,1e100,1e200 preserves both balanced
survivors; e.g. K=1e200 gives F'=B'=1e-100. A continuous symmetric two-target
problem has B(t)=1/(1+2t);40/80/160/320 steps give errors
0.00602293/0.00303133/0.00152072/0.000761639 at t=1, ratios approaching2.

A physical DD-two-branch plus DT-both-projectile-orientations composition
uses the new beam spectra and this depletion operator. It produces
2.1229200419828784e16 events/m3, max product-species residual7.18e-17 and
energy residual1.94e-17 including neutrons/spills and both fuel debits.
The DD branch ratio matches its coefficient ratio within4.54e-17.
This is a full source/depletion accounting test, not production coupled FP,
atomic host acceptance, or a four-fuel EXL convergence claim.

Final selected GNU4/4 and Intel/r8 4/4 suites pass, including the previous
single-target regression and the new Fortran ABI/shape/zero-edge checks.
Standalone installed strict C11 and Fortran consumers pass with exported
pb11::pb11 and pb11::fusion_target_network_fortran targets. Evidence, including
the failed initial stiff experiment, is in docs/validation/target-network/.
