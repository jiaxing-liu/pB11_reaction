# Finite thermal-fuel reaction trial

`fusion_thermal_burn.h` advances the six thermal species and five event channels
in `fusion_network.h` through a reaction-only backward-Euler stage. Unlike the
instantaneous stoichiometry API, it enforces competition for shared finite fuel
inventories. It is an independently usable component of the burn operator;
product spectra, collisions, heating, thermal handoff and host state ownership
must still be composed before claiming a complete reactor calculation.

## Model and derivation

Let Kc be the supplied frozen reactivity, dt the trial duration, and nc' the
new thermal densities. For channel c consuming reactants a,b,

```
Rc = dt * Kc * na' * nb' / (1 + delta_ab)
ni' = ni_old - sum_c stoichiometric_loss_ic * Rc.
```

Rc is an event AMOUNT in m^-3, not a rate. Both DD branches carry the pair
factor1/2, but each event removes TWO deuterons. Products born in this trial
enter a separate FAST inventory. In particular, DD-born T and He3 are not
instantly put into the same thermal pool and burned again in the same stage.
Previously thermalized T/He3 already present in the old pool can react.

The p/B pair is disjoint from the other thermal reactants. With minor initial
population m, initial difference d, and A=dt*KpB, its surviving minor population
is the cancellation-free positive quadratic root

```
x = 2m / [1+Ad + sqrt((1+Ad)^2 + 4Am)].
```

The other pool remains x+d. This covers equal reactant populations without a
separate subtraction of nearly equal numbers. Zero density/rate is handled
without division by zero.

Write D'=D0*y, aDD=dt*(KDD1+KDD2)*D0, aT=dt*KDT*D0,
aH=dt*KDHe3*D0, bT=T0/D0 and bH=He30/D0. Then

```
y + aDD*y^2 + bT*aT*y/(1+aT*y) + bH*aH*y/(1+aH*y) = 1
T'   = T0/(1+aT*y)
He3' = He30/(1+aH*y).
```

The left side is strictly increasing for y>=0. The unique root is bracketed
by 1/(1+aDD+bT*aT+bH*aH) and the positive DD-only root (or1). The implementation
uses geometric bracket contraction across disparate scales followed by
safeguarded Newton steps. It never advances the channels sequentially against
separate copies of D. Event amounts are computed from the implicit populations,
which preserves tiny reaction signals otherwise lost in `old-new` subtraction.

## Energy and accepted-state ownership

The caller supplies two reaction-conditioned mean reactant energies per event
for each channel, in joules. One possible source is the two-Maxwellian pair
moments in `fusion_beam.h`, divided by their matching nonzero reactivity.
Those means are generally NOT the ordinary3kT/2 Maxwellian mean: the reacting
subset is weighted by sigma*v. Its nuclear-data window diagnostics and model
validity remain the caller's responsibility; this solver cannot certify K.
A disabled channel K=0 requires zero energy means.

Each thermal species loses the sum of its selected reactant energies. Both DD
means debit D. Remaining energy must be nonnegative; otherwise the entire trial
returns NUMERICAL_FAILURE with cleared outputs. The caller reduces the step or
updates/iterates its frozen background. No negative-energy clipping is allowed.
If returned density would round to zero while finite thermal energy remains,
the state is also rejected rather than creating an unusable zero-density bath.

The returned product counts are FAST births, including DD-produced p/T/He3,
DT-produced alpha, DHe3-produced alpha/p, and three pB alpha per event. Neutrons
are separate. The source owner must assign product kinetic energy from removed
reactant energy plus mass-consistent Q, including neutron escape, exactly once.
No fusion Q, deposited heat or thermal ash is added by this stage itself.

All densities and energies returned are candidate states. Repeated calls from
the same accepted old arrays produce the same trial; discarding it is rollback.
Host acceptance must atomically commit the thermal depletion, corresponding
fast births, neutron energy and cumulative counters. This API has no hidden
state, I/O or timestep-dependent rate counters. Coefficients are frozen only
for this first-order substep; evolving-temperature convergence must be checked
by the composing model, not inferred from positivity of a single trial.

## Numerical checks

The ledger uses event amounts rounded once to double precision for product
stoichiometry and energy debits. Species number and energy residuals compare
the actually returned doubles against the accepted old state, with a relative
2e-12 rejection threshold on nonzero inventory/removed-amount scales. Positive
subnormal amounts may round to zero; finite unrepresentable quantities reject.

Root's100-digit mpmath audit independently evaluated the raw returned-state
BE equations, species balances, baryon and charge balances for500 deterministic
random cases over density10^-20..10^80, K10^-70..10^30 and dt10^-30..10^30.
All passed: maximum relative BE-rate residual2.83e-16, species balance1.09e-16,
baryon/charge9.33e-17. This is algebraic/floating-point validation of supplied
coefficients, not a nuclear-data or time-continuum validation. The script and
result are retained alongside the unit-test evidence.

## Interface and time-discretization acceptance

GNU, Intel default and Intel-r8 configurations each pass33/33 CTest cases;
the no-Fortran build passes21/21. The installed C11 example compiles and links
with the exported header/archive. BALDUR Intel-r8 rebuild/link also succeeds;
these operators are still not activated in host source assembly.

Unit tests cover the pB quadratic root, DD branch counting, joint DT/DHe3
competition, the instantaneous limit, finite-step convergence, conditional
energy removal, rejection and unchanged old state. For fixed DD coefficients,
refining dt from.1 to.05s reduces final density error relative to initial D
from.0164939 to.00844893 (ratio1.95219), approaching the independent continuous
solution D(t)=D0/(1+(KDD1+KDD2)*D0*t). The new C++ test's failure checks remain
active in Release builds.

Two initial test-fixture issues were corrected without changing the solver:
(1) the instantaneous stoichiometry helper must receive events/dt before its
rate outputs are multiplied by dt; (2) a guessed5e-6 instantaneous-limit
threshold was too small for the chosen total competing reaction stiffness.
The latter now uses the derived depletion-only bound fa+fb, where
fs=dt*initial_reactant_loss_s/ns, and also requires first-order reduction when
dt is halved. The failed second test log is retained; the first observed
failure was `FAIL: pure pB11 reactant stoichiometry` in the same helper.
No physical or conservation acceptance threshold was loosened.

Reproduce the high-precision algebra audit (mpmath is validation-only):

```
c++ -O2 -std=c++17 -shared -fPIC -Iinclude src/fusion_thermal_burn.cpp \
    src/fusion_network.cpp -o /tmp/libthermal_burn_check.so
python tools/check_thermal_burn.py /tmp/libthermal_burn_check.so
```

Fortran uses `fusion_thermal_burn_fortran` with explicit ONLY imports when
combined with other modules, and exported CMake target
`pb11::fusion_thermal_burn_fortran`. Its public IDs are the C zero-based IDs;
add1 when indexing ordinary Fortran arrays. `examples/thermal_burn_consumer.c`
is a supplied-coefficient demonstration, not an EXL prediction.
