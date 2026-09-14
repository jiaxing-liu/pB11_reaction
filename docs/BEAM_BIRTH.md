# Fast-projectile / Maxwellian-target product birth

`fusion_c_beam_birth_grid` supplies scalar-energy product source coefficients
for a monoenergetic projectile and a stationary, isotropic Maxwellian target.
The projectile may occupy either canonical reactant slot. This extends the
independent source library; it does not advance fuel, apply densities, choose
self-pair factors, or activate fast-fuel reactions in BALDUR.

## Distribution and energy ownership

Use canonical projectile/target masses mp, mt, reduced mass mu, projectile
speed v=sqrt(2*Ep/mp), and target thermal speed u=sqrt(2*Tt/mt). For relative
speed w, x=w/u and s=v/u, integration of the target Gaussian over relative
velocity directions gives

```
P_w(w) = 4*x^2/(sqrt(pi)*u) * exp(-(x-s)^2)
         * [-expm1(-4*x*s)/(4*x*s)]
kappa = 2*x*s
P(c | w) = kappa*exp(kappa*(c-1))/[-expm1(-2*kappa)]
```

The bracket tends to1 and P(c|w) tends to1/2 at zero beam speed. These are the
same Gaussian distribution and classical conventions used by the independently
implemented beam-rate kernel in BEAM_AND_THERMAL_MOMENTS.md and fusion_beam.cpp.
That kernel analytically integrates conditional energy moments; birth instead
integrates the conditional direction explicitly and reconstructs both momenta.
No mean target velocity or energy replaces this distribution.

For kappa>=1 use z=kappa*(1-c), whose density is exp(-z)/[-expm1(-2*kappa)]
on[0,2*kappa]. The caller supplies an explicit maximum z (8..80), and its
omitted probability is integrated and returned. For kappa<1, the full cosine
interval is integrated directly. This avoids trying to resolve an exponentially
narrow forward peak with unsplit Gauss nodes. Relative-energy knots include
nuclear joins/resonances, population knots and the shifted Gaussian peak.

The reaction coefficient uses sigma(w)*w times THIS joint distribution. Product
birth, target energy debit and projectile debit use the same weights. The fast
particle debit uses the exact caller energy Ep; conversion of its momentum to
double introduces only rounding, checked by the unchanged1e-10 conservation
gate. Both debit-array entries always use canonical channel order, including
when the projectile is slot1. The reference fields are reordered accordingly.

The existing classical aggregate parent, exact product kinematics, pB amplitude
weights and isotropic outgoing-CM event closure are reused. Event preparation,
emission and conservative output assembly are shared with the thermal sources.
There is no new calibrated angular nuclear model or beam-pitch transport.
For the classical input convention the parent invariant satisfies

```
E_parent^2 - c^2*P_parent^2 = M^2*c^4 + mp*mt*c^2*w^2 + K_total^2,
K_total = M*v^2/2 + mt*w^2/2 - mt*v*w*c.
```

Thus the upper retained w and c=-1 bound the available CM energy; the existing
pB12MeV domain is checked conservatively there and again at emitted nodes.

## Cutoffs and cold targets

relative_max_J bounds relative kinetic energy. No missing probability is
renormalized. relative_retained_probability is the radial probability integral;
retained_pair_probability includes actual angular quadrature; and
angular_omitted_pair_probability is the analytically omitted angular probability
integrated over the retained radial range. Numerical sums may slightly exceed1.
These are probability diagnostics, not missing-reaction or product-energy error
bounds. Full matching continuation-model rate/debit discrepancies are also
returned. The nested spectrum.cm_* fields are zero/not applicable.

A zero-temperature target has one exact relative speed and a stationary target.
Its coefficient is sigma*v if it lies inside the supplied relative cutoff;
otherwise the retained source is zero with a nonzero reference/discrepancy where
appropriate. A zero-energy cold projectile gives zero reactions. All options
are still validated. Negative/nonfinite input and unsupported domains reject;
all non-null outputs clear on errors. Valid cell count1..100000; no state or I/O.

## Coupling next

The spectrum coefficients and target debit can feed fusion_c_target_burn_trial:
use K and Mtarget from this SAME source quadrature, then multiply each bin's
birth/spill by accepted_loss/K. Shared-target competition and the atomic
nuclear_born/fast_consumed ledger still need explicit composition. This new API
alone does not complete fast-fuel burn, fast-fast reactions, or the host goal.

## D076 validation

Final GNU selected suite passes 7/7 and Intel/r8 passes 5/5, including cold
DT/pB projectiles in both canonical slots, warm DT selected-target moments,
spills and number/energy conservation, invalid-output clearing, Fortran ABI,
and a two-bin shared-target depletion/birth composition. The complete host
Intel build passes. Four pre-refactor pB outputs (0.002,0.006,0.1,2 keV) retain
byte-identical result structures and all 11200 grid coefficients.

Root angular refinement uses 128 cells on 0..25 MeV, relative order16,
event orders8/8 and maximum angular exponent40. DT Ep50/Tt10 keV gives
8->16 number/energy spectrum L1 differences 8.47566e-7/8.86327e-7;
pB proton Ep675/Tt30 keV gives 2.12489e-7/2.37693e-7. These sums include
all species and spills. At order16 the matching independent rate/debit
errors are below5e-16 at these two points. This is angular refinement only,
not certification of every cutoff, nuclear option or grid.

The two-bin DT composition produces 3.185065884590081e15 events/m3;
fast/target kinetic debits are 46.23858023083423 and7.526826776185522 J/m3.
Number and energy closure residuals are below1.4e-16 relative. It demonstrates
consistent source/debit weights for one shared target, not a production
multi-channel reaction operator or atomic host trial.

Evidence is archived in docs/validation/beam-birth/. Existing outgoing isotropic
closure and nuclear data assumptions remain unchanged; no new angular-data
validation or missing-tail error bound is claimed.
