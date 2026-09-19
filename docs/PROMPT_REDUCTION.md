# Prompt one-particle weight reduction

`fusion_c_prompt_reduce` is a pure C-ABI extensive accounting primitive. Each
input supplies number weight N, single-particle kinetic energy K[J], orbit API
status and numerical outcome. Energy weight is derived as N*K. No independent
energy weight can drift away from the number measure. The caller must ensure K
matches its orbit marker's mass and u; this reducer does not receive those fields.
Group compatible source epochs, species/channels and requested horizons before
interpreting a reduction. This is not a sampler or a source normalization rule.

Categories are indexed by the orbit outcome enum: unresolved0, retained1,event2.
For number and energy independently, the classification interval is
[event/total,(event+unresolved)/total]. It describes only supplied marker weights;
it is not a quadrature convergence bound, a confidence interval, a physical
model uncertainty bound, or a material-wall loss fraction. Floating rounding
still applies. Empty or zero-weight ensembles have undefined fractions; their
separate defined flags are0 and numeric0 fields are placeholders, not0%loss.
A positive population with K=0 can have a defined number fraction and undefined
energy fraction.

All callback/API errors propagate and clear the whole result, including errors
on zero-weight markers. They cannot be recast as physical unresolved weight.
Negative/nonfinite inputs and invalid outcomes reject. Category sums use
compensated long-double accumulation; unrepresentable positive category energy,
total or fraction causes numerical failure rather than silent underflow or
clipping. Input order is deterministic, no state is advanced, no allocation/I/O.
Tests cover mixed/pure categories, differing K, zero measures, API errors,
deterministic repetition, invalid data, overflow and representable/unrepresentable
subnormals. C11 header and ASan/UBSan checks accompany CTest.

## Source ownership discovered during integration audit

Existing zero-drift thermal APIs provide one-particle LAB energy marginals under
uniform outgoing rotations. Their isotropic one-particle angular completion may
be used under that declared model, but does not reconstruct correlated reaction
events or pB three-alpha momentum correlations. Beam APIs retain energy marginals
only; a uniform angular completion is not justified. Joint LAB momentum must be
retained from the production integration before beam orbit sampling.

Accepted coupled nuclear-birth totals combine thermal and fast contributions:
per-bin amount is fast_birth + sum_ch(grid_ch/rate_ch * thermal_events_ch).
The physical ledger adds original below-grid birth N/E under the optional floor;
the kinetic source instead places its number at the first center, with numerical
energy correction separately recorded. Above-grid charged spill still rejects.
Do not normalize a thermal shape to mixed species totals or use the floor's
mapped center energy as physical nuclear birth energy. Below-grid mean energy
is insufficient to reconstruct its orbit distribution. Existing accepted CSVs
cannot recover the separate thermal/fast spectral amounts from total moments.
A future additive trial-output API must preserve these separate source packets
without callbacks that mutate state during rejected trials.

Missing spectral/spatial/angular measure and source tails are caller accounting
responsibilities; this reducer cannot infer them. Do not convert an aggregate
unknown spectrum to a single mean-energy orbit. Collisionless prompt applicability
still requires actual slowing/scattering time-scale checks; no FP birth or
escape term is changed by this primitive.
