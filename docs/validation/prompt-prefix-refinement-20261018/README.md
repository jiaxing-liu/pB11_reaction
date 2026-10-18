# D886: continue proving the inside prefix for a retained ancestor candidate

Root traced the actual failed segments (D884). The old traversal stopped when
a small leaf reached requested fractional width even if the retained safe
outside ancestor was still much wider. Its final search could shrink the
outside endpoint but could not prove/advance the discarded inside prefix.

The correction keeps subdividing only in that case. It preserves left-first
proof order,safe inside/outside witnesses,ambiguity,no-outside-field policy,
52-depth stack and exact caller traversal budgets. No geometry allowance,
polygon resolution,physical tolerance or numerical source threshold changed.
The analytic circle test now requires meeting0.001 fractional tolerance;
its1e-5 below-band query remains explicitly unresolved. This strengthens the
regression instead of retaining the prior avoidable wide-candidate outcome.

Isolated existing segment,analytic prompt-orbit and prepared-boundary tests
pass. All four historical unresolved markers become stable model-boundary
EVENT at level4. A separate stricter min_levels7 run returns the same four
classifications; maximum coarse/fine bracket spread is recorded in summary
and remains below the unchanged1e-12s event tolerance. Both runs preserve the
original number/kinetic weights and have zero outside-field requests. There
is no inferred material-wall or total-source loss fraction.

Root applied this source/test increment alongside D885 prepared-transfer to
the external working tree for D887 complete regression. The future identity
is25e6491c...,while current D880 binaries/caches stay54b35a3c. The staged
geometry-only prototype linked the explicitly copied modified segment source
with the old static library,without invoking any caches; its source hashes
and retained raw outputs prevent confusion with a deployed library revision.
