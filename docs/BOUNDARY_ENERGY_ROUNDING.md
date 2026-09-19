# Monoenergetic boundary outflow in binary64

The lower/upper kinetic energy-domain boundary is an immutable single energy,
not a spectrum average. Its extensive particle count N defines energy N*edge.
Previously energy-work and geometry allowed IEEE casts to underflow, but source
state validation demanded a purely relative energy relation. A valid positive
subnormal N with U rounded to zero therefore failed staging with INVALID_ARGUMENT.

The public42double transport ledger and C/Fortran ABI are unchanged. For a
positive exact product below DBL_MIN only, U is the unique binary64 nearest-even
representation of N*edge. N is never clipped. The helper uses integer binary64
significands and a fixed128bit product (Boost is already a dependency). It is
independent of long-double precision, FTZ, and the active floating rounding mode.
It handles half-quantum ties and tiny products rounding up to DBL_MIN.

Geometry canonicalizes after both component aggregation and extensive-volume
conversion. Source state validates this unique tiny U, then canonicalizes the
cumulative U from the accumulated N. Normal products retain the prior arithmetic
and1e-10 relative validation, without adding an absolute epsilon. Zero N/edge
continues the prior zero contract. Normal1ULP perturbations were not universally
rejected by the old relative contract; this change does not claim otherwise.

Tiny boundary energy in the existing inventory balance uses N*edge in its wider
temporary accumulator instead of interpreting diagnostic U=0 as physical zero.
Those existing balance sums still depend on the supported long-double arithmetic;
they are not a new arbitrary-precision serialized energy accumulator. N and the
immutable edge remain available across restart. No numerical heat is invented.

Validation includes exact midpoint/adjacent-subnormal/min-normal cases in all
four standard rounding modes;18public geometry→stage→commit→pack→unpack scenarios;
wrong adjacent tiny U, false zero, wrong normal U, NaN and negative values;
accumulation from rounded-zero energy to nonzero subnormal and tiny-to-normal;
and existing energy-work, geometry, volume, state and numerical-state tests.
The host evidence additionally compares12,000 deterministic products against
Python exact Fraction conversion, including6,222tiny products. No new physical
convergence or EXL full-window claim follows from this representation repair.
