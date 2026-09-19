# Immutable prepared model-boundary queries

The additive C API prepares a fixed Fourier LCFS once and exposes read-only
point and segment queries. It preserves the existing model-boundary classifier,
error allowances, ambiguity handling and left-first segment traversal. It is
not a material-wall or physical-loss model.

The legacy point implementation remains unchanged and allocation-free. Segment
traversal is factored into one internal template used by both legacy and prepared
queries; the legacy callback still uses the legacy classifier. The prepared
context owns polygon vertices and geometry-only validation results. It stores
long-double vertices/bounds, preserving the original intermediate precision;
on x86-64 the65536-vertex maximum is approximately4MiB, not2MiB. No caller
coefficient pointer, mutable scratch storage, query-dependent distance or winding
is retained. Query distances and winding use the original arithmetic order.

Prepare clears the output handle and publishes only after validation/allocation
succeeds. Query errors clear outputs. NULL destruction is safe. Other handles
must be live library handles; arbitrary/freed pointers and double destruction
are outside the ownership contract. Concurrent read-only queries are supported;
the caller must finish all queries before destruction. Allocation failure returns
EXCEPTION and is not injected in the current tests. No hidden/global cache.

Validation covers circle/ellipse/five-harmonic shapes, inside/outside/near-boundary
queries, segment clear/candidate/unresolved budgets, invalid inputs, maximum
vertex count, original-coefficient mutation after prepare, and concurrent reads.
A separately compiled pre-change source reference verifies1008 point and81
segment results field-by-field bitwise (excluding C-struct padding). Sanitizers,
C11 compilation/link/run and existing boundary/segment/prompt tests pass.
Actual BALDUR boundary point benchmarks also match all return fields bitwise.
Performance timings are workload/platform measurements, not speed guarantees.
The running host/orbit jobs use their frozen binaries and have not switched APIs.
