# D889: thread-safe coupled-trial first call

The mutable lazy timing flag was a data race when distinct zones entered the
coupled C API concurrently for the first time. Replace it with a C++11
function-local constant initialized by a lambda. Environment semantics stay
unchanged: unset/empty disables timing; every nonempty value, including `0`,
enables it. The environment must remain unchanged during evaluation.

Four separate regression processes cover unset, empty, `1`, and `0`.
Each releases eight workers before any coupled call, uses immutable shared
inputs and distinct outputs, then compares every defined output field with a
serial call made after joining. Nonzero electron/ion collision heat and changed
fast populations prevent a vacuous no-op comparison. Input snapshots must
remain unchanged. Enabled timing must report all nine calls.

This removes one prerequisite for zone parallelism; it does not deploy a
parallel host or establish full-library thread safety under every API usage.
No ThreadSanitizer claim is made. Full formal suite results are recorded in
receipt.json after completion. Per-test budget is 12 hours, not a short cutoff.

The active D880 matched cohort keeps the frozen 54b35a3c binary and cache.
This modified library has a distinct identity 87747ef25aab70db3470b62e5681d0481bf66304d6646e8bbfa59f04cd43337c.
Do not reuse or relabel old caches under it.

The current native case has a 12-hour budget, versus approximately 3–4 hours
projected for its 20 ms window. Dependency/resource admission allows 48 hours.
Observation waits do not restart the process. Retain one complete native case
at a time and memory/disk guards.
