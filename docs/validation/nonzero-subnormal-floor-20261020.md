# Positive subnormal birth-floor moments — 2026-10-20

A coupled DD development calculation exposed a representation error in the
selection of the existing extended-precision birth-floor path. The last accepted
host time was11.0090625ms; a trial at11.010625ms returned OUT_OF_RANGE in the
floor projection. This was an incomplete20ms calculation, not a timeout.

For a helium4 trace source, the physical per-step number and energy densities
were `0xb.3de44d597f08da4p-991` m^-3 and
`0x8.661bc1d8cd4027dp-1078` J/m^3. The physical mapped energy was
`0xd.ef810e082992a52p-1078` J/m^3, so the original energy was in range.
Binary64 energy rounded upward to denorm_min, exceeding the mapped bound.
The old selector detected only moments rounding to zero; this moment remained
strictly positive after conversion.

The fix routes positive subnormal number or energy amounts through the existing
long-double path before binary64 conversion. It also covers a value just below
the smallest normal double that rounds upward to that normal value. Public ABI,
source models, grids, source-discrepancy tolerances, floor displacement and joint
energy-borrowing limits remain unchanged. No amount is replaced by zero.

Validation:

- A real captured-call replay using the original numerical library reproduces
  status3; an isolated precision-routing prototype returns0 on the same packet.
- The added wrapped-provider regression fails with the original library and
  passes with the correction. Its deliberately inconsistent synthetic nuclear
  yield still causes the later NUMERICAL_FAILURE, clearing every output; it is
  a representation/error-path test, not a physical source model.
- The measured raw-moment witness verifies physical range and reconstructed
  energy correction through the existing precise floor projection.
- A fresh Release build, including Intel ifx Fortran bindings and diagnostics,
  succeeds. All27 selected floor, kinetic, coupled, restart and Fortran tests pass.
  Build uses one worker; each stage permits30days. Source snapshot1331files.
- New genuine cache identity:
  `559bdb1ba34892c17024d3ae6a52fbc5481111fc204b69c4943a718de93907fc`.
  Both thermal and beam APIs report it and reject old df95 packed tables with
  status2 and null handles. No old identity is substituted into the new library.

Coordinator evidence is kept in the BALDUR feature branch under
`docs/exl50u-program/checkpoint-current-threefuel-capture-20261019/terminal-dd-failure`
and `checkpoint-nonzero-subnormal-floor-fix-20261020`.
These gates qualify the local numerical correction and library build. Fresh
cache/host pairing and full-window native acceptance remain required. They do
not establish full1s four-fuel convergence or complete the device campaign.
