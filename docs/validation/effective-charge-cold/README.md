# Explicit charge factors and cold-trace numerics

New additive C/Fortran effective-charge entry point leaves the existing physical
charge validation strict. He4 effective factor 4.04 is passed without clamping;
this represents caller-declared legacy regularization, not a physical ion state.

GNU: 61/61 tests excluding the formerly mislabeled coupled rejection fixture,
then the corrected complete coupled test passed. Intel-r8: three native tests
passed; three binding tests initially failed to load libimf.so because oneAPI
was not sourced. After sourcing it, all three passed. The corrected native
coupled test also passed. Strict C11 consumer compilation/link/null-output call
passed. The actual host driver exercises the new Fortran variant nontrivially.

Fail-before tests expose bath-heat and residual underflow; actual cold-host
runs additionally expose internal-transfer ledger underflow and a purely
subnormal T-component. FP now propagates measured state/heat underflow rounding
through its balance operator. No density/energy floor, negative clipping or
relaxed ordinary relative tolerance. The combined S+T balance remains checked.

The former negative-energy test was a hot alpha in a cold bath; it actually
failed from underflow. It now explicitly checks positive heat and energy closure.
The replacement overdraw test runs the same two-component collision/transfer
operator on hot baths and a huge cold alpha inventory, verifies that its finite
heat exceeds both reservoirs, then verifies coupled rejection and cleared outputs.
This intentionally extreme numerical fixture is not a physical EXL plasma.

Host fixed-grid injection and rejection evidence is paired in the BALDUR repo
under docs/exl50u-program/host-driver. No full nuclear-burning host trajectory,
fast radial transport or production restart is certified by this increment.
