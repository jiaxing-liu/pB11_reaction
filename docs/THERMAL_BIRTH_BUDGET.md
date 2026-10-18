# Explicit bounded pB thermal-birth quadrature

`fusion_thermal_birth_budget.h` adds a stateless independent-library entry,
`fusion_c_pb_thermal_birth_budget_grid`, for equal-temperature isotropic
zero-drift proton/B11 Maxwellians. It returns source coefficients in number*m^3/s
on the existing seven-species grid plus explicit below/above number and energy
ledgers. The caller multiplies by n_H*n_B exactly once. It does not evolve fuel,
fast populations, ash, transport or accepted time.

The existing full-source C and Fortran interfaces, their physical options and
default behavior remain available unchanged. Unequal-temperature sources and
beams continue to use their existing explicit APIs; this new entry does not
silently replace them. The new Fortran target is
`pb11::fusion_thermal_birth_budget_fortran`, module of the same unprefixed name.
Its `pb_thermal_birth_budget_grid` wrapper takes `(cells,7)` birth storage.

The policy has `version=1`, `reserved=0`, separate number/energy normalized L1
tolerances in [0,2], and independent coarse nq/ncos between4 and the respective
full orders. Every field is validated, including inactive controls. Either zero
tolerance retains every represented parent at full order and reproduces the
original full-source result. Values2 permit an all-coarse contract test; this
is not a recommended physical accuracy setting. There is no new default.

For example, when the full source options explicitly use orders at least8:

```c
fusion_thermal_birth_budget_policy_v1 policy = {1,8,8,0,1e-4,1e-4};
fusion_thermal_birth_budget_report_v1 report;
int status = fusion_c_pb_thermal_birth_budget_grid(kT_J, &options, &policy,
    cells, edges_J, birth, &out, &report);
```

All positive represented relative-energy/CM parents from the full reference
remain. Parents selected by their reaction and energy weights receive the full
or a cheaper positive conservative decay quadrature, including alpha0's angular
rule. Reaction rates, reactant debits and nuclear Q retain the original scalar
quadrature and traversal order. No source weight is clipped, moment-rescaled or
renormalized, and no population or model-data continuation is silently changed.

If L is the coarse subset, its ideal number and energy normalized L1 bounds
against the same finite full-order reference are

```
b_N = min(2, 2*sum_L(w)/sum_all(w))
b_E = min(2, 2*sum_L(w*E_lab)/sum_all(w*E_lab)).
```

Here E_lab=Q+E_relative+E_CM under this declared classical equal-temperature
parent convention. Positive-mixture total-mass bounds include all grid cells
and both spill ledgers. Report fields expose the two bounds, parent and
relative-node counts, actual full/coarse orders, reaction/energy weights and
whether any parent used the coarse rule. A zero source returns zero weights
and bounds, with valid quadrature metadata.

The separate1e-12 roundoff allowances are engineering allowances verified in
the stated comparisons, not directed interval-arithmetic certificates. The
bounds do not cover the continuous angular integral, finite relative/CM tails,
physical nuclear uncertainty, energy-grid refinement or subsequent nonlinear
plasma response. Every production use still needs appropriate convergence and
physical-model validation. A resolved spectrum at one cold snapshot is not a
whole-trajectory or all-model qualification.

All three output pointers are required. On failure, nonnull output structs are
cleared; nonnull birth clears only for cells1..100000, because other extents
are unknown. Null policy/options/edges reject with INVALID_ARGUMENT, null
required output with NULL_OUTPUT. Unsupported policy version/reserved fields,
nonfinite controls, bad orders/grid and the original physical-domain failures
reject. Arrays and structs must not overlap. Results are staged locally and
published only after complete source and report checks succeed. Calls share
no mutable library globals or hidden sampling/environment policy.

This increment does not enable the policy in birth-table caching, coupled
source evaluation or BALDUR. Future integration must serialize the explicit
accuracy controls and report scope into table/config identities and diagnostics.
The normal source-based kernel identity changes honestly; previously created
cache files must not be relabeled. Running frozen host/library pairs retain
their original identity and remain separate from this new library build.

For an Intel Fortran executable consuming the installed static package, select
its final link driver explicitly, as the library's own Fortran tests do:

```cmake
find_package(pb11_reaction CONFIG REQUIRED)
add_executable(example example.f90)
set_target_properties(example PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(example PRIVATE pb11::fusion_thermal_birth_budget_fortran)
```

The consumer project enables Fortran and CXX and uses a module-compatible
Fortran compiler. Letting CMake choose the CXX driver for an Intel Fortran
object can fail on systems with mismatched default PIE settings; the isolated
installation check retains this failed consumer and the corrected example.
