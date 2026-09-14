# DD and D–He3 companion source-table comparisons

These are local numerical stress tests, not EXL predictions or a fair fuel comparison. Both use 800 energy cells, 16 steps over 10 ms, the existing common-temperature birth model at nq=ncos8, and table bounds 19.5–23 keV with 0.1% sampled interpolation gates. The direct and table paths have otherwise identical inputs. Exact commands, binary hash, full CSV/logs and per-column full-history differences are retained. No thermal handoff projections occur.

DD evaluates both primary branches (channels1/2); the driver also constructs secondary DT/DHe3 tables for its enabled network. Final table/direct differences are fast energy0.01159%, electron heat0.01051%, network-ion heat0.01144%, carbon heat0.01833% (0.0109775 J/m3), and Q0.01050%. Maximum local energy residual is1.9423e-16.

D–He3 final differences are fast energy0.01744%, electron heat0.01613%, network-ion heat0.01726%, carbon heat0.01829% (0.0129664 J/m3), and Q0.01624%. Maximum local energy residual is1.3629e-16. Full-history electron-heat discrepancy peaks at0.01622%; terminal comparisons do not replace history checks.

The previous direct files predate the additional separate heat columns. All16 rows and every shared column agree exactly with the newly generated direct files. This checks the source-provider refactor independently of interpolated-mode sensitivity. Zero reference values are reported explicitly rather than assigned an arbitrary relative-error floor. These runs do not establish timestep/grid/source-model convergence; those remain separate checks.
