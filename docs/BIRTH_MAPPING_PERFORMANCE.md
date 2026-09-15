# Reuse boost coefficients within one product emission

The source mapper now computes beta-dependent gamma and cancellation-safe
(gamma-1) once per `emit_products` call and shares those unchanged values
across all emitted products/events. It does not change quadrature nodes,
weights, event energies, source normalization or accumulation order. No cache
or mutable state survives the call; C/Fortran interfaces are unchanged.

Validation in the BALDUR coordinator D131 evidence:

- One captured high-energy pB kernel: all22776returnedbytes exact before/after.
- Four alternating process pairs: median3.177086909s original vs3.044571902s
  optimized (4.17percent reduction). This is one kernel on a shared machine,
  not a whole-discharge acceleration claim.
- Full captured400cell coupled request with143unique direct kernels:
  all40276returnedbytes exact against the frozen unmodified reference.
- Existing beam and thermal birth regression executables both pass. Formal
  GNU CMake build and ten source/table/Fortran CTests pass.

Diagnostic inclusive timers found52272emissions and21326976box mappings for
that single high-energy kernel. Timer overhead is material; nested recorded
times are not additive native-time fractions. Repeated uninstrumented paired
measurements above establish the limited performance benefit.

This source edit changes the conservative beam-cache kernel identity. Earlier
cache files may therefore be rejected despite unchanged numerical outputs.
Retain frozen-library provenance for archived cached-run comparisons.
