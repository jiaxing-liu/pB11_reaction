# Explicit regular-core surrogate, revision1 (D213)

`fusion_c_flux_core_eval` is a separate opt-in C ABI. It does not modify
`fusion_c_flux_field_eval` or claim to repair the physical equilibrium. Native
quasi-Hermite axis slopes need not satisfy Cartesian regularity. This model uses
axis R/F and first-harmonic slopes plus value/slope constraints at an explicit
zero-based matching node k (1<=k<nodes-1). Compare multiple matching radii.

Let xc=xi[k], t=xi/xc. For R0 and F use f=fa+alpha*t²+beta*t⁴, where
alpha=2*(fc-fa)-xc*f'c/2, beta=xc*f'c/2-(fc-fa).
For harmonic m>=2 use g=t^m*(alpha+beta*t²), beta=(xc*g'c-m*gc)/2,
alpha=gc-beta. C uses this latter formula with m=1. For geometry m=1 use
xc*t*(a+b*t²+d*t⁴), with supplied axis slope a, G=gc/xc,
b=(5G-g'c-4a)/2 and d=(g'c-3G+2a)/2.
Geometry boundary slopes are native; F/C use the right-hand v1 linear slopes.

All axis harmonic values and C(0) must be exactly zero. Axis first-harmonic
mapping must be nondegenerate. For q>=xc the function delegates to v1 unchanged.
At q=0 physical B=(0,0,Fa/Ra) is unique. J and theta derivatives vanish;
reported xi derivatives still depend on the polar direction, as they should.
Matching makes coefficient geometry C1 and B continuous. Field gradients are NOT
generally C1: second geometric derivatives are not matched. This matters for a
future guiding-center interface requiring field gradients.

This explicitly changes near-axis representation; it is not a uniquely determined
physical core. The finite-radius match and normalized coefficients avoid large
xi^-m factors but do not guarantee global nesting or absence of overshoot. A
local fold against axis orientation returns NUMERICAL_FAILURE. Such failures are
never loss or boundary exit. Input/error semantics otherwise follow the headers.
No RZ inverse or marker advancement is present yet.

Validation includes circular/elliptic analytic geometry/field including exact
axis and approach to1e-10, three-harmonic regular polynomial reproduction,
matching continuity, outer bitwise parity, and invalid/degenerate inputs.
The original field-v1 and new core CTest targets both pass; ASan/UBSan passes.
The first CMake build attempted an unknown target before reconfiguration; the
subsequent explicit configure/build runs both tests. Do not count the initial
old one-test CTest output as proof of the new target.

Host D209 snapshot sensitivity at k1/2,3136points each: sampled J positive,
no sampled polygon self-crossing, F/C no sign reversal, outer bitwise parity.
Maximum relative |B| difference from v1 is0.4261%/0.3997%; field-line pitch
integral differences0.4423%/0.4189%. These are model differences, not accepted
physical error bounds or a grid convergence proof. Core self-intersection tests
are finite samples and do not certify global inter-surface nesting. General
shaped-field divergence, inverse-map round trips and orbit sensitivity remain.

Build using an explicitly reconfigured CMake tree; run
`ctest --test-dir <build> -R '^fusion_flux_(core|field)_tests$' --output-on-failure`.
The host archives case-dependent sensitivity scripts/data. Source identities in
new builds change normally; existing frozen campaign executables are untouched.
