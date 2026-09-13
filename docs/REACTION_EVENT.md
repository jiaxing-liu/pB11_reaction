# Conservative reaction-parent and laboratory-event bridge

The new C/Fortran interface `fusion_reaction_event.h` connects channel-ordered
reactant momenta to a complete deterministic laboratory event for all five
network channels. It resolves how host kinetic budgets enter product
kinematics. It does not supply reaction/angle probabilities or a thermally
integrated birth source, and it does not activate BALDUR source ownership.

## Two explicit input energy conventions

Masses and Q come from the existing canonical nuclear data API. Inputs are
Cartesian momenta pa,pb in kg m/s, in the channel's reactant-ID order.
The classical convention assigns `Ti=pi²/(2mi)` to each reactant, matching
the energy convention of the transport/rate kernel. These individual inputs
are explicitly NOT on-shell relativistic four-vectors. This is a conservative
extension of a classical kinetic budget to on-shell products, not a claim
that the classical reactants satisfy exact relativistic dynamics.

The on-shell convention instead assigns
`Ti=sqrt(mi²c⁴+pi²c²)-mi c²` for the SAME momenta. It is a separate consistency
and sensitivity choice: replacing classical rate energy moments by these
values without integrating the corresponding selected-energy distribution
would be inconsistent. Neither convention silently changes the rate model.
The diagnostic `classical_minus_on_shell_J` is computed stably as
`sum(Ton_i²/(2mi c²))`; it is a model-convention difference, not heat or a
missing source to deposit. The API imposes no hidden classical-validity
threshold; the caller must quantify its relevance at actual case conditions.

## Aggregate parent and final-state conservation

Let product rest-energy sum be R, chosen input kinetic sum be T, total product
laboratory kinetic budget L=Q+T and total incoming momentum P=pa+pb. The
aggregate has laboratory energy R+L and momentum P. Define

`Delta = L(2R+L)-c²|P|²`,
`A = Delta/[sqrt(R²+Delta)+R]`,
`V = c² P/(R+L)`.

A is the available kinetic energy of the products in the aggregate CM. This
form avoids subtracting two nearly equal rest energies. The classical
relative energy `mu*|pa/ma-pb/mb|²/2` is returned separately for the existing
cross-section convention. In general it is not interchangeable with A-Q
when a chosen classical budget is combined with an exact Lorentz boost.
All values use the same rounded canonical Q returned by the nuclear API.

For DD/DT/DHe3, use the existing exact two-body kinematics at A with canonical
product masses, then boost by V. For pB, use the exact three-equal-particle
sequential event at explicit intermediate q and secondary cosine. Rotate the
CM event to the supplied primary direction/secondary azimuth, then boost.
Products are on shell under BOTH input conventions. Final ledger checks
sum(K_product)=L and sum(p_product)=P, including explicit neutron products.
Unused third output slot is zero for two-product channels. The pB-specific
q/cosine/azimuth arguments must all be zero for those channels.

The coordinate basis is deterministic: t=normalized cross(z,n), using x in
place of z for |nz|>=.9, b=n cross t. Azimuth rotates t toward b. This is a
coordinate convention, not an angular probability or detector acceptance.
The CM energy marginals in fusion_pb_birth.h cannot supply these correlated
angles by themselves. An NR amplitude model also still needs an explicit
mapping to consistent on-shell event coordinates; no NR momenta are silently
passed as valid on-shell inputs.

## Numerical evidence

The focused C++ test checks400 events across all five channels and both
conventions, canonical product masses, product mass shells, independent
energy/momentum sums, stationary limits, unused slots and error clearing.
An independent80-digit calculation checks50 parent cases using direct
invariant subtraction and direct on-shell energy square roots, rather than
the production stable excess formulas. Maximum relative discrepancy is
9.7788e-17. Some numerical stress cases have momenta unsuitable for a
classical physical interpretation; passing those cases validates arithmetic,
not the classical approximation at such energies.

The interface returns all classical/on-shell/chosen reactant energies, the
relative classical energy, CM available energy, boost, total momentum,
product count, selected convention and the lab kinetic budget. This makes
future distribution integration and source-ledger acceptance auditable.

Remaining work is the full reaction-selected joint energy/CM distribution,
explicit source/ghost/angular models and sensitivity, integration of complete
laboratory births, burn/collision/handoff with evolving baths, host atomic
acceptance/I007, and EXL four-fuel0.6nG full-window validation.

Final matrix: GNU/Intel default/Intel-r8 each48/48; C++-only29/29. Installed
C11 and Intel-r8 consumer programs and host Intel-r8 build pass. The initial
misleading-indentation compiler warning was corrected before this matrix;
no test or numerical-reference failure occurred. No new host mode is activated.
