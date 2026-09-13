# Laursen (2016) angular-correlation audit

**Scope.** This note audits the apparent conflict between Laursen et al., Eur. Phys. J. A **52**, 271 (2016), p. 3 Eq. (2),

\[
W_{l=2}(\theta_2)=1.12+0.80\sin^2(2\theta_2),
\]

and the direct specialization of their p. 3 Eq. (3) followed by the p. 3 Eq. (4) sum over the initial projection. It records the source conditions and the unresolved interpretation; it does not choose a production model.

Sources checked:

- Laursen et al., “Complete kinematical study of the 3α breakup of the 16.11 MeV state in 12C,” DOI [10.1140/epja/i2016-16271-2](https://doi.org/10.1140/epja/i2016-16271-2), local PDF `/tmp/m3-spectrum-sources/laursen-2016.pdf` (SHA-256 `88ac6168d8ca3dfdc9ab0b94e55d62f65f306c8b5ff534bd8fa00be74173c7d9`).
- L. C. Biedenharn and M. E. Rose, “Theory of Angular Correlation of Nuclear Radiations,” Rev. Mod. Phys. **25**, 729 (1953), DOI [10.1103/RevModPhys.25.729](https://doi.org/10.1103/RevModPhys.25.729), local PDF `/tmp/m3-spectrum-sources/biedenharn-rose-1953.pdf` (SHA-256 `6bab32f3325e18129909a7ae5f451f07081c9a2f32ae4e93f48a19a162f7e22c`).

## Independent CG calculation

The standard Condon–Shortley Clebsch–Gordan convention was evaluated independently with the exact Racah factorial sum (integer angular momenta; no floating-point fitting). For the parent setup

\[
\sum_{m=-2}^{2}\left|\langle 2,m;2,0\mid 2,m\rangle
Y_2^m(\theta_2,0)\right|^2,
\]

the exact coefficients are:

| \(m\) | \(\langle2,m;2,0\mid2,m\rangle^2\) |
|---:|---:|
| −2 | \(2/7\) |
| −1 | \(1/14\) |
| 0 | \(2/7\) |
| +1 | \(1/14\) |
| +2 | \(2/7\) |

The order in Laursen Eq. (3) is \(\langle2,0;2,m\mid2,m\rangle\). Interchanging the two \(j=2\) entries contributes \((-1)^{2+2-2}=+1\), so it has the same values here. Using \(x=\cos\theta_2\) and the standard normalized \(Y_2^m\), with the common \(1/\pi\) factor suppressed,

\[
\pi W(x)=\frac{5}{14}-\frac{45}{56}x^2+\frac{45}{56}x^4,
\]

and therefore

\[
\boxed{\frac{W(\theta_2)}{W(0)}
 =1-\frac94x^2+\frac94x^4
 =1-\frac{9}{16}\sin^2(2\theta_2).}
\]

This confirms the independent calculation in `/tmp/check_cg.py`. It is also independently checked by the exact CG orthogonality sums: for each fixed \(m\), summing the squared coefficients over all allowed total \(J\) gives one. The result is convention-independent at the level of the intensity: Condon–Shortley phase changes cannot change these squared terms.

At \(\theta_2=\pi/4\), this direct result is \(W/W(0)=7/16=0.4375\). Treating the quoted Eq. (2) decimals as the displayed values gives

\[
\frac{W_{\rm Eq.(2)}(\theta_2)}{W_{\rm Eq.(2)}(0)}
 =1+\frac57\sin^2(2\theta_2),
\]

which is \(12/7\simeq1.7143\) at \(\theta_2=\pi/4\). Thus the two expressions have opposite curvature and opposite extrema; this is not a CG sign or normalization issue.

## What Laursen actually specifies

The notation block on printed pp. 2–3 defines \(j,j'\) as total angular momenta, \(m,m'\) as projections, and \(l,l'\) as orbital angular momenta. Unprimed quantities belong to \(^{12}{\rm C}\to\alpha_1+{}^8{\rm Be}\); primed quantities belong to the \(^{8}{\rm Be}\to\alpha_2+\alpha_3\) decay. The state is \(J^\pi=2^+\), the excited \(^{8}{\rm Be}\) channel is \(j'=2\), and its two-alpha decay has \(l'=2\). For the first decay the allowed \(l\) values are \(0,2,4\); Eq. (2) assumes \(l=2\). The paper defines \(\theta_2\) as the angle of \(\alpha_2\) relative to \(\alpha_1\), measured in the \(^{8}{\rm Be}\) rest frame.

The placement of the formulas matters:

1. Eq. (2) is in §2.2.1 (“No symmetrisation”), after the statement that a single first-decay orbital momentum is assumed. It is attributed to Ref. [23], Biedenharn–Rose.
2. Eq. (3) is the unsymmetrized-channel amplitude \(f_{1,23}\) used later in the Bose-symmetrized calculation. It contains the CG factor \((l\,m-m';\,j'\,m'\mid j\,m)\), \(Y_l^{m-m'}(\Theta_1,\Phi_1)\), and \(Y_{l'}^{m'}(\theta_2,\phi_2)\).
3. Eq. (4) says to square the three permutation amplitudes and sum over \(m\), i.e. average over initial spin directions. If permutations are omitted, the paper says Eq. (2.3) is recovered (the equation label is a typesetting reference to the preceding unsymmetrized expression).

Laursen does **not** state in the Eq. (2) paragraph an alternative initial-spin density matrix, a special primary direction, a different \(m\) weighting, or an alternate definition of the angle that would account for the sign change. The later experimental discussion explicitly uses an unpolarized initial state for the kinematic distributions.

## Independent source check of the angular-correlation convention

Biedenharn–Rose pp. 742–743, Eqs. (46)–(50), describes the direction–direction correlation for an unperturbed intermediate state. It states that the sums can be made incoherent by taking the first propagation direction as the quantization axis, and gives the pure-multipole expression in terms of the first-transition CG population and the second-transition angular function. The adjacent text identifies this form with pure multipoles and unpolarized radiation. Specializing that prescription to initial \(j_1=2\), intermediate \(j=2\), final \(j_2=0\), and \(L_1=L_2=2\) gives the same structure used above: the first direction sets \(q_1=0\), while the second \(L=2\) intensity is summed over the intermediate \(m\). This supports the parent calculation as a valid specialization of the cited angular-correlation formalism.

This source check does not supply a derivation of Laursen’s numerical \(1.12,0.80\) pair. The general Biedenharn–Rose formulas (including their pp. 739–740 Eq. (40) and pp. 742–743 Eqs. (46)–(50)) therefore establish the convention and the missing assumptions, but do not by themselves identify a corrected Laursen equation without reconstructing every tensor-parameter normalization.

## Findings and unresolved items

**Confirmed.** Under the explicit setup “primary direction is the quantization \(z\)-axis; \(j=j'=l=l'=2\); unpolarized initial projection; no Bose permutations,” the exact result is \(1-(9/16)\sin^2(2\theta_2)\). The calculation is independently reproducible and agrees with the structure of Biedenharn–Rose’s direction–direction prescription.

**Confirmed.** Laursen’s printed Eq. (2) is \(1.12+0.80\sin^2(2\theta_2)\), and its stated local conditions do not specify a different density matrix or angular convention that explains the discrepancy. The difference is therefore substantive at the source level, not an arithmetic phase/normalization mistake.

**Unresolved.** The available Laursen text does not show whether Eq. (2) was transcribed from a different Biedenharn–Rose tensor-correlation normalization, contains an unreported condition, or contains a sign/coefficient error. A bounded check of the Springer article landing page and title/DOI searches found no separate erratum or correction; this is an absence-of-hit result, not proof that no later author communication exists. The paper’s footnote 1 concerns a wrong sign in an earlier R-matrix denominator, not the angular Eq. (2).

**Handoff constraint.** Do not silently replace Eq. (2) with the CG result, or vice versa. Preserve both as separately sourced alternatives until the parent model decision explicitly selects one and records the missing convention/possible source error.


## Root implementation decision

The reusable kernel will implement the explicit CG/spherical-harmonic amplitude
with coherent permutations (Laursen Eq3-4/Kuhlwein Eq1-3), and test that algebra
independently. Eq2 is retained as a documented, incompatible standalone model
for comparison; it is not substituted into a coherent amplitude. This does
not certify the resulting spectrum against experiment. Model validation must
compare reconstructed Dalitz/single-alpha distributions and report the
source discrepancy. No claim of an author-confirmed erratum is made.
