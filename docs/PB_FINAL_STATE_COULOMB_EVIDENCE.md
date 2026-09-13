# Final-state Coulomb correction: Laursen 2016 and Refsgaard 2018

This memo records the printed equations and explicit radius choices in the
two supplied sources.  It distinguishes an amplitude-level substitution from
the later probability/decay-weight calculation.  The Refsgaard paper studies
the Hoyle-state decay and is a formalism/caveat source; it is not a new
measurement of the p+11B channels.

## Laursen et al. 2016

Source: K. L. Laursen et al., *Eur. Phys. J. A* **52**, 271 (2016),
[arXiv:1604.01244](https://arxiv.org/abs/1604.01244),
[doi:10.1140/epja/i2016-16271-2](https://doi.org/10.1140/epja/i2016-16271-2).
The relevant equations are on printed p. 3 (PDF page 3), in Secs. 2.2.1-
2.2.3.

### What is an amplitude and what is a probability

Laursen defines a particular-permutation sequential decay amplitude
`f_{1,23}` in Eq. (3).  Its numerator contains the primary and secondary
width factors, and the footnote to Eq. (5) states that the penetrability is
included through

\[
\Gamma=2P(E)\gamma^2.
\]

The Bose-symmetrised decay weight is then Eq. (4),

\[
|f|^2=\sum_m\left|f_{1,23}+f_{2,31}+f_{3,12}\right|^2,
\]

after symmetrising, squaring, and averaging over the initial spin directions.
Thus Eq. (5) is printed as a replacement of the penetrability factor **in
the amplitude Eq. (3)**.  The source does not define it as a separate factor
to be multiplied onto the already squared `|f|^2`; the corrected amplitudes
are what enter Eq. (4).  Because Eq. (3) carries the widths under its
amplitude prefactor, the correction must retain this source ordering when
interpreting the formula.

### Exact Laursen Eq. (5)

The printed Eq. (5) is

\[
P_l\!\left(\frac{3}{2}E_1\right)
\;\longrightarrow\;
\left(\frac{E_1}{E_{12}E_{13}}\right)^{1/2}
\frac{P_l\!\left(\frac{3}{2}E_1\right)}
     {\widetilde P_l\!\left(\frac{3}{2}E_1\right)}
\widetilde P'_{l_{12}}(E_{12})\,
\widetilde P'_{l_{13}}(E_{13}) .
\tag{Laursen 5}
\]

Here `E1` is the kinetic energy of the primary α in the `12C` rest frame;
`E12` and `E13` are the relative energies of the primary α with each
secondary α; and `l12`, `l13` are the corresponding α1-α2 and α1-α3 orbital
angular momenta.  `P_l` is the original primary-channel penetrability at the
first-step energy `3E1/2`.  A tilde means that the penetrability is evaluated
at the enlarged channel radius.  The prime on the two tilde factors marks the
secondary pair factors in the source notation.

The construction is physically described as follows.  The primary α and the
intermediate `8Be` first move apart until a separation `r0`; the intermediate
then breaks up, and the primary α tunnels through the combined α1-α2 and
α1-α3 Coulomb barriers to infinity.  Laursen says that the tunnelling
probabilities combine multiplicatively, which motivates the product in Eq.
(5).  This explanatory statement does not change the equation's placement:
the printed replacement is made before the amplitude is squared.

### Laursen radius and angular choices

On printed p. 3, immediately after Eq. (5), Laursen states:

* the tilde penetrabilities are evaluated at an enlarged radius
  `tilde a = r0`;
* for the calculations in that paper, `tilde a = 10 fm`;
* `l12 = l13 = 2` is assumed for the α1-α2 and α1-α3 systems.

The same page gives the ordinary R-matrix channel radii as

\[
a=a_0\left(4^{1/3}+8^{1/3}\right),\qquad
a'=a_0\left(4^{1/3}+4^{1/3}\right),
\]

with `a0 = 1.42 fm`; numerically the primary radius is about `a = 5.1 fm`
and the secondary radius is `a' = 4.5 fm`.  This `a0` is a nucleon-radius
parameter.  It is distinct from the `r0` used in Sec. 2.2.3 for the breakup
separation.

Their naive separation estimate uses asymptotic primary α-`8Be` speed
`v ≈ 0.068c` and the `8Be(2+)` lifetime
`tau ≈ 0.47 × 10^-22 s`, giving `v tau ≈ 9.6 fm`.  Starting at `a = 5.1 fm`
therefore gives `r0 ≈ 15 fm`.  The paper nevertheless adopts
`tilde a = 10 fm` for the present calculations.  These are two distinct
numbers in the source: `10 fm` is the calculation choice, while `15 fm` is a
rough kinematic estimate.

Laursen Table 1 (printed p. 4) confirms that the Coulomb correction of Eq.
(5) is included in their symmetrised sequential models M3 (`l=0`) and M4
(`l=2`).  It is not included in the uncorrected sequential M2 model or the
democratic M1 model.  The table does not introduce a fitted radius or a
radius uncertainty.

## Refsgaard et al. 2018

Source: J. Refsgaard et al., *Phys. Lett. B* **779**, 414-419 (2018),
[doi:10.1016/j.physletb.2018.02.031](https://doi.org/10.1016/j.physletb.2018.02.031).
The amplitude and FSCI replacement are on printed p. 415 (PDF page 2), the
variable-radius construction is on printed p. 416 (PDF page 3), and the
radius sensitivity/caveats are on printed pp. 416-417 (PDF pages 3-4).

### Refsgaard Model I and Eq. (4)

Refsgaard Eq. (1) defines a sequential amplitude containing

\[
\gamma_c\left(\frac{2P_{l_1}}{\rho_1}\right)^{1/2}
\exp[i(\omega_{l_1}-\phi_{l_1})]F_c(E_{23}),
\]

and Eq. (3) forms the total decay weight by Bose symmetrisation and modulus
square.  Their Model I includes the ordinary primary and secondary
penetrabilities.  They state that this is correct only if the `alpha1+8Be`
and `alpha2+alpha3` pairs can propagate to infinity in their relative
coordinates.  A short-lived intermediate `8Be` invalidates that picture.

Their finite-lifetime replacement, Eq. (4), is printed as

\[
\frac{P_{l_1}}{\rho_1}
\;\longrightarrow\;
\frac{P_{l_1}}{\rho_1}
\left[
\frac{\widetilde\rho_1\,
      \widetilde P_{l_2}(E_{12})\,
      \widetilde P_{l_2}(E_{13})}
     {\widetilde P_{l_1}\,
      \widetilde\rho_{12}\,\widetilde\rho_{13}}
\right].
\tag{Refsgaard 4}
\]

The tilde functions are the usual R-matrix functions evaluated at a chosen
separation `tilde r`; `Eij` is the relative energy of αi and αj.  Refsgaard
explicitly says that this replaces the penetration factor of the α1-`8Be`
pair by the product of factors for α1-α2 and α1-α3, treating each α pair
symmetrically.  Since Eq. (4) is a substitution in their amplitude Eq. (1),
it is likewise applied before their Eq. (3) modulus square.  The authors
mention a different modification in their Ref. [29], but say its
interpretation as transmission probabilities is less clear than Eq. (4).

They call the constant-radius version Model II.  In the Table 4 calculations
they use `tilde r = 16 fm`.  In the surrounding discussion they describe an
average around `15 fm` as a reasonable value for higher-lying decays, and a
footnote says that the `tilde r = 10 fm` quoted in earlier Refs. [16, 27] is
too small because of a calculation error; better agreement is obtained with
somewhat larger `tilde r`.

### Refsgaard variable-radius Model III

Refsgaard obtains a variable separation from the intermediate lifetime.  Their
Eq. (5) is

\[
\tau_2=\hbar\frac{d\delta_2}{dE_{23}}+\frac{a_2}{v_{23}},
\tag{Refsgaard 5}
\]

where `delta2` is the αα scattering phase shift, `a2` is the secondary
channel radius, and `v23` is the secondary relative velocity, approximated by
its asymptotic value.  Their Eq. (6) then gives

\[
\widetilde r=a_1+v_1\tau_2 .
\tag{Refsgaard 6}
\]

The source chooses `a1 = 5.1 fm` and `a2 = 4.5 fm`, corresponding to the
same `r0 = 1.42 fm` nucleon-radius parameter used to construct the channel
radii.  The variable-radius model assumes the relative kinetic energy inside
the Coulomb barrier can be represented by its asymptotic value.  Refsgaard
explicitly warns that this classical kinetic energy is not well defined while
tunnelling.

For the Hoyle-state example through `8Be(0+)`, they say that on resonance the
intermediate state can travel about `10^6 fm` before breaking up, so Model I's
infinite-propagation approximation is expected to work there.  For the
higher-lying `12.71 MeV 1+` example through `8Be(2+)`, the variable `tilde r`
varies around an average of about `15 fm`; this is the context for the
constant `15-16 fm` Model II comparison.  These radius statements are model
comparisons, not direct radius measurements.

### Refsgaard radius sensitivity and failure modes

Refsgaard gives several explicit sensitivity checks for the Hoyle-state
calculation:

* With `a2 = 4.5 fm`, only `57(2)%` of the `8Be` ground-state strength lies
  in the observed narrow ground-state peak.  With `a2 = 7.0 fm`, the peak
  area rises to `86(1)%` (printed p. 416).
* Varying channel radii from `1.42 fm` to `2 fm` changes the calculated
  fractional intensity outside the ground-state peak by about `±10%`.
* Model III fails for the Hoyle-state decay even though it is based on the
  R-matrix formalism.  Refsgaard attributes this in part to the wide
  approximately `35 fm` Coulomb barrier and the use of asymptotic velocity
  in a classically forbidden region; they say the effective tunnelling speed
  may be larger, so the Eq. (6) values of `tilde r` may be underestimated.
  Larger `tilde r` would reduce the three-body Coulomb effect and improve
  Model III's agreement.

The source therefore does not establish one universal final-state Coulomb
radius.  It gives a constant-radius prescription (`16 fm` in its Table 4), a
variable-radius prescription from Eq. (6), and explicit evidence that the
result can be sensitive to channel-radius choices.

## Cross-source evidence boundary

Laursen's p+11B `16.11 MeV` calculation uses Eq. (5) at
`tilde a = 10 fm`, with `l12=l13=2`, while quoting a rough separation estimate
near `15 fm`.  Refsgaard's later analysis identifies the same type of
finite-lifetime correction as a substitution in the decay amplitude, compares
constant `15-16 fm` and variable `tilde r` treatments, and warns that a
10-fm value used in earlier work was too small because of a calculation
error.  Refsgaard also documents that no single radius is generally valid,
especially when the intermediate lifetime and the Coulomb barrier vary across
phase space.

These sources support reporting the equation, its amplitude-before-square
placement, and the stated radius choices.  They do not support treating a
radius as an experimentally measured quantity or selecting one value for all
energies without an additional model decision.
