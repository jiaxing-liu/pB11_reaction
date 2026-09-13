# Kuhlwein (2021) \(l_1/l_3\) mixture convention

Source checked: Kuhlwein et al., arXiv:2109.07886v1 (16 Sep 2021), local files
/tmp/m3-spectrum-sources/kuhlwein-2021-preprint.pdf and .txt. Printed page numbers
below refer to the PDF.

## What the paper explicitly says

- On p. 2, Eq. (1), each channel amplitude \(f_c^{m_a}(123)\) contains the
  Clebsch–Gordan factor, two spherical harmonics, the primary reduced-width
  amplitude \(\gamma_c\), the primary penetrability/phase factors, and the
  secondary factor \(F_c(E_{23})\). The channel label \(c\) specifies
  \(\{l_1,\lambda_b\}\) (Table 1, p. 2).
- On p. 2, Eq. (3), the total decay weight is formed by summing over \(m_a\)
  and \(c\) after adding the three alpha permutations coherently.
- Immediately before Eq. (4), p. 2, the paper says contributions with
  \(l_1=1\) and \(l_1=3\) can be added coherently and introduces an adjustable
  relative phase \(\delta\). Eq. (4) is
  \[
  |f|^2=\sum_{m_a}\left|\sqrt{k}\,f_{l_1=1}
       +\sqrt{1-k}\,e^{i\delta}f_{l_1=3}\right|^2.
  \]
- On p. 4 the authors say the upper panels of Fig. 2 show the \(l_1=1\),
  \(l_1=3\), and mixed predictions; Table 2 gives fits, including a fit to
  Eq. (4) with \(\delta=0\). They report that the best model allows both
  angular momenta, with or without the extra phase.
- Table 2 is on p. 5. Its columns are \(l_1=1\) [%], \(l_1=3\) [%],
  \(\chi^2\), and \(\delta\). The rows are:
  \[
  (100,0,22705,\mathrm{n/a}),\quad
  (0,100,22808,\mathrm{n/a}),\quad
  (85(1),15(1),15360,0),\quad
  (76(5),25(5),12297,67(1)\%\cdot2\pi).
  \]
  The table caption says only that the first two columns give the percentage
  of the \(l_1=1\) and \(l_1=3\) components; uncertainties are statistical.
- The p. 5 conclusion repeats that the model with an admixture of \(l_1=1\)
  and \(l_1=3\) describes the data well.

## What is and is not fixed by this text

The displayed \(76(5)\%\) / \(25(5)\%\) row is consistent with reading the
fit coefficient in Eq. (4) as approximately \(k=0.76\), with the \(l_1=3\)
coefficient approximately \(1-k\). The table values are rounded (they sum to
101% as printed), so the paper does not provide a more precise identity
between the tabulated percentages and a numerical \(k\).

The paper does **not** state that each \(f_{l_1}\) is first normalized by a
separate complete three-body phase-space integral. It gives no condition of
the form
\[
\int d\Phi_3\,|f_{l_1}|^2=1,
\]
no component integrals, no normalization constants, and no code or numerical
normalization procedure. It also does not state whether the \(\gamma_c\)
factors already contained in Eq. (1) are held fixed, absorbed into the
definition of each \(f_{l_1}\), or rescaled when Eq. (4) is fitted.

Therefore the source supports only this bounded statement: \(k\) is the
coefficient used in Eq. (4), and Table 2 reports fitted component percentages
associated with that mixture. The source does not establish whether those
percentages were obtained after separately normalizing the two full
three-body amplitudes, after a common simulation normalization, or under
another scale convention. Selecting one of those conventions would go beyond
the paper text.

## Quantitative follow-up

PUBLISHED_SPECTRUM_CHECK.md now compares original vector curves with the
library at128/256/512 quadrature. With declared Laursen radial parameters,
the refined pure-l1/pure-l3/unit-k=.76 mixture differ from the paper's model
curves by half-L1 distances0.0081/0.0461/0.0656 respectively. Reinterpreting
.76 as a raw-amplitude coefficient maps to unit k~.81916 and gives0.0500
for the mixture. Thus basis rescaling alone does not reproduce the plotted
model. The source's missing numeric radial settings must not be conflated
with the normalization ambiguity. These comparisons do not fit selected
detector data or establish an experimental uncertainty interval.
