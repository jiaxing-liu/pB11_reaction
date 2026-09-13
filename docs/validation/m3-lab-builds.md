# Laboratory product build evidence

GNU and Intel ifx default/-r8 each pass 19/19; no-Fortran passes 13/13.
Build commands follow m3-products-builds.md, now including laboratory targets.
Intel uses -check all -traceback, with -r8 in the r8 build.
Install no-Fortran build with cmake --install BUILD --prefix PREFIX.
Compile the versioned m3-lab-smoke.c using cc -std=c11 -Wall -Wextra -Werror
-I PREFIX/include; link its object with c++ PREFIX/lib/libpb11.a -lm and run.
The caller checks total CM and boosted lab energies at their actual SI scale.

The separate R-matrix reference is not a production spectrum implementation.
Reproduce with Python and pinned mpmath==1.3.0:
python docs/validation/rmatrix_boundary_reference.py
Explicit numerical constants and 40-digit precision are recorded by the script.
It checks the Coulomb Wronskian, calculates B=S(E0) for Laursen's parameter
convention and estimates the local width by linearizing the real resonance
denominator: Gamma_local=2*gamma_squared*P/(1+gamma_squared*dS/dE).
Two centered finite-difference steps verify the derivative's convergence.
The result 1.477082 MeV matches Refsgaard's tabulated 1.477(13) MeV width.
This validates the parameter/convention mapping, not the full three-alpha
spectrum, interference, incident-energy domain or branching model.
