#pragma once

// Astrophysical S-factor from Tentori & Belloni (2023), equations (2)-(5).
// Input: centre-of-mass energy in MeV. Output: MeV barn.
// The published fit is defined on 0 <= E <= 9.76 MeV; invalid inputs return NaN.
double pb11_sfactor(double E_MeV);

// Fusion cross section from Tentori & Belloni (2023), equation (1).
// Input: centre-of-mass energy in MeV. Output: barn.
// The E -> 0 limit is returned as exactly zero; invalid inputs return NaN.
double pb11_cross_section(double E_MeV);

// Maxwellian thermal reactivity from direct numerical evaluation of equation (6).
// Input: kT in keV. Output: m^3/s. Positive finite temperatures are accepted.
double pb11_reactivity_integral(double T_keV);

// Independent analytic approximation from equations (7)-(15).
// Input: kT in keV. Output: m^3/s.
// The published range is 10 <= kT <= 500 keV; inputs outside it return NaN.
double pb11_reactivity_fast(double T_keV);
