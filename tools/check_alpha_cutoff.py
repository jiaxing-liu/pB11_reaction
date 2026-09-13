#!/usr/bin/env python3
"""Reference the error from dropping one Coulomb radial amplitude.

This is a deliberately small, independent numerical probe.  It does not
generate production tables or claim a bound on the full spectrum.  Every
reported maximum is the maximum over the finite A, l, and energy samples
listed in the JSON output.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import mpmath as mp


MP_DPS = 50
mp.mp.dps = MP_DPS

ALPHA_INVERSE = mp.mpf("137.035999084")
ALPHA = 1 / ALPHA_INVERSE
HBARC_MEV_FM = mp.mpf("197.3269804")
ALPHA_MASS_MEV_C2 = mp.mpf("3727.3794118")
PRIMARY_ZPRODUCT = mp.mpf("8")
PRIMARY_MU = 2 * ALPHA_MASS_MEV_C2 / 3
PRIMARY_RADIUS_FM = mp.mpf("5.1")
SECONDARY_ZPRODUCT = mp.mpf("4")
SECONDARY_MU = ALPHA_MASS_MEV_C2 / 2
SECONDARY_RADIUS_FM = mp.mpf("4.5")
Q_RESONANCE_MEV = mp.mpf("3.129")
GAMMA_SQUARED_MEV = mp.mpf("1.075")
WRONSKIAN_TOL = mp.mpf("1e-30")

CUTS_MEV = tuple(mp.mpf(value) for value in ("0.001", "0.002", "0.004", "0.01"))
A_VALUES_MEV = tuple(mp.mpf(value) for value in ("8.68", "8.829", "9.3", "12"))
E_FRACTION_LABELS = ("cut", "cut_over_2", "cut_over_10")
PRIMARY_L = (1, 2, 3)
SECONDARY_L = 2

FAMILY = {
    "primary": {
        "Zproduct": PRIMARY_ZPRODUCT,
        "mu": PRIMARY_MU,
        "radius": PRIMARY_RADIUS_FM,
        "ls": PRIMARY_L,
    },
    "secondary": {
        "Zproduct": SECONDARY_ZPRODUCT,
        "mu": SECONDARY_MU,
        "radius": SECONDARY_RADIUS_FM,
        "ls": (SECONDARY_L,),
    },
}


def mp_text(value: mp.mpf) -> str:
    return mp.nstr(value, 60)


def finite_float(value: mp.mpf, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"{label} is not finite as double")
    return result


def atomic_json_write(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


class CoulombCache:
    """High-precision F/G, fixed-eta derivatives, and U/S/P by energy."""

    def __init__(self) -> None:
        self.values: dict[str, dict[str, dict[int, dict[str, mp.mpf]]]] = {
            "primary": {},
            "secondary": {},
        }
        self.wronskian_checks = 0

    def evaluate(
        self, family_name: str, energy_MeV: mp.mpf
    ) -> dict[int, dict[str, mp.mpf]]:
        energy_MeV = mp.mpf(energy_MeV)
        key = mp_text(energy_MeV)
        cached = self.values[family_name].get(key)
        if cached is not None:
            return cached
        if energy_MeV <= 0:
            raise RuntimeError(f"nonpositive {family_name} energy {energy_MeV}")

        family = FAMILY[family_name]
        Zproduct = family["Zproduct"]
        mu = family["mu"]
        radius = family["radius"]
        ls = family["ls"]
        eta = Zproduct * ALPHA * mp.sqrt(mu / (2 * energy_MeV))
        rho = radius * mp.sqrt(2 * mu * energy_MeV) / HBARC_MEV_FM
        function_ls = range(max(ls) + 2)
        regular = {l: mp.coulombf(l, eta, rho) for l in function_ls}
        irregular = {l: mp.coulombg(l, eta, rho) for l in function_ls}

        result: dict[int, dict[str, mp.mpf]] = {}
        for l in ls:
            factor = mp.sqrt((l + 1) ** 2 + eta**2) / (l + 1)
            slope = (l + 1) / rho + eta / (l + 1)
            F = regular[l]
            G = irregular[l]
            Fp = slope * F - factor * regular[l + 1]
            Gp = slope * G - factor * irregular[l + 1]
            wronskian = Fp * G - F * Gp
            self.wronskian_checks += 1
            if abs(wronskian - 1) > WRONSKIAN_TOL:
                raise RuntimeError(
                    "fixed-eta derivative Wronskian failure "
                    f"family={family_name} l={l} E={energy_MeV} "
                    f"value={wronskian}"
                )
            denominator = F * F + G * G
            P = rho / denominator
            U = mp.sqrt(P / rho)
            S = rho * (F * Fp + G * Gp) / denominator
            result[l] = {"P": P, "S": S, "U": U}
        self.values[family_name][key] = result
        return result


def radial_amplitude(
    cache: CoulombCache, primary_energy: mp.mpf, q_energy: mp.mpf, l: int
) -> mp.mpf:
    primary = cache.evaluate("primary", primary_energy)[l]["U"]
    secondary_values = cache.evaluate("secondary", q_energy)[SECONDARY_L]
    q_P = secondary_values["P"]
    q_S = secondary_values["S"]
    S_resonance = cache.evaluate("secondary", Q_RESONANCE_MEV)[
        SECONDARY_L
    ]["S"]
    real_denominator = (
        Q_RESONANCE_MEV
        - q_energy
        - GAMMA_SQUARED_MEV * (q_S - S_resonance)
    )
    imaginary_denominator = GAMMA_SQUARED_MEV * q_P
    denominator_abs = mp.sqrt(
        real_denominator**2 + imaginary_denominator**2
    )
    return (
        2
        * mp.sqrt(GAMMA_SQUARED_MEV)
        * primary
        * secondary_values["U"]
        / denominator_abs
    )


def sample_one_energy(
    cache: CoulombCache, energy: mp.mpf, orientation: str
) -> dict[str, object]:
    samples: list[dict[str, object]] = []
    for A in A_VALUES_MEV:
        if orientation == "primary_at_E":
            primary_energy = energy
            q_energy = A - energy
        elif orientation == "secondary_at_E":
            primary_energy = A - energy
            q_energy = energy
        else:
            raise ValueError(orientation)
        if primary_energy <= 0 or q_energy <= 0:
            raise RuntimeError(
                f"invalid split A={A} E={energy} orientation={orientation}"
            )
        for l in PRIMARY_L:
            value = radial_amplitude(cache, primary_energy, q_energy, l)
            samples.append(
                {
                    "A_MeV": finite_float(A, "A"),
                    "l_primary": l,
                    "primary_energy_MeV": finite_float(
                        primary_energy, "primary energy"
                    ),
                    "q_energy_MeV": finite_float(q_energy, "q energy"),
                    "radial_abs": finite_float(value, "radial amplitude"),
                }
            )
    maximum = max(samples, key=lambda sample: float(sample["radial_abs"]))
    return {
        "E_MeV": finite_float(energy, "sample energy"),
        "orientation": orientation,
        "max_radial_abs": maximum["radial_abs"],
        "argmax": maximum,
        "samples": samples,
    }


def summarize_trend(rows: list[dict[str, object]]) -> dict[str, object]:
    maximum = max(rows, key=lambda row: float(row["max_radial_abs"]))
    values = [float(row["max_radial_abs"]) for row in rows]
    # rows are intentionally ordered cut, cut/2, cut/10.
    return {
        "sampled_maximum": maximum["max_radial_abs"],
        "sampled_maximum_E_MeV": maximum["E_MeV"],
        "values_in_decreasing_E_order": values,
        "decreases_at_every_sample": all(
            values[index + 1] <= values[index]
            for index in range(len(values) - 1)
        ),
        "last_over_first": values[-1] / values[0] if values[0] else None,
    }


def run() -> dict[str, object]:
    started = time.monotonic()
    cache = CoulombCache()
    cuts: list[dict[str, object]] = []
    for cut in CUTS_MEV:
        energies = (cut, cut / 2, cut / 10)
        orientations: dict[str, list[dict[str, object]]] = {}
        for orientation in ("primary_at_E", "secondary_at_E"):
            rows = [
                sample_one_energy(cache, energy, orientation)
                for energy in energies
            ]
            orientations[orientation] = rows
        cuts.append(
            {
                "cut_MeV": finite_float(cut, "cut"),
                "sample_energies_MeV": [
                    finite_float(energy, "sample energy") for energy in energies
                ],
                "primary_at_E": orientations["primary_at_E"],
                "secondary_at_E": orientations["secondary_at_E"],
                "trend": {
                    "primary_at_E": summarize_trend(
                        orientations["primary_at_E"]
                    ),
                    "secondary_at_E": summarize_trend(
                        orientations["secondary_at_E"]
                    ),
                },
            }
        )

    return {
        "generator": "tools/check_alpha_cutoff.py",
        "mpmath_dps": MP_DPS,
        "parameters": {
            "alpha_inverse": finite_float(ALPHA_INVERSE, "alpha inverse"),
            "hbarc_MeV_fm": finite_float(HBARC_MEV_FM, "hbar c"),
            "alpha_mass_MeV_c2": finite_float(
                ALPHA_MASS_MEV_C2, "alpha mass"
            ),
            "primary": {
                "Zproduct": finite_float(PRIMARY_ZPRODUCT, "primary Z"),
                "mu_MeV_c2": finite_float(PRIMARY_MU, "primary mu"),
                "radius_fm": finite_float(
                    PRIMARY_RADIUS_FM, "primary radius"
                ),
                "l_values": list(PRIMARY_L),
            },
            "secondary": {
                "Zproduct": finite_float(SECONDARY_ZPRODUCT, "secondary Z"),
                "mu_MeV_c2": finite_float(SECONDARY_MU, "secondary mu"),
                "radius_fm": finite_float(
                    SECONDARY_RADIUS_FM, "secondary radius"
                ),
                "l_values": [SECONDARY_L],
            },
            "resonance_energy_MeV": finite_float(
                Q_RESONANCE_MEV, "resonance energy"
            ),
            "gamma_squared_MeV": finite_float(
                GAMMA_SQUARED_MEV, "gamma squared"
            ),
        },
        "formula": {
            "P": "rho/(F^2+G^2)",
            "U": "sqrt(P/rho)",
            "radial_abs": (
                "2*sqrt(1.075)*Up(Ep)*Us(q)/abs("
                "3.129-q-1.075*(S(q)-S(3.129))-i*1.075*P(q))"
            ),
            "orientations": {
                "primary_at_E": "Ep=E, q=A-E",
                "secondary_at_E": "q=E, Ep=A-E",
            },
        },
        "cuts_MeV": [finite_float(cut, "cut") for cut in CUTS_MEV],
        "A_values_MeV": [
            finite_float(value, "A value") for value in A_VALUES_MEV
        ],
        "wronskian_tolerance": finite_float(
            WRONSKIAN_TOL, "Wronskian tolerance"
        ),
        "wronskian_checks": cache.wronskian_checks,
        "cache_entries": {
            family: len(values) for family, values in cache.values.items()
        },
        "elapsed_seconds": time.monotonic() - started,
        "sampled_maximum_warning": (
            "These are maxima over the listed finite E, A, and l samples. "
            "They are not strict mathematical upper bounds and do not "
            "validate a complete alpha spectrum."
        ),
        "cuts": cuts,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    repository = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--output",
        type=Path,
        default=repository / "docs/validation/alpha_cutoff_reference.json",
    )
    args = parser.parse_args()
    payload = run()
    atomic_json_write(args.output, payload)
    print(f"wrote {args.output}")
    print(f"elapsed_seconds={payload['elapsed_seconds']:.3f}")
    print(f"cache_entries={payload['cache_entries']}")
    for cut in payload["cuts"]:
        primary = cut["trend"]["primary_at_E"]["sampled_maximum"]
        secondary = cut["trend"]["secondary_at_E"]["sampled_maximum"]
        print(
            f"cut={cut['cut_MeV']:.6g} "
            f"primary_at_E_max={primary:.12g} "
            f"secondary_at_E_max={secondary:.12g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
