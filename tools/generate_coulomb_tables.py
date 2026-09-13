#!/usr/bin/env python3
"""Generate Chebyshev tables for the finite-radius Coulomb channels.

The production library does not depend on mpmath.  This script is the
high-precision reference generator: it evaluates the Coulomb functions at
40 decimal digits, checks their Wronskian, fits each log-energy segment with
17 Chebyshev-root samples, and independently validates seven non-node points.
Completed initial segments are cached separately so an interrupted run can be
resumed without accepting an incomplete table.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
from pathlib import Path
import tempfile
import time
from typing import Any

import mpmath as mp


MP_DPS = 40
mp.mp.dps = MP_DPS
N_COEFF = 17
MAX_DEPTH = 10
WRONSKIAN_TOL = mp.mpf("1e-25")
LOGP_TOL = mp.mpf("1e-7")
S_TOL = mp.mpf("1e-8")
PHASE_TOL = mp.mpf("1e-7")
VALIDATION_FRACTIONS = (
    mp.mpf("0.07"),
    mp.mpf("0.19"),
    mp.mpf("0.31"),
    mp.mpf("0.43"),
    mp.mpf("0.57"),
    mp.mpf("0.71"),
    mp.mpf("0.89"),
)
REFERENCE_ENERGIES = (
    mp.mpf("0.001"),
    mp.mpf("0.002"),
    mp.mpf("0.01"),
    mp.mpf("0.1"),
    mp.mpf("1.0"),
    mp.mpf("4.0"),
    mp.mpf("12.0"),
)
ENERGY_MIN = mp.mpf("0.001")
ENERGY_MAX = mp.mpf("12.0")
SEGMENT_FACTOR = mp.mpf("2.0")
CACHE_VERSION = 3

ALPHA_INVERSE = mp.mpf("137.035999084")
ALPHA = 1 / ALPHA_INVERSE
HBARC_MEV_FM = mp.mpf("197.3269804")
ALPHA_MASS_MEV_C2 = mp.mpf("3727.3794118")

# Keep this metadata in ordinary Python scalars so it is easy to serialize
# and safe to send to worker processes.  All numerical evaluation below
# converts the constants to mpmath values at MP_DPS precision.
CHANNELS = (
    {
        "channel": 0,
        "family": "alpha_be8",
        "l": 1,
        "Zproduct": "8",
        "mu_MeV_c2": "2*3727.3794118/3",
        "radius_fm": "5.1",
    },
    {
        "channel": 1,
        "family": "alpha_be8",
        "l": 2,
        "Zproduct": "8",
        "mu_MeV_c2": "2*3727.3794118/3",
        "radius_fm": "5.1",
    },
    {
        "channel": 2,
        "family": "alpha_be8",
        "l": 3,
        "Zproduct": "8",
        "mu_MeV_c2": "2*3727.3794118/3",
        "radius_fm": "5.1",
    },
    {
        "channel": 3,
        "family": "alpha_alpha",
        "l": 0,
        "Zproduct": "4",
        "mu_MeV_c2": "3727.3794118/2",
        "radius_fm": "4.5",
    },
    {
        "channel": 4,
        "family": "alpha_alpha",
        "l": 2,
        "Zproduct": "4",
        "mu_MeV_c2": "3727.3794118/2",
        "radius_fm": "4.5",
    },
)

FAMILIES: dict[str, dict[str, Any]] = {
    "alpha_be8": {
        "Zproduct": mp.mpf("8"),
        "mu_MeV_c2": 2 * ALPHA_MASS_MEV_C2 / 3,
        "radius_fm": mp.mpf("5.1"),
        "ls": (1, 2, 3),
    },
    "alpha_alpha": {
        "Zproduct": mp.mpf("4"),
        "mu_MeV_c2": ALPHA_MASS_MEV_C2 / 2,
        "radius_fm": mp.mpf("4.5"),
        "ls": (0, 2),
    },
}


class ReferenceFailure(RuntimeError):
    """Raised when a reference special-function check is not trustworthy."""


def mp_string(value: mp.mpf) -> str:
    return mp.nstr(value, 60)


def atomic_json_write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2)
            stream.write("\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def family_l_values(family: str) -> tuple[int, ...]:
    return FAMILIES[family]["ls"]


def family_constants(family: str) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    spec = FAMILIES[family]
    return spec["Zproduct"], spec["mu_MeV_c2"], spec["radius_fm"]


def energy_boundaries() -> list[mp.mpf]:
    result = [ENERGY_MIN]
    while result[-1] < ENERGY_MAX:
        result.append(min(ENERGY_MAX, result[-1] * SEGMENT_FACTOR))
    return result


def channel_for_l(family: str, l: int) -> int:
    for channel in CHANNELS:
        if channel["family"] == family and channel["l"] == l:
            return int(channel["channel"])
    raise KeyError((family, l))


def real_reference(value: Any, label: str) -> mp.mpf:
    """Convert a Coulomb result to real mp arithmetic without hiding errors."""
    value = mp.mpc(value)
    if abs(mp.im(value)) > mp.mpf("1e-30"):
        raise ReferenceFailure(f"{label} unexpectedly has imaginary part {value}")
    return mp.re(value)


class FamilyEvaluator:
    """Evaluate and cache all requested l values for one physical family."""

    def __init__(self, family: str):
        mp.mp.dps = MP_DPS
        self.family = family
        self.Zproduct, self.mu, self.radius = family_constants(family)
        self.ls = family_l_values(family)
        self.max_l = max(self.ls)
        self.cache: dict[str, dict[int, tuple[mp.mpf, ...]]] = {}
        # The derivative recurrence for l uses l+1.  Evaluating all l through
        # max_l+1 avoids duplicate Coulomb calls for adjacent channels.
        self.function_ls = tuple(range(self.max_l + 2))

    def evaluate(self, energy_MeV: mp.mpf) -> dict[int, tuple[mp.mpf, ...]]:
        energy_MeV = mp.mpf(energy_MeV)
        key = mp_string(energy_MeV)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        if energy_MeV <= 0:
            raise ReferenceFailure(f"nonpositive energy {energy_MeV}")

        eta = self.Zproduct * ALPHA * mp.sqrt(self.mu / (2 * energy_MeV))
        rho = self.radius * mp.sqrt(2 * self.mu * energy_MeV) / HBARC_MEV_FM
        regular = {
            l: real_reference(
                mp.coulombf(l, eta, rho), f"F_l family={self.family} l={l}"
            )
            for l in self.function_ls
        }
        irregular = {
            l: real_reference(
                mp.coulombg(l, eta, rho), f"G_l family={self.family} l={l}"
            )
            for l in self.function_ls
        }

        result: dict[int, tuple[mp.mpf, ...]] = {}
        for l in self.ls:
            F = regular[l]
            G = irregular[l]
            # Fixed-eta rho derivatives from the exact adjacent-l recurrence:
            # y_l' = ((l+1)/rho + eta/(l+1))*y_l
            #        - sqrt((l+1)^2+eta^2)/(l+1)*y_{l+1}.
            factor = mp.sqrt((l + 1) ** 2 + eta**2) / (l + 1)
            slope = (l + 1) / rho + eta / (l + 1)
            Fp = slope * F - factor * regular[l + 1]
            Gp = slope * G - factor * irregular[l + 1]
            wronskian = Fp * G - F * Gp
            if abs(wronskian - 1) > WRONSKIAN_TOL:
                raise ReferenceFailure(
                    "Coulomb Wronskian failure: "
                    f"family={self.family} l={l} E={mp_string(energy_MeV)} "
                    f"value={mp_string(wronskian)}"
                )

            denominator = F * F + G * G
            if denominator <= 0:
                raise ReferenceFailure(
                    f"nonpositive Coulomb denominator family={self.family} "
                    f"l={l} E={mp_string(energy_MeV)}"
                )
            log_probability = mp.log(rho / denominator)
            S = rho * (F * Fp + G * Gp) / denominator
            omega = mp.im(mp.loggamma(l + 1 + 1j * eta))
            phase = omega - mp.atan2(F, G)
            result[l] = (
                log_probability,
                S,
                mp.cos(phase),
                mp.sin(phase),
            )
        self.cache[key] = result
        return result


def chebyshev_nodes(xlo: mp.mpf, xhi: mp.mpf) -> tuple[list[mp.mpf], list[mp.mpf]]:
    midpoint = (xlo + xhi) / 2
    half_width = (xhi - xlo) / 2
    ts = [
        mp.cos(mp.pi * (j + mp.mpf("0.5")) / N_COEFF)
        for j in range(N_COEFF)
    ]
    xs = [midpoint + half_width * t for t in ts]
    return ts, xs


def dct_coefficients(samples: list[mp.mpf]) -> list[mp.mpf]:
    coefficients: list[mp.mpf] = []
    for k in range(N_COEFF):
        total = mp.fsum(
            samples[j]
            * mp.cos(mp.pi * (j + mp.mpf("0.5")) * k / N_COEFF)
            for j in range(N_COEFF)
        )
        coefficient = 2 * total / N_COEFF
        # The generated convention is c0 + sum(k=1..16) c_k T_k.
        coefficients.append(coefficient / 2 if k == 0 else coefficient)
    return coefficients


def clenshaw(coefficients: list[mp.mpf], t: mp.mpf) -> mp.mpf:
    b1 = mp.mpf("0")
    b2 = mp.mpf("0")
    for k in range(N_COEFF - 1, 0, -1):
        b0 = 2 * t * b1 - b2 + coefficients[k]
        b2 = b1
        b1 = b0
    return t * b1 - b2 + coefficients[0]


def validation_metrics(
    evaluator: FamilyEvaluator,
    l: int,
    xlo: mp.mpf,
    xhi: mp.mpf,
    coefficients: list[list[mp.mpf]],
) -> dict[str, mp.mpf]:
    metrics = {
        "logP_abs": mp.mpf("0"),
        "S_abs": mp.mpf("0"),
        "S_scaled": mp.mpf("0"),
        "phase_complex": mp.mpf("0"),
    }
    for fraction in VALIDATION_FRACTIONS:
        x = xlo + (xhi - xlo) * fraction
        t = 2 * (x - xlo) / (xhi - xlo) - 1
        reference = evaluator.evaluate(mp.exp(x))[l]
        approximation = [clenshaw(row, t) for row in coefficients]
        metrics["logP_abs"] = max(
            metrics["logP_abs"], abs(approximation[0] - reference[0])
        )
        s_error = abs(approximation[1] - reference[1])
        metrics["S_abs"] = max(metrics["S_abs"], s_error)
        metrics["S_scaled"] = max(
            metrics["S_scaled"], s_error / (1 + abs(reference[1]))
        )
        metrics["phase_complex"] = max(
            metrics["phase_complex"],
            mp.sqrt(
                (approximation[2] - reference[2]) ** 2
                + (approximation[3] - reference[3]) ** 2
            ),
        )
    return metrics


def metrics_pass(metrics: dict[str, mp.mpf]) -> bool:
    return (
        metrics["logP_abs"] < LOGP_TOL
        and metrics["S_scaled"] < S_TOL
        and metrics["phase_complex"] < PHASE_TOL
    )


def max_metrics(
    left: dict[str, mp.mpf], right: dict[str, mp.mpf]
) -> dict[str, mp.mpf]:
    return {key: max(left[key], right[key]) for key in left}


def fit_family_segment(
    evaluator: FamilyEvaluator,
    xlo: mp.mpf,
    xhi: mp.mpf,
    depth: int,
) -> list[dict[str, Any]]:
    ts, xs = chebyshev_nodes(xlo, xhi)
    node_values = [evaluator.evaluate(mp.exp(x)) for x in xs]
    segment_by_l: dict[int, dict[str, Any]] = {}
    failed: list[tuple[int, dict[str, mp.mpf]]] = []
    for l in evaluator.ls:
        coefficients = [
            dct_coefficients([node_values[j][l][quantity] for j in range(N_COEFF)])
            for quantity in range(4)
        ]
        metrics = validation_metrics(evaluator, l, xlo, xhi, coefficients)
        segment_by_l[l] = {
            "coefficients": coefficients,
            "metrics": metrics,
        }
        if not metrics_pass(metrics):
            failed.append((l, metrics))

    if failed:
        if depth >= MAX_DEPTH:
            details = "; ".join(
                f"l={l} logP={mp_string(m['logP_abs'])} "
                f"Sscaled={mp_string(m['S_scaled'])} "
                f"phase={mp_string(m['phase_complex'])}"
                for l, m in failed
            )
            raise ReferenceFailure(
                f"Chebyshev validation failed at maximum depth {MAX_DEPTH} "
                f"for logE=[{mp_string(xlo)},{mp_string(xhi)}]: {details}"
            )
        midpoint = (xlo + xhi) / 2
        left = fit_family_segment(evaluator, xlo, midpoint, depth + 1)
        right = fit_family_segment(evaluator, midpoint, xhi, depth + 1)
        return left + right

    # Keep the cache in a JSON-friendly, precision-preserving form.  The
    # final include is rounded to the target C++ double separately.
    return [
        {
            "lo": mp_string(xlo),
            "hi": mp_string(xhi),
            "by_l": {
                str(l): {
                    "coefficients": [
                        [mp_string(value) for value in row]
                        for row in segment_by_l[l]["coefficients"]
                    ],
                    "metrics": {
                        key: mp_string(value)
                        for key, value in segment_by_l[l]["metrics"].items()
                    },
                }
                for l in evaluator.ls
            },
        }
    ]


def cache_signature(family: str, index: int, xlo: str, xhi: str) -> dict[str, Any]:
    return {
        "version": CACHE_VERSION,
        "family": family,
        "index": index,
        "xlo": xlo,
        "xhi": xhi,
        "mp_dps": MP_DPS,
        "n_coeff": N_COEFF,
        "max_depth": MAX_DEPTH,
        "validation_fractions": [mp_string(x) for x in VALIDATION_FRACTIONS],
    }


def load_task_cache(
    path: Path, family: str, index: int, xlo: str, xhi: str
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as stream:
            payload = json.load(stream)
    except (OSError, json.JSONDecodeError):
        return None
    if payload.get("signature") != cache_signature(family, index, xlo, xhi):
        return None
    if not isinstance(payload.get("segments"), list) or not payload["segments"]:
        return None
    return payload


def run_segment_task(task: tuple[str, int, str, str, str]) -> dict[str, Any]:
    family, index, xlo_text, xhi_text, cache_text = task
    mp.mp.dps = MP_DPS
    cache_path = Path(cache_text)
    cached = load_task_cache(cache_path, family, index, xlo_text, xhi_text)
    if cached is not None:
        return cached

    evaluator = FamilyEvaluator(family)
    segments = fit_family_segment(
        evaluator, mp.mpf(xlo_text), mp.mpf(xhi_text), depth=0
    )
    payload = {
        "signature": cache_signature(family, index, xlo_text, xhi_text),
        "segments": segments,
    }
    atomic_json_write(cache_path, payload)
    return payload


def mp_metrics_from_segment(
    segment: dict[str, Any], l: int
) -> dict[str, mp.mpf]:
    return {
        key: mp.mpf(value)
        for key, value in segment["by_l"][str(l)]["metrics"].items()
    }


def coefficient_matrix_from_segment(
    segment: dict[str, Any], l: int
) -> list[list[mp.mpf]]:
    return [
        [mp.mpf(value) for value in row]
        for row in segment["by_l"][str(l)]["coefficients"]
    ]


def finite_double(value: mp.mpf, label: str) -> float:
    converted = float(value)
    if not math.isfinite(converted):
        raise ReferenceFailure(f"{label} is not representable as a finite double")
    return converted


def c_double_literal(value: mp.mpf, label: str) -> str:
    converted = finite_double(value, label)
    if converted == 0.0:
        return "0.0"
    return format(converted, ".17g")


def write_include(
    path: Path, channel_segments: dict[int, list[dict[str, Any]]]
) -> None:
    lines = [
        "// Generated by tools/generate_coulomb_tables.py; do not edit.",
        "// Coulomb reference: mpmath 1.3.0, 40 decimal digits.",
        "namespace fusion_table_data {",
        "struct Segment {",
        "    int channel;",
        "    double lo, hi;",
        "    double coeff[4][17];",
        "};",
        "constexpr Segment segments[] = {",
    ]
    ordered: list[tuple[int, dict[str, Any]]] = []
    for channel in sorted(channel_segments):
        for segment in channel_segments[channel]:
            ordered.append((channel, segment))
    for position, (channel, segment) in enumerate(ordered):
        comma = "," if position + 1 < len(ordered) else ""
        lo = mp.mpf(segment["lo"])
        hi = mp.mpf(segment["hi"])
        lines.append(
            f"    {{{channel}, {c_double_literal(lo, 'lo')}, "
            f"{c_double_literal(hi, 'hi')}, {{"
        )
        matrix = coefficient_matrix_from_segment(
            segment, next(
                item["l"] for item in CHANNELS if item["channel"] == channel
            )
        )
        for row_index, row in enumerate(matrix):
            row_comma = "," if row_index + 1 < len(matrix) else ""
            values = ", ".join(
                c_double_literal(
                    value, f"channel={channel} coefficient={row_index},{column}"
                )
                for column, value in enumerate(row)
            )
            lines.append(f"        {{{values}}}{row_comma}")
        lines.append(f"    }}}}{comma}")
    lines.extend(("};", "}  // namespace fusion_table_data", ""))
    atomic_text_write(path, "\n".join(lines))


def atomic_text_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(content)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def reference_points() -> dict[int, list[dict[str, float]]]:
    evaluators = {family: FamilyEvaluator(family) for family in FAMILIES}
    references: dict[int, list[dict[str, float]]] = {}
    for channel in CHANNELS:
        family = str(channel["family"])
        l = int(channel["l"])
        evaluator = evaluators[family]
        points: list[dict[str, float]] = []
        for energy in REFERENCE_ENERGIES:
            values = evaluator.evaluate(energy)[l]
            points.append(
                {
                    "E_MeV": finite_double(energy, "reference energy"),
                    "logP": finite_double(values[0], "reference logP"),
                    "S": finite_double(values[1], "reference S"),
                    "phase_cos": finite_double(values[2], "reference phase cos"),
                    "phase_sin": finite_double(values[3], "reference phase sin"),
                }
            )
        references[int(channel["channel"])] = points
    return references


def write_validation_json(
    path: Path,
    channel_segments: dict[int, list[dict[str, Any]]],
    elapsed_seconds: float,
) -> None:
    references = reference_points()
    channel_payload: list[dict[str, Any]] = []
    for channel in CHANNELS:
        channel_id = int(channel["channel"])
        l = int(channel["l"])
        aggregate = {
            "logP_abs": mp.mpf("0"),
            "S_abs": mp.mpf("0"),
            "S_scaled": mp.mpf("0"),
            "phase_complex": mp.mpf("0"),
        }
        for segment in channel_segments[channel_id]:
            aggregate = max_metrics(
                aggregate, mp_metrics_from_segment(segment, l)
            )
        channel_payload.append(
            {
                "channel": channel_id,
                "family": channel["family"],
                "l": l,
                "Zproduct": finite_double(
                    mp.mpf(str(channel["Zproduct"])), "Zproduct"
                ),
                "mu_MeV_c2": finite_double(
                    family_constants(str(channel["family"]))[1], "reduced mass"
                ),
                "radius_fm": finite_double(
                    family_constants(str(channel["family"]))[2], "radius"
                ),
                "segment_count": len(channel_segments[channel_id]),
                "max_validation_error": {
                    key: finite_double(value, f"channel {channel_id} {key}")
                    for key, value in aggregate.items()
                },
                "thresholds": {
                    "logP_abs": finite_double(LOGP_TOL, "logP threshold"),
                    "S_scaled": finite_double(S_TOL, "S threshold"),
                    "phase_complex": finite_double(PHASE_TOL, "phase threshold"),
                },
                "reference_points": references[channel_id],
            }
        )

    payload = {
        "generator": "tools/generate_coulomb_tables.py",
        "mpmath_dps": MP_DPS,
        "energy_range_MeV": [
            finite_double(ENERGY_MIN, "minimum energy"),
            finite_double(ENERGY_MAX, "maximum energy"),
        ],
        "interpolation": {
            "variable": "log(E_MeV)",
            "segment_factor": finite_double(SEGMENT_FACTOR, "segment factor"),
            "nodes": N_COEFF,
            "coefficient_convention": "c0 + sum(k=1..16) c[k] T[k](t), c0 halved",
            "validation_fractions": [
                finite_double(value, "validation fraction")
                for value in VALIDATION_FRACTIONS
            ],
            "max_bisection_depth": MAX_DEPTH,
        },
        "constants": {
            "alpha_inverse": finite_double(ALPHA_INVERSE, "alpha inverse"),
            "alpha": finite_double(ALPHA, "alpha"),
            "hbarc_MeV_fm": finite_double(HBARC_MEV_FM, "hbar c"),
            "alpha_mass_MeV_c2": finite_double(
                ALPHA_MASS_MEV_C2, "alpha mass"
            ),
            "wronskian_tolerance": finite_double(
                WRONSKIAN_TOL, "Wronskian tolerance"
            ),
        },
        "elapsed_seconds": elapsed_seconds,
        "channels": channel_payload,
    }
    atomic_json_write(path, payload)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    repository = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(5, os.cpu_count() or 1),
        help="worker processes (default: min(5, CPU count))",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("/tmp/pb11-coulomb-table-cache"),
        help="per-energy-segment resume cache",
    )
    parser.add_argument(
        "--output-inc",
        type=Path,
        default=repository / "src/fusion_coulomb_tables.inc",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=repository / "docs/validation/coulomb_table_validation.json",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1 or args.jobs > 5:
        raise SystemExit("--jobs must be between 1 and 5")
    mp.mp.dps = MP_DPS
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    boundaries = energy_boundaries()
    tasks: list[tuple[str, int, str, str, str]] = []
    for family in FAMILIES:
        for index, (elo, ehi) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            xlo = mp_string(mp.log(elo))
            xhi = mp_string(mp.log(ehi))
            cache = args.cache_dir / f"{family}_segment_{index:03d}.json"
            tasks.append((family, index, xlo, xhi, str(cache)))

    started = time.monotonic()
    print(
        f"generating {len(tasks)} family segments with {args.jobs} workers "
        f"(cache={args.cache_dir})",
        flush=True,
    )
    completed: list[tuple[str, int, dict[str, Any]]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(run_segment_task, task): task for task in tasks
        }
        try:
            for future in concurrent.futures.as_completed(futures):
                task = futures[future]
                result = future.result()
                completed.append((task[0], task[1], result))
                print(
                    f"[{len(completed)}/{len(tasks)}] {task[0]} "
                    f"segment {task[1]} -> {len(result['segments'])} "
                    f"table segments; elapsed {time.monotonic() - started:.1f}s",
                    flush=True,
                )
        except BaseException:
            for future in futures:
                future.cancel()
            raise

    channel_segments: dict[int, list[dict[str, Any]]] = {
        int(channel["channel"]): [] for channel in CHANNELS
    }
    for family, index, payload in sorted(completed, key=lambda item: (item[0], item[1])):
        del index
        for segment in payload["segments"]:
            for channel in CHANNELS:
                if channel["family"] != family:
                    continue
                channel_segments[int(channel["channel"])].append(segment)
    for channel in channel_segments:
        channel_segments[channel].sort(key=lambda segment: mp.mpf(segment["lo"]))
        if not channel_segments[channel]:
            raise ReferenceFailure(f"no generated segments for channel {channel}")

    write_include(args.output_inc, channel_segments)
    write_validation_json(
        args.output_json, channel_segments, time.monotonic() - started
    )
    print(f"wrote {args.output_inc}", flush=True)
    print(f"wrote {args.output_json}", flush=True)
    print(
        "segment counts: "
        + ", ".join(
            f"ch{channel}={len(channel_segments[channel])}"
            for channel in sorted(channel_segments)
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReferenceFailure as error:
        print(f"reference generation failed: {error}", flush=True)
        raise SystemExit(2)
