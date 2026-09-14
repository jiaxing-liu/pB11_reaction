#!/usr/bin/env python3
"""Archive and plot the completed birth-table validation studies.

The inputs are a frozen allow-list of completed microbenchmarks and coupled
thermal runs.  CSV and log files are copied byte-for-byte into the study
directory and their SHA-256 digests are recorded before any numeric analysis.
The direct DT CSV is required to contain all 64 accepted steps; an incomplete
run raises an error rather than being treated as a completed reference.

The plotting environment used for the published artifact is::

    /home/cloud/research/pB-baldur/.venv-plots/bin/python \
        tools/plot_birth_table_validation.py

Use ``--input-dir`` to select another directory containing the same frozen
input names, or ``--output-dir`` to write a separate artifact directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import statistics
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = Path("/tmp")
DEFAULT_OUTPUT_DIR = REPO_ROOT / "docs" / "validation" / "birth-table" / "studies"

BIRTH_ERROR_FIELDS = (
    "rate_error",
    "debit_error",
    "max_number_L1",
    "max_energy_L1",
)
HEAT_FIELDS = {
    "electron": "electron_heat_J_m3",
    "network_ion": "network_ion_heat_J_m3",
    "carbon": "inert_heat_J_m3",
}
COUPLED_FIELDS = (
    "Te_keV",
    "Ti_keV",
    "Ue_J_m3",
    "Ui_J_m3",
    "fast_energy_J_m3",
    "neutron_energy_J_m3",
    "Q_J_m3",
    "energy_relative_residual",
    "inert_heat_J_m3",
    "events_pb",
    "events_dd_tp",
    "events_dd_he3n",
    "events_dt",
    "events_dhe3",
    "thermal_H",
    "thermal_D",
    "thermal_T",
    "thermal_He3",
    "thermal_He4",
    "thermal_B11",
    "fast_He4",
    "thermalized_He4",
    "projections",
    "electron_heat_J_m3",
    "network_ion_heat_J_m3",
    "removed_thermal_energy_J_m3",
)
TERMINAL_PLOT_FIELDS = (
    "Te_keV",
    "Ti_keV",
    "Ue_J_m3",
    "Ui_J_m3",
    "fast_energy_J_m3",
    "Q_J_m3",
)

BIRTH_SPECS = {
    "dt": "birth-table-dt",
    "pb": "birth-table-pb",
}
COUPLED_SPECS = {
    "dt": {
        "lookup": "coupled-dt-s64-n800-table001",
        "direct": "coupled-dt-s64-n800-table-era-direct",
    },
    "pb": {
        "lookup": "coupled-pb-s64-n800-table001",
        "direct": "coupled-pb-s64-n800-q8",
    },
}

BIRTH_LOG_RE = re.compile(
    r"channel=(?P<channel>\d+)\s+cells=(?P<cells>\d+)\s+"
    r"range_keV=(?P<lower>[-+0-9.eE]+):(?P<upper>[-+0-9.eE]+)\s+"
    r"tolerance=(?P<tolerance>[-+0-9.eE]+)\s+knots=(?P<knots>\d+)\s+"
    r"direct_evaluations=(?P<direct_evaluations>\d+)\s+"
    r"build_seconds=(?P<build_seconds>[-+0-9.eE]+)\s+"
    r"accepted_max_rate=(?P<accepted_max_rate>[-+0-9.eE]+)\s+"
    r"debit=(?P<accepted_max_debit>[-+0-9.eE]+)\s+"
    r"number_L1=(?P<accepted_max_number_L1>[-+0-9.eE]+)\s+"
    r"energy_L1=(?P<accepted_max_energy_L1>[-+0-9.eE]+)"
)
TABLE_LOG_RE = re.compile(
    r"table channel=(?P<channel>\d+)\s+knots=(?P<knots>\d+)\s+"
    r"direct_evaluations=(?P<direct_evaluations>\d+)\s+"
    r"build_seconds=(?P<build_seconds>[-+0-9.eE]+)\s+"
    r"max_sampled_rate=(?P<max_sampled_rate>[-+0-9.eE]+)\s+"
    r"max_sampled_debit=(?P<max_sampled_debit>[-+0-9.eE]+)\s+"
    r"max_sampled_number_L1=(?P<max_sampled_number_L1>[-+0-9.eE]+)\s+"
    r"max_sampled_energy_L1=(?P<max_sampled_energy_L1>[-+0-9.eE]+)"
)
COUPLED_LOG_RE = re.compile(
    r"fuel=(?P<fuel>\w+)\s+steps=(?P<steps>\d+)\s+cells=(?P<cells>\d+)\s+"
    r"duration=(?P<duration>[-+0-9.eE]+)\s+"
    r"max_total_energy_relative_residual=(?P<residual>[-+0-9.eE]+)\s+"
    r"projections=(?P<projections>\d+)\s+"
    r"neutron_number=(?P<neutron_number>[-+0-9.eE]+)"
)


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"required {label} is missing: {path}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_exact(source: Path, destination: Path, label: str) -> dict[str, Any]:
    require_file(source, label)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    source_hash = sha256_file(source)
    destination_hash = sha256_file(destination)
    if source_hash != destination_hash or source.stat().st_size != destination.stat().st_size:
        raise IOError(f"byte-for-byte copy failed for {label}: {source}")
    return {
        "file": destination.relative_to(DEFAULT_OUTPUT_DIR).as_posix()
        if destination.is_relative_to(DEFAULT_OUTPUT_DIR)
        else destination.name,
        "source": str(source),
        "bytes": destination.stat().st_size,
        "sha256": destination_hash,
    }


def finite(value: float) -> bool:
    return math.isfinite(value)


def load_birth_csv(path: Path) -> list[dict[str, float]]:
    require_file(path, "birth-table CSV")
    rows: list[dict[str, float]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no CSV header")
        required = {"T_keV", *BIRTH_ERROR_FIELDS, "direct_seconds", "table_seconds"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        for line_number, raw in enumerate(reader, start=2):
            if not raw or all(value in (None, "") for value in raw.values()):
                continue
            row: dict[str, float] = {}
            for field in reader.fieldnames:
                try:
                    row[field] = float(raw[field] or "")
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"invalid numeric value in {path}:{line_number}") from exc
                if not finite(row[field]):
                    raise ValueError(f"non-finite value in {path}:{line_number}:{field}")
            rows.append(row)
    if len(rows) != 4:
        raise ValueError(f"{path}: expected four off-construction samples, got {len(rows)}")
    temperatures = [row["T_keV"] for row in rows]
    if any(temperature <= 0 for temperature in temperatures) or any(
        b <= a for a, b in zip(temperatures, temperatures[1:])
    ):
        raise ValueError(f"{path}: sample temperatures are not strictly increasing")
    if any(row["direct_seconds"] <= 0 or row["table_seconds"] <= 0 for row in rows):
        raise ValueError(f"{path}: timing samples must be positive")
    return rows


def parse_birth_log(path: Path) -> dict[str, Any]:
    require_file(path, "birth-table log")
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    match = next((BIRTH_LOG_RE.fullmatch(line) for line in lines if line.startswith("channel=")), None)
    if match is None:
        raise ValueError(f"unrecognized birth-table log format: {path}")
    if "off-construction checks passed" not in lines:
        raise ValueError(f"birth-table log does not report completed checks: {path}")
    values: dict[str, Any] = {
        "channel": int(match["channel"]),
        "cells": int(match["cells"]),
        "range_keV": [float(match["lower"]), float(match["upper"])],
        "tolerance": float(match["tolerance"]),
        "knots": int(match["knots"]),
        "direct_evaluations": int(match["direct_evaluations"]),
        "build_seconds": float(match["build_seconds"]),
        "accepted_max_rate": float(match["accepted_max_rate"]),
        "accepted_max_debit": float(match["accepted_max_debit"]),
        "accepted_max_number_L1": float(match["accepted_max_number_L1"]),
        "accepted_max_energy_L1": float(match["accepted_max_energy_L1"]),
        "raw": path.read_text().strip(),
    }
    for key, value in values.items():
        if key == "raw":
            continue
        numeric_values = value if isinstance(value, list) else [value]
        if not all(finite(float(item)) for item in numeric_values):
            raise ValueError(f"non-finite value in birth-table log: {path}")
    return values


def load_coupled_csv(path: Path) -> list[dict[str, float | int]]:
    require_file(path, "coupled CSV")
    rows: list[dict[str, float | int]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no CSV header")
        required = {"step", "time_s", *COUPLED_FIELDS}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        for line_number, raw in enumerate(reader, start=2):
            if not raw or all(value in (None, "") for value in raw.values()):
                continue
            row: dict[str, float | int] = {}
            try:
                row["step"] = int(raw["step"] or "")
                for field in reader.fieldnames:
                    if field == "step":
                        continue
                    row[field] = float(raw[field] or "")
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid numeric value in {path}:{line_number}") from exc
            for field, value in row.items():
                if field != "step" and not finite(float(value)):
                    raise ValueError(f"non-finite value in {path}:{line_number}:{field}")
            rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no accepted coupled rows")
    steps = [int(row["step"]) for row in rows]
    if steps != list(range(1, len(rows) + 1)):
        raise ValueError(f"{path}: step sequence is not contiguous from one")
    times = [float(row["time_s"]) for row in rows]
    if any(time <= 0 for time in times) or any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError(f"{path}: accepted times are not strictly increasing")
    return rows


def parse_coupled_log(path: Path) -> dict[str, Any]:
    require_file(path, "coupled log")
    raw = path.read_text().strip()
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    table_match = next((TABLE_LOG_RE.fullmatch(line) for line in lines if line.startswith("table channel=")), None)
    completion_match = next((COUPLED_LOG_RE.fullmatch(line) for line in lines if line.startswith("fuel=")), None)
    if completion_match is None:
        raise ValueError(f"unrecognized coupled completion log format: {path}")
    completion: dict[str, Any] = {
        "fuel": completion_match["fuel"],
        "steps": int(completion_match["steps"]),
        "cells": int(completion_match["cells"]),
        "duration_s": float(completion_match["duration"]),
        "max_total_energy_relative_residual": float(completion_match["residual"]),
        "projections": int(completion_match["projections"]),
        "neutron_number_m3": float(completion_match["neutron_number"]),
        "raw": raw,
    }
    if table_match is None:
        table: dict[str, Any] | None = None
    else:
        table = {
            "channel": int(table_match["channel"]),
            "knots": int(table_match["knots"]),
            "direct_evaluations": int(table_match["direct_evaluations"]),
            "build_seconds": float(table_match["build_seconds"]),
            "max_sampled_rate": float(table_match["max_sampled_rate"]),
            "max_sampled_debit": float(table_match["max_sampled_debit"]),
            "max_sampled_number_L1": float(table_match["max_sampled_number_L1"]),
            "max_sampled_energy_L1": float(table_match["max_sampled_energy_L1"]),
        }
    numeric_values = [value for value in completion.values() if isinstance(value, (int, float))]
    if table is not None:
        numeric_values.extend(value for value in table.values() if isinstance(value, (int, float)))
    if not all(finite(float(value)) for value in numeric_values):
        raise ValueError(f"non-finite value in coupled log: {path}")
    return {"completion": completion, "table": table}


def archive_input_pair(
    input_dir: Path,
    output_dir: Path,
    stem: str,
    label: str,
) -> dict[str, Any]:
    csv_source = input_dir / f"{stem}.csv"
    log_source = input_dir / f"{stem}.log"
    csv_artifact = copy_exact(csv_source, output_dir / csv_source.name, f"{label} CSV")
    log_artifact = copy_exact(log_source, output_dir / log_source.name, f"{label} log")
    return {"csv": csv_artifact, "log": log_artifact}


def relative_error(table_value: float, direct_value: float) -> float | None:
    difference = abs(table_value - direct_value)
    denominator = abs(direct_value)
    if denominator == 0.0:
        return 0.0 if difference == 0.0 else None
    return difference / denominator


def compare_field(
    table_rows: list[dict[str, float | int]],
    direct_rows: list[dict[str, float | int]],
    field: str,
) -> dict[str, Any]:
    per_step: list[dict[str, float | int | None]] = []
    for table_row, direct_row in zip(table_rows, direct_rows):
        table_value = float(table_row[field])
        direct_value = float(direct_row[field])
        per_step.append(
            {
                "step": int(direct_row["step"]),
                "time_s": float(direct_row["time_s"]),
                "table": table_value,
                "direct": direct_value,
                "absolute_error": abs(table_value - direct_value),
                "relative_error": relative_error(table_value, direct_value),
            }
        )
    relative_values = [
        float(entry["relative_error"])
        for entry in per_step
        if entry["relative_error"] is not None
    ]
    terminal = per_step[-1]
    return {
        "max_absolute_error": max(float(entry["absolute_error"]) for entry in per_step),
        "max_relative_error": max(relative_values) if relative_values else None,
        "terminal": terminal,
        "per_step": per_step,
    }


def compare_coupled(
    table_rows: list[dict[str, float | int]],
    direct_rows: list[dict[str, float | int]],
) -> dict[str, Any]:
    if len(table_rows) != 64 or len(direct_rows) != 64:
        raise ValueError(
            "coupled comparison requires 64 accepted rows in both lookup and direct CSVs"
        )
    for table_row, direct_row in zip(table_rows, direct_rows):
        if int(table_row["step"]) != int(direct_row["step"]):
            raise ValueError("lookup/direct step sequences differ")
        table_time = float(table_row["time_s"])
        direct_time = float(direct_row["time_s"])
        if abs(table_time - direct_time) > 1.0e-12 * max(abs(table_time), abs(direct_time), 1.0):
            raise ValueError("lookup/direct time grids differ")
    return {
        "steps": len(table_rows),
        "fields": {
            field: compare_field(table_rows, direct_rows, field) for field in COUPLED_FIELDS
        },
    }


def birth_summary(
    rows_by_fuel: dict[str, list[dict[str, float]]],
    logs_by_fuel: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for fuel, rows in rows_by_fuel.items():
        log = logs_by_fuel[fuel]
        sample_records: list[dict[str, float]] = []
        for row in rows:
            maximum = max(row[field] for field in BIRTH_ERROR_FIELDS)
            sample_records.append(
                {
                    "T_keV": row["T_keV"],
                    **{field: row[field] for field in BIRTH_ERROR_FIELDS},
                    "max_error": maximum,
                    "direct_seconds": row["direct_seconds"],
                    "table_seconds": row["table_seconds"],
                    "direct_lookup_ratio": row["direct_seconds"] / row["table_seconds"],
                }
            )
        max_errors = {field: max(row[field] for row in rows) for field in BIRTH_ERROR_FIELDS}
        log_max_errors = {
            "rate_error": log["accepted_max_rate"],
            "debit_error": log["accepted_max_debit"],
            "max_number_L1": log["accepted_max_number_L1"],
            "max_energy_L1": log["accepted_max_energy_L1"],
        }
        summary[fuel] = {
            "samples": sample_records,
            "sample_count": len(rows),
            "temperature_range_keV": [rows[0]["T_keV"], rows[-1]["T_keV"]],
            "off_construction_max_errors": max_errors,
            "off_construction_global_max_error": max(max_errors.values()),
            "off_construction_gate": 1.0e-3,
            "off_construction_gate_pass": max(max_errors.values()) <= 1.0e-3,
            "log_accepted_max_errors": log_max_errors,
            "timing_ratio": {
                "samples": [record["direct_lookup_ratio"] for record in sample_records],
                "min": min(record["direct_lookup_ratio"] for record in sample_records),
                "max": max(record["direct_lookup_ratio"] for record in sample_records),
                "median": statistics.median(
                    record["direct_lookup_ratio"] for record in sample_records
                ),
                "direct_seconds_total": sum(record["direct_seconds"] for record in sample_records),
                "table_seconds_total": sum(record["table_seconds"] for record in sample_records),
                "birth_table_build_seconds": log["build_seconds"],
            },
            "construction_log": log,
        }
    return summary


def add_terminal_aliases(comparison: dict[str, Any]) -> dict[str, Any]:
    heat: dict[str, Any] = {}
    for name, field in HEAT_FIELDS.items():
        heat[name] = comparison["fields"][field]
    comparison["independent_heat_errors"] = heat
    return comparison


def binary_record(path: Path, label: str) -> dict[str, str]:
    require_file(path, label)
    return {"path": str(path), "sha256": sha256_file(path)}


def write_sha256sums(output_dir: Path, artifact_paths: list[Path]) -> None:
    lines = []
    for path in sorted(artifact_paths):
        lines.append(f"{sha256_file(path)}  {path.relative_to(output_dir).as_posix()}")
    (output_dir / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")


def plot_validation(
    output_dir: Path,
    birth_rows: dict[str, list[dict[str, float]]],
    birth_summary_data: dict[str, Any],
    comparisons: dict[str, dict[str, Any]],
) -> None:
    fuel_colors = {"dt": "#1769aa", "pb": "#b24745"}
    heat_colors = {"electron": "#1769aa", "network_ion": "#2f8f5b", "carbon": "#b24745"}
    fuel_labels = {"dt": "DT", "pb": "pB"}
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 9.0), constrained_layout=True)
    figure.suptitle("Birth-table validation and coupled lookup parity", fontsize=15)

    # Panel 1: the maximum of the four independent off-construction metrics at
    # each sampled temperature, with the declared 0.1% gate.
    axis = axes[0, 0]
    for fuel in ("dt", "pb"):
        rows = birth_rows[fuel]
        axis.plot(
            list(range(1, len(rows) + 1)),
            [max(row[field] for field in BIRTH_ERROR_FIELDS) for row in rows],
            marker="o",
            linewidth=1.8,
            color=fuel_colors[fuel],
            label=(
                f"{fuel_labels[fuel]} ({rows[0]['T_keV']:.4g}-{rows[-1]['T_keV']:.4g} keV)"
            ),
        )
    axis.axhline(1.0e-3, color="#444444", linestyle="--", linewidth=1.0, label="0.1% reference")
    axis.set_yscale("log")
    axis.set_xlabel("off-construction sample index")
    axis.set_xticks([1, 2, 3, 4])
    axis.set_ylabel("max off-construction relative error")
    axis.set_title("Off-construction interpolation error")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=8)

    # Panel 2: direct evaluation versus lookup, with construction costs kept
    # in a separate annotation so they are not mistaken for per-call timing.
    axis = axes[0, 1]
    for fuel in ("dt", "pb"):
        samples = birth_summary_data[fuel]["samples"]
        axis.plot(
            list(range(1, len(samples) + 1)),
            [sample["direct_lookup_ratio"] for sample in samples],
            marker="o",
            linewidth=1.8,
            color=fuel_colors[fuel],
            label=(
                f"{fuel_labels[fuel]} ({samples[0]['T_keV']:.4g}-{samples[-1]['T_keV']:.4g} keV)"
            ),
        )
    axis.set_yscale("log")
    axis.set_xlabel("off-construction sample index")
    axis.set_xticks([1, 2, 3, 4])
    axis.set_ylabel("direct time / lookup time")
    axis.set_title("Lookup timing ratio")
    costs = "; ".join(
        f"{fuel_labels[fuel]} table build {birth_summary_data[fuel]['timing_ratio']['birth_table_build_seconds']:.3g} s"
        for fuel in ("dt", "pb")
    )
    axis.text(
        0.97,
        0.50,
        costs,
        transform=axis.transAxes,
        fontsize=8,
        ha="right",
        va="center",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "#bbbbbb"},
    )
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=8)

    # Panel 3: three separate cumulative heat ledgers.  Solid lines are DT,
    # dashed lines are pB, and color identifies the heat destination.
    axis = axes[1, 0]
    heat_labels = {"electron": "electron", "network_ion": "ions", "carbon": "carbon"}
    for fuel, linestyle in (("dt", "-"), ("pb", "--")):
        comparison = comparisons[fuel]
        for name, field in HEAT_FIELDS.items():
            entries = comparison["fields"][field]["per_step"]
            values = [
                max(float(entry["relative_error"] or 0.0) * 100.0, 1.0e-14)
                for entry in entries
            ]
            axis.plot(
                [float(entry["time_s"]) * 1.0e3 for entry in entries],
                values,
                linestyle=linestyle,
                linewidth=1.5,
                color=heat_colors[name],
                label=f"{fuel_labels[fuel]} {heat_labels[name]}",
            )
    axis.axhline(0.1, color="#444444", linestyle=":", linewidth=1.0, label="0.1% reference")
    axis.set_yscale("log")
    axis.set_xlabel("time (ms)")
    axis.set_ylabel("absolute relative error (%)")
    axis.set_title("Cumulative heat parity: lookup vs direct")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=7, ncol=2)

    # Panel 4: terminal parity for a compact set of nonzero physical fields.
    axis = axes[1, 1]
    positions = list(range(len(TERMINAL_PLOT_FIELDS)))
    width = 0.36
    dt_values = []
    pb_values = []
    for field in TERMINAL_PLOT_FIELDS:
        dt_error = comparisons["dt"]["fields"][field]["terminal"]["relative_error"]
        pb_error = comparisons["pb"]["fields"][field]["terminal"]["relative_error"]
        if dt_error is None or pb_error is None or dt_error <= 0.0 or pb_error <= 0.0:
            raise ValueError(f"terminal plot field {field} must have nonzero finite parity errors")
        dt_values.append(float(dt_error) * 100.0)
        pb_values.append(float(pb_error) * 100.0)
    axis.bar([position - width / 2 for position in positions], dt_values, width, label="DT", color=fuel_colors["dt"])
    axis.bar([position + width / 2 for position in positions], pb_values, width, label="pB", color=fuel_colors["pb"])
    axis.set_yscale("log")
    axis.set_xticks(positions)
    axis.set_xticklabels(
        [
            "Te",
            "Ti",
            "Ue",
            "Ui",
            "fast E",
            "Q",
        ]
    )
    axis.set_ylabel("absolute terminal relative error (%)")
    axis.set_title("Terminal physical-quantity parity")
    axis.grid(True, which="both", axis="y", alpha=0.25)
    axis.legend(fontsize=8)

    png_path = output_dir / "birth-table-validation.png"
    pdf_path = output_dir / "birth-table-validation.pdf"
    figure.savefig(png_path, dpi=180)
    figure.savefig(pdf_path)
    plt.close(figure)


def render_readme(
    birth_data: dict[str, Any],
    coupled_data: dict[str, Any],
    binary_data: dict[str, Any],
) -> str:
    fuel_display = {"dt": "DT", "pb": "pB"}
    lines = [
        "# Birth-table validation studies",
        "",
        "This directory archives the completed numerical acceleration checks.",
        "CSV and log artifacts are copied byte-for-byte from `/tmp`; their SHA-256",
        "digests are in `summary.json` and `SHA256SUMS.txt`.",
        "",
        "## Inputs and scope",
        "",
        "- The birth-table microbenchmarks use four off-construction samples for DT",
        "  and pB. The declared acceptance gate is 0.1% for each rate/debit/number/energy",
        "  metric.",
        "- The coupled comparisons retain all 64 accepted time points for each fuel and",
        "  compare the table run with its direct-evaluation reference.",
        "- `electron_heat_J_m3`, `network_ion_heat_J_m3`, and `inert_heat_J_m3` are",
        "  reported independently; the last is the carbon-bath heat ledger.",
        "- This is a numerical acceleration and parity study. It is not an EXL study",
        "  or a fair performance ranking of fuels.",
        "",
        "## Build provenance",
        "",
        f"- Birth-table study binary: `{binary_data['birth_table']['path']}`",
        f"  (SHA-256 `{binary_data['birth_table']['sha256']}`).",
        f"- Coupled driver binary: `{binary_data['coupled_driver']['path']}`",
        f"  (SHA-256 `{binary_data['coupled_driver']['sha256']}`).",
        "- The microbenchmark binary predates the two sampled-direct diagnostic fields;",
        "  the successful interpolation values are unchanged. No source hash is used",
        "  as compile provenance.",
        "",
        "## Results",
        "",
    ]
    for fuel in ("dt", "pb"):
        birth = birth_data[fuel]
        timing = birth["timing_ratio"]
        lines.append(
            f"- **{fuel_display[fuel]}**: off-construction max "
            f"`{birth['off_construction_global_max_error']:.6g}`, "
            f"direct/lookup timing ratio "
            f"`{timing['min']:.4g}`–`{timing['max']:.4g}`, "
            f"table construction `{timing['birth_table_build_seconds']:.6g} s`."
        )
        comparison = coupled_data[fuel]
        heat_text = ", ".join(
            f"{name} {100.0 * float(comparison['independent_heat_errors'][name]['terminal']['relative_error'] or 0.0):.6g}%"
            for name in HEAT_FIELDS
        )
        lines.append(f"  Terminal independent heat relative errors: {heat_text}.")
    lines.extend(
        [
            "",
            "The four-panel figure is `birth-table-validation.png` (also available as",
            "`birth-table-validation.pdf`). Panel 2 uses a logarithmic timing axis and",
            "annotates construction cost separately; panels 3 and 4 show parity errors",
            "against the direct coupled runs.",
        ]
    )
    return "\n".join(lines) + "\n"


def build_artifact(args: argparse.Namespace) -> None:
    input_dir = args.input_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    birth_artifacts: dict[str, Any] = {}
    birth_rows: dict[str, list[dict[str, float]]] = {}
    birth_logs: dict[str, dict[str, Any]] = {}
    for fuel, stem in BIRTH_SPECS.items():
        birth_artifacts[fuel] = archive_input_pair(input_dir, output_dir, stem, f"{fuel} birth-table")
        birth_rows[fuel] = load_birth_csv(input_dir / f"{stem}.csv")
        birth_logs[fuel] = parse_birth_log(input_dir / f"{stem}.log")

    coupled_artifacts: dict[str, Any] = {}
    coupled_rows: dict[str, dict[str, list[dict[str, float | int]]]] = {}
    coupled_logs: dict[str, dict[str, Any]] = {}
    coupled_comparisons: dict[str, dict[str, Any]] = {}
    for fuel, spec in COUPLED_SPECS.items():
        coupled_artifacts[fuel] = {}
        coupled_rows[fuel] = {}
        coupled_logs[fuel] = {}
        for role, stem in spec.items():
            label = f"{fuel} coupled {role}"
            coupled_artifacts[fuel][role] = archive_input_pair(input_dir, output_dir, stem, label)
            rows = load_coupled_csv(input_dir / f"{stem}.csv")
            log = parse_coupled_log(input_dir / f"{stem}.log")
            if log["completion"]["steps"] != 64 or len(rows) != 64:
                raise ValueError(
                    f"{label} is incomplete: log steps={log['completion']['steps']}, rows={len(rows)}"
                )
            if log["completion"]["fuel"] != fuel:
                raise ValueError(f"{label} fuel mismatch in completion log")
            coupled_rows[fuel][role] = rows
            coupled_logs[fuel][role] = log
        coupled_comparisons[fuel] = add_terminal_aliases(
            compare_coupled(coupled_rows[fuel]["lookup"], coupled_rows[fuel]["direct"])
        )

    birth_data = birth_summary(birth_rows, birth_logs)
    coupled_data: dict[str, Any] = {}
    for fuel in ("dt", "pb"):
        coupled_data[fuel] = {
            "lookup_log": coupled_logs[fuel]["lookup"],
            "direct_log": coupled_logs[fuel]["direct"],
            "comparison": coupled_comparisons[fuel],
            "independent_heat_errors": coupled_comparisons[fuel]["independent_heat_errors"],
        }

    binary_data = {
        "birth_table": binary_record(input_dir / "study_birth_table", "birth-table study binary"),
        "coupled_driver": binary_record(
            input_dir / "study_coupled_thermal_table", "coupled thermal study binary"
        ),
    }
    summary: dict[str, Any] = {
        "schema_version": 1,
        "analysis": {
            "script": str(Path(__file__).resolve()),
            "input_dir": str(input_dir),
            "output_dir": str(output_dir),
            "off_construction_gate": 1.0e-3,
            "coupled_steps": 64,
            "comparison_convention": "absolute difference and absolute difference divided by direct value",
        },
        "provenance": {
            "study_binaries": binary_data,
            "microbenchmark_note": (
                "The microbenchmark binary predates two sampled-direct diagnostic fields; "
                "successful interpolation values are unchanged."
            ),
            "source_hash_policy": "No source hash is used as compile provenance; source hashes, if computed elsewhere, are analysis-time only.",
        },
        "birth_table": {
            "artifacts": birth_artifacts,
            "fuels": birth_data,
        },
        "coupled": {
            "artifacts": coupled_artifacts,
            "fuels": coupled_data,
        },
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    archived_files = [
        path
        for path in output_dir.iterdir()
        if path.is_file() and path.suffix in {".csv", ".log"}
    ]
    write_sha256sums(output_dir, archived_files)
    plot_validation(output_dir, birth_rows, birth_data, coupled_comparisons)
    (output_dir / "README.md").write_text(render_readme(birth_data, coupled_data, binary_data))

    print(f"archived birth-table/coupled studies to {output_dir}")
    for fuel in ("dt", "pb"):
        print(
            f"{fuel}: off_max={birth_data[fuel]['off_construction_global_max_error']:.12g} "
            f"timing_ratio={birth_data[fuel]['timing_ratio']['min']:.6g}:"
            f"{birth_data[fuel]['timing_ratio']['max']:.6g} "
            f"terminal_heat_rel_percent="
            + ",".join(
                f"{name}={100.0 * float(coupled_data[fuel]['independent_heat_errors'][name]['terminal']['relative_error'] or 0.0):.12g}"
                for name in HEAT_FIELDS
            )
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


if __name__ == "__main__":
    build_artifact(parse_args())
