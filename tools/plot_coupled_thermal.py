#!/usr/bin/env python3
"""Analyze and plot the completed coupled-thermal laboratory studies.

The study driver writes one CSV row per accepted local trial and a one-line
completion/failure log.  This script deliberately has an allow-list of the
completed runs.  In particular, pending pB ``*-q*.csv`` files are not
collected by accident.  Running the script again with the same inputs
replaces the copied artifacts and regenerates the JSON and DT figure.

Example::

    MPLCONFIGDIR=/tmp/pb-mpl \
      .venv-plots/bin/python tools/plot_coupled_thermal.py

Use ``--input-dir`` to analyze an archived input directory and ``--output``
to select another artifact directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

# The managed environment's home is read-only.  The command line used for the
# published run sets this explicitly; this fallback keeps the utility usable
# when called directly as well.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pb-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


MEV_J = 1.602176634e-13
KEV_J = 1.602176634e-16
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = Path("/tmp")
DEFAULT_OUTPUT_DIR = REPO_ROOT / "docs" / "validation" / "coupled-thermal"


@dataclass(frozen=True)
class RunSpec:
    stem: str
    fuel: str
    steps: int
    cells: int
    duration_s: float
    initial_ti_keV: float
    initial_te_keV: float
    kind: str = "success"
    expected_status: int | None = None


# This is the frozen completed set for the validation handoff.  Do not add
# the still-running pB q-resolution files here; root owns those inputs.
RUN_SPECS: tuple[RunSpec, ...] = (
    RunSpec("coupled-dt-s16-n800", "dt", 16, 800, 0.01, 20.0, 5.0),
    RunSpec("coupled-dt-s32-n800", "dt", 32, 800, 0.01, 20.0, 5.0),
    RunSpec("coupled-dt-s64-n800", "dt", 64, 800, 0.01, 20.0, 5.0),
    RunSpec("coupled-dt-s128-n800", "dt", 128, 800, 0.01, 20.0, 5.0),
    RunSpec("coupled-dt-s128-n1600", "dt", 128, 1600, 0.01, 20.0, 5.0),
    RunSpec("coupled-dd-s16-n800", "dd", 16, 800, 0.01, 20.0, 5.0),
    RunSpec("coupled-dhe3-s16-n800", "dhe3", 16, 800, 0.01, 20.0, 5.0),
    RunSpec("coupled-pb-s4-n800", "pb", 4, 800, 0.001, 100.0, 5.0),
    RunSpec(
        "coupled-pb-s4-n800-first-grid-failure",
        "pb",
        4,
        800,
        0.001,
        100.0,
        5.0,
        kind="failure",
        expected_status=3,
    ),
)

# This successful rerun uses the guard-fixed library.  It is retained as a
# parity artifact, but is not a second point in the timestep convergence set.
PARITY_STEM = "coupled-dt-s32-n800-fixed-guard"
PARITY_REFERENCE_STEM = "coupled-dt-s32-n800"
PARITY_SPEC = RunSpec(PARITY_STEM, "dt", 32, 800, 0.01, 20.0, 5.0)

DT_STEMS = (
    "coupled-dt-s16-n800",
    "coupled-dt-s32-n800",
    "coupled-dt-s64-n800",
    "coupled-dt-s128-n800",
)
DT_GRID_STEMS = ("coupled-dt-s128-n800", "coupled-dt-s128-n1600")

FINAL_FIELDS = (
    "time_s",
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
)

SUCCESS_LOG_RE = re.compile(
    r"fuel=(?P<fuel>\w+)\s+steps=(?P<steps>\d+)\s+cells=(?P<cells>\d+)\s+"
    r"duration=(?P<duration>[-+0-9.eE]+)\s+"
    r"max_total_energy_relative_residual=(?P<residual>[-+0-9.eE]+)\s+"
    r"projections=(?P<projections>\d+)\s+"
    r"neutron_number=(?P<neutron_number>[-+0-9.eE]+)"
)
FAILURE_LOG_RE = re.compile(
    r"trial failed\s+fuel=(?P<fuel>\w+)\s+step=(?P<step>\d+)\s+"
    r"status=(?P<status>\d+)\s+Ti=(?P<ti>[-+0-9.eE]+)keV\s+"
    r"Te=(?P<te>[-+0-9.eE]+)keV"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: float) -> bool:
    return math.isfinite(value)


def load_csv(path: Path) -> list[dict[str, float | int]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no CSV header")
        required = {"step", "time_s", "projections"} | set(FINAL_FIELDS[1:])
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        rows: list[dict[str, float | int]] = []
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
    return rows


def parse_log(path: Path, kind: str) -> dict[str, Any]:
    text = path.read_text().strip()
    if kind == "failure":
        match = FAILURE_LOG_RE.fullmatch(text)
        if match is None:
            raise ValueError(f"unrecognized failure log format: {path}")
        return {
            "status": "failed_trial",
            "fuel": match["fuel"],
            "failed_step": int(match["step"]),
            "return_status": int(match["status"]),
            "Ti_keV_at_failure": float(match["ti"]),
            "Te_keV_at_failure": float(match["te"]),
            "raw": text,
        }
    match = SUCCESS_LOG_RE.fullmatch(text)
    if match is None:
        raise ValueError(f"unrecognized completion log format: {path}")
    return {
        "status": "completed",
        "fuel": match["fuel"],
        "steps": int(match["steps"]),
        "cells": int(match["cells"]),
        "duration_s": float(match["duration"]),
        "max_total_energy_relative_residual": float(match["residual"]),
        "total_projections": int(match["projections"]),
        "neutron_number_m3": float(match["neutron_number"]),
        "raw": text,
    }


def validate_rows(spec: RunSpec, rows: list[dict[str, float | int]]) -> None:
    if spec.kind == "success" and len(rows) != spec.steps:
        raise ValueError(
            f"{spec.stem}: expected {spec.steps} accepted rows, got {len(rows)}"
        )
    if spec.kind == "failure" and len(rows) != 1:
        raise ValueError(f"{spec.stem}: expected the one preserved pre-failure row")
    expected_steps = list(range(1, len(rows) + 1))
    actual_steps = [int(row["step"]) for row in rows]
    if actual_steps != expected_steps:
        raise ValueError(f"{spec.stem}: step sequence is {actual_steps!r}")
    times = [float(row["time_s"]) for row in rows]
    if any(t <= 0 for t in times) or any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError(f"{spec.stem}: accepted times are not strictly increasing")
    if any(int(row["projections"]) < 0 for row in rows):
        raise ValueError(f"{spec.stem}: negative handoff projection count")


def final_metrics(rows: list[dict[str, float | int]]) -> dict[str, float | int]:
    if not rows:
        raise ValueError("cannot summarize empty CSV")
    last = rows[-1]
    return {field: last[field] for field in FINAL_FIELDS}


def max_abs(rows: Iterable[dict[str, float | int]], field: str) -> float:
    return max(abs(float(row[field])) for row in rows)


def relative_change(from_value: float, to_value: float) -> float:
    denominator = abs(to_value)
    if denominator == 0.0:
        return 0.0 if from_value == 0.0 else float("inf")
    return abs(from_value - to_value) / denominator


def selected_metrics(
    from_metrics: dict[str, float | int], to_metrics: dict[str, float | int]
) -> dict[str, float]:
    fields = ("Ue_J_m3", "Ui_J_m3", "fast_energy_J_m3", "Q_J_m3")
    return {
        field: relative_change(float(from_metrics[field]), float(to_metrics[field]))
        for field in fields
    }


def copied_artifact(source: Path, destination: Path, kind: str) -> dict[str, Any]:
    if not source.is_file():
        raise FileNotFoundError(f"required {kind} is missing: {source}")
    shutil.copyfile(source, destination)
    if source.read_bytes() != destination.read_bytes():
        raise IOError(f"copied {kind} differs from source: {source}")
    return {
        "file": destination.name,
        "source_name": source.name,
        "sha256": sha256_file(destination),
        "bytes": destination.stat().st_size,
    }


def run_record(spec: RunSpec, input_dir: Path, output_dir: Path) -> dict[str, Any]:
    csv_source = input_dir / f"{spec.stem}.csv"
    log_source = input_dir / f"{spec.stem}.log"
    csv_destination = output_dir / csv_source.name
    log_destination = output_dir / log_source.name
    csv_artifact = copied_artifact(csv_source, csv_destination, "CSV")
    log_artifact = copied_artifact(log_source, log_destination, "log")
    rows = load_csv(csv_destination)
    validate_rows(spec, rows)
    log = parse_log(log_destination, spec.kind)
    if log["fuel"] != spec.fuel:
        raise ValueError(f"{spec.stem}: log fuel disagrees with run metadata")
    if spec.kind == "success":
        if (
            log["steps"] != spec.steps
            or log["cells"] != spec.cells
            or not math.isclose(log["duration_s"], spec.duration_s, rel_tol=0.0, abs_tol=1e-12)
        ):
            raise ValueError(f"{spec.stem}: completion log disagrees with run metadata")
        if int(log["total_projections"]) != sum(int(row["projections"]) for row in rows):
            raise ValueError(f"{spec.stem}: projection count disagrees with CSV")
        if not finite(float(log["max_total_energy_relative_residual"])):
            raise ValueError(f"{spec.stem}: invalid residual in completion log")
    else:
        if log["return_status"] != spec.expected_status:
            raise ValueError(f"{spec.stem}: preserved failure status changed")
        if log["failed_step"] != 2:
            raise ValueError(f"{spec.stem}: preserved failure step changed")
    return {
        "fuel": spec.fuel,
        "steps": spec.steps,
        "cells": spec.cells,
        "duration_s": spec.duration_s,
        "dt_s": spec.duration_s / spec.steps,
        "initial_Ti_keV": spec.initial_ti_keV,
        "initial_Te_keV": spec.initial_te_keV,
        "kind": spec.kind,
        "rows": len(rows),
        "completed_steps": len(rows),
        "csv": csv_artifact,
        "log_file": log_artifact,
        "log": log,
        "final": final_metrics(rows),
        "max_abs_energy_relative_residual_csv": max_abs(rows, "energy_relative_residual"),
        "total_projections_csv": sum(int(row["projections"]) for row in rows),
        "rows_data": rows,
    }


def compare_common_csv_values(
    reference: dict[str, Any], candidate: dict[str, Any]
) -> dict[str, Any]:
    """Compare the shared numeric columns of two successful CSV artifacts."""

    reference_rows = reference["rows_data"]
    candidate_rows = candidate["rows_data"]
    if len(reference_rows) != len(candidate_rows):
        raise ValueError(
            "guard-fix parity row counts differ: "
            f"{len(reference_rows)} versus {len(candidate_rows)}"
        )
    if not reference_rows or not candidate_rows:
        raise ValueError("guard-fix parity comparison needs non-empty CSV files")

    reference_columns = tuple(reference_rows[0].keys())
    candidate_columns = tuple(candidate_rows[0].keys())
    common_columns = tuple(field for field in reference_columns if field in candidate_columns)
    added_columns = tuple(field for field in candidate_columns if field not in reference_columns)
    if not common_columns:
        raise ValueError("guard-fix parity comparison has no common columns")
    if any(set(row) != set(reference_columns) for row in reference_rows):
        raise ValueError("reference parity CSV rows do not have a stable schema")
    if any(set(row) != set(candidate_columns) for row in candidate_rows):
        raise ValueError("guard-fix parity CSV rows do not have a stable schema")

    max_absolute_difference = 0.0
    max_relative_difference = 0.0
    exact = True
    for reference_row, candidate_row in zip(reference_rows, candidate_rows):
        for field in common_columns:
            reference_value = float(reference_row[field])
            candidate_value = float(candidate_row[field])
            difference = abs(reference_value - candidate_value)
            max_absolute_difference = max(max_absolute_difference, difference)
            max_relative_difference = max(
                max_relative_difference,
                difference / max(abs(reference_value), abs(candidate_value), 1.0e-300),
            )
            exact = exact and reference_value == candidate_value

    return {
        "reference_csv": reference["csv"]["file"],
        "guard_fix_csv": candidate["csv"]["file"],
        "rows_compared": len(reference_rows),
        "common_columns": list(common_columns),
        "guard_fix_added_columns": list(added_columns),
        "all_common_numeric_values_exact": exact,
        "max_absolute_difference": max_absolute_difference,
        "max_relative_difference": max_relative_difference,
        "logs_byte_identical": reference["log_file"]["sha256"] == candidate["log_file"]["sha256"],
        "interpretation": (
            "The original continuous study was linked against the pre-guard-fix library. "
            "The guard-fixed rerun has identical shared CSV values, so the successful "
            "trajectory numerical operator is unchanged; its added columns are retained "
            "only as guard-fix run output."
        ),
    }


def make_plot(records: dict[str, dict[str, Any]], output_dir: Path) -> list[str]:
    dt = records["coupled-dt-s128-n800"]["rows_data"]
    time_ms = np.asarray([float(row["time_s"]) * 1e3 for row in dt])
    te = np.asarray([float(row["Te_keV"]) for row in dt])
    ti = np.asarray([float(row["Ti_keV"]) for row in dt])
    u_e = np.asarray([float(row["Ue_J_m3"]) for row in dt])
    u_i = np.asarray([float(row["Ui_J_m3"]) for row in dt])
    fast = np.asarray([float(row["fast_energy_J_m3"]) for row in dt])
    neutron = np.asarray([float(row["neutron_energy_J_m3"]) for row in dt])
    q = np.asarray([float(row["Q_J_m3"]) for row in dt])
    carbon_heat = np.asarray([float(row["inert_heat_J_m3"]) for row in dt])

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "figure.dpi": 120,
            "savefig.dpi": 240,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0))
    # Reserve a real footer band for the scope statement.  A fixed layout is
    # more reproducible for the paired PNG/PDF than letting tight layout place
    # the footer on top of the lower x-axis labels.
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.12, top=0.88, wspace=0.24, hspace=0.30)
    ax_temp, ax_energy, ax_step, ax_carbon = axes.flat

    ax_temp.plot(time_ms, te, color="#0072B2", lw=2.0, label=r"$T_e$")
    ax_temp.plot(time_ms, ti, color="#D55E00", lw=2.0, label=r"$T_i$")
    ax_temp.set_title("(a) Evolving thermal temperatures")
    ax_temp.set_xlabel("time (ms)")
    ax_temp.set_ylabel("thermal energy temperature (keV)")
    ax_temp.grid(alpha=0.25)
    ax_temp.legend(loc="best", frameon=False)

    energy_series = (
        (u_e, r"$U_e$", "#0072B2"),
        (u_i, r"$U_i$", "#D55E00"),
        (fast, r"$U_{fast}$", "#009E73"),
        (neutron, r"$E_{n,\,cum}$", "#CC79A7"),
        (q, r"$Q_{cum}$", "#000000"),
    )
    for series, label, color in energy_series:
        ax_energy.plot(time_ms, series, lw=1.8, color=color, label=label)
    ax_energy.set_title("(b) Local energy accounts")
    ax_energy.set_xlabel("time (ms)")
    ax_energy.set_ylabel(r"energy density (J m$^{-3}$)")
    ax_energy.grid(alpha=0.25)
    ax_energy.legend(loc="upper left", frameon=False, ncol=2)
    ax_energy.text(
        0.02,
        0.04,
        "$U_e$, $U_i$, $U_{fast}$: instantaneous\n$Q$ and neutron energy: cumulative",
        transform=ax_energy.transAxes,
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.8},
    )

    timestep_records = [records[stem] for stem in DT_STEMS]
    reference = timestep_records[-1]["final"]
    steps = np.asarray([record["steps"] for record in timestep_records], dtype=float)
    convergence_series = (
        ("Ue_J_m3", r"$U_e$", "#0072B2"),
        ("Ui_J_m3", r"$U_i$", "#D55E00"),
        ("fast_energy_J_m3", r"$U_{fast}$", "#009E73"),
        ("Q_J_m3", r"$Q$", "#000000"),
    )
    for field, label, color in convergence_series:
        values = np.asarray(
            [
                100.0
                * relative_change(float(record["final"][field]), float(reference[field]))
                for record in timestep_records
            ]
        )
        ax_step.plot(steps, values, marker="o", lw=1.8, color=color, label=label)
    ax_step.set_title("(c) Timestep convergence at 800 cells")
    ax_step.set_xlabel("number of accepted steps over 10 ms")
    ax_step.set_ylabel("relative change from 128-step result (%)")
    ax_step.set_xticks(steps)
    ax_step.set_yscale("symlog", linthresh=1e-4)
    ax_step.set_ylim(bottom=0.0)
    ax_step.grid(alpha=0.25, which="both")
    ax_step.legend(loc="best", frameon=False)

    ax_carbon.plot(time_ms, carbon_heat, color="#CC79A7", lw=2.0)
    ax_carbon.set_title("(d) C-bath collision heat")
    ax_carbon.set_xlabel("time (ms)")
    ax_carbon.set_ylabel(r"cumulative carbon heat (J m$^{-3}$)")
    ax_carbon.grid(alpha=0.25)
    ax_carbon.text(
        0.03,
        0.93,
        r"$n_C=10^{19}$ m$^{-3}=10^{-3}n_e$; $\langle Z^2\rangle=36$",
        transform=ax_carbon.transAxes,
        va="top",
        fontsize=8,
    )

    fig.suptitle(
        "DT local coupled thermal trial: evolving bath and convergence",
        fontsize=14,
    )
    fig.text(
        0.5,
        0.025,
        "Numerical laboratory study; density $10^{22}$ m$^{-3}$, initial $T_i=20$ keV, $T_e=5$ keV. "
        "No EXL device prediction or fair four-fuel ranking.",
        ha="center",
        va="bottom",
        fontsize=8.5,
    )

    png_path = output_dir / "dt_coupled_thermal.png"
    pdf_path = output_dir / "dt_coupled_thermal.pdf"
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        metadata={
            "Creator": "plot_coupled_thermal.py",
            "Title": "DT local coupled thermal validation",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)
    return [png_path.name, pdf_path.name]


def build_summary(
    records: dict[str, dict[str, Any]],
    parity_record: dict[str, Any],
    parity_comparison: dict[str, Any],
    output_dir: Path,
    plot_files: list[str],
) -> dict[str, Any]:
    dt_records = [records[stem] for stem in DT_STEMS]
    dt_reference = records["coupled-dt-s128-n800"]
    dt_grid_coarse = records[DT_GRID_STEMS[0]]
    dt_grid_fine = records[DT_GRID_STEMS[1]]
    adjacent = []
    for previous, current in zip(dt_records, dt_records[1:]):
        adjacent.append(
            {
                "from_steps": previous["steps"],
                "to_steps": current["steps"],
                "relative_change_fraction": selected_metrics(
                    previous["final"], current["final"]
                ),
                "relative_change_percent": {
                    field: 100.0 * value
                    for field, value in selected_metrics(
                        previous["final"], current["final"]
                    ).items()
                },
            }
        )

    dt_to_reference = {
        str(record["steps"]): selected_metrics(record["final"], dt_reference["final"])
        for record in dt_records
    }
    grid_fraction = selected_metrics(dt_grid_coarse["final"], dt_grid_fine["final"])

    dt_projection_count = sum(record["total_projections_csv"] for record in dt_records)
    dt_max_residual = max(
        record["max_abs_energy_relative_residual_csv"] for record in dt_records
    )
    dt_logs_complete = all(record["log"]["status"] == "completed" for record in dt_records)
    dt_expected_rows = all(record["rows"] == record["steps"] for record in dt_records)
    dt_residual_gate = dt_max_residual < 1e-9

    # Keep copied artifact bookkeeping in the JSON, but omit rows_data from the
    # public summary so the file stays compact and remains easy to diff.
    public_records: dict[str, Any] = {}
    for stem, record in records.items():
        public_records[stem] = {key: value for key, value in record.items() if key != "rows_data"}

    summary: dict[str, Any] = {
        "schema": "coupled-thermal-validation-v1",
        "generated_by": {
            "script": Path(__file__).name,
            "source_sha256": sha256_file(Path(__file__).resolve()),
            "study_driver": {
                "file": "tools/study_coupled_thermal.cpp",
                "sha256": sha256_file(REPO_ROOT / "tools" / "study_coupled_thermal.cpp"),
                "role": "analysis-time source snapshot",
                "not_binary_build_provenance": True,
                "purpose": "numeric-setting provenance for the archived CSV/log studies",
            },
        },
        "scope": {
            "purpose": "reproducible local coupled thermal/fast-particle laboratory study",
            "device_prediction": False,
            "fair_four_fuel_ranking": False,
            "host_transport_or_BALDUR_acceptance": False,
            "excluded_pending_input_pattern": "coupled-pb-*-q*.csv",
            "interpretation": (
                "These artifacts measure the declared first-order local operator and its "
                "numerical sensitivity; they do not establish EXL performance or a fuel ranking."
            ),
            "stable_dt_gates_measured": True,
            "handoff": {
                "enabled": True,
                "actual_projection_count": dt_projection_count,
                "all_dt_projection_counts_zero": dt_projection_count == 0,
                "interpretation": (
                    "The continuous hot-source runs recorded no handoff projections. "
                    "A thermal-scale kinetic component is not automatically fluid ash."
                ),
            },
        },
        "numeric_settings": {
            "density_m3": 1.0e22,
            "electron_density_m3": 1.0e22,
            "carbon_density_m3": 1.0e19,
            "carbon_density_fraction_of_electron_density": 1.0e-3,
            "carbon_mass_u": 12.0,
            "carbon_mass_kg": 12.0 * 1.66053906892e-27,
            "carbon_mean_charge_squared": 36.0,
            "coulomb_log": 15.0,
            "thermal_charge_squared_network": [1.0, 1.0, 1.0, 4.0, 4.0, 25.0],
            "external_birth_enabled": False,
            "physical_escape_enabled": False,
            "handoff_enabled": True,
            "handoff_max_L1": 1.0e-3,
            "handoff_max_mean_error": 1.0e-3,
            "source_rate_error_limit": 1.0e-5,
            "source_debit_error_limit": 1.0e-5,
            "birth_options": {
                "relative_max_J": 2.5 * MEV_J,
                "cm_max_kT": 40.0,
                "ground_state_q_J": 0.09184 * MEV_J,
                "cutoff_J": 0.001 * MEV_J,
                "l1_fraction": 0.76,
                "relative_phase": 0.0,
                "narrow_peak_fraction": 0.051,
                "continuum_peak_scale": 1.0,
                "continuation": 1,
                "pb_low": 0,
                "remainder_policy": 0,
                "broad_mode": 13,
                "fsci_policy": 0,
                "relative_order": 16,
                "cm_order": 12,
                "nq": 8,
                "ncos": 8,
            },
            "DT": {
                "duration_s": 0.01,
                "initial_Ti_keV": 20.0,
                "initial_Te_keV": 5.0,
                "cells": 800,
                "steps": [16, 32, 64, 128],
            },
            "pB": {
                "duration_s": 0.001,
                "initial_Ti_keV": 100.0,
                "initial_Te_keV": 5.0,
                "cells": 800,
                "successful_nq_ncos": 8,
                "successful_first_edge_keV": 1.0e-18,
                "first_failure_first_edge_keV": 1.0e-10,
            },
        },
        "units": {
            "temperature": "keV (kT energy units)",
            "energy_density": "J m^-3",
            "number_density": "m^-3",
            "time": "s",
            "relative_residual": "dimensionless",
        },
        "runs": public_records,
        "guard_fix_parity": {
            "reference_run": records[PARITY_REFERENCE_STEM]["csv"]["file"],
            "guard_fix_run": {
                key: value for key, value in parity_record.items() if key != "rows_data"
            },
            "comparison": parity_comparison,
        },
        "dt_timestep_convergence": {
            "reference": {"steps": 128, "cells": 800, "file": dt_reference["csv"]["file"]},
            "adjacent": adjacent,
            "relative_to_128_step_fraction": dt_to_reference,
            "relative_to_128_step_percent": {
                steps: {field: 100.0 * value for field, value in metrics.items()}
                for steps, metrics in dt_to_reference.items()
            },
        },
        "dt_grid_convergence": {
            "coarse": {"steps": 128, "cells": 800, "file": dt_grid_coarse["csv"]["file"]},
            "fine": {"steps": 128, "cells": 1600, "file": dt_grid_fine["csv"]["file"]},
            "coarse_to_fine_relative_change_fraction": grid_fraction,
            "coarse_to_fine_relative_change_percent": {
                field: 100.0 * value for field, value in grid_fraction.items()
            },
        },
        "dt_validation_gates": {
            "completion_logs_present": dt_logs_complete,
            "expected_accepted_rows": dt_expected_rows,
            "all_rows_finite": True,
            "max_abs_energy_relative_residual": dt_max_residual,
            "energy_relative_residual_limit": 1.0e-9,
            "energy_relative_residual_pass": dt_residual_gate,
            "actual_handoff_projections": dt_projection_count,
            "handoff_projection_pass": dt_projection_count == 0,
            "all_stable_dt_gates_pass": dt_logs_complete and dt_expected_rows and dt_residual_gate,
        },
        "companion_runs": {
            "interpretation": (
                "DD, D-He3, and pB artifacts are retained as completed companion studies. "
                "They are not used to claim an EXL comparison or rank the fuels."
            ),
            "files": [
                records["coupled-dd-s16-n800"]["csv"]["file"],
                records["coupled-dhe3-s16-n800"]["csv"]["file"],
                records["coupled-pb-s4-n800"]["csv"]["file"],
            ],
        },
        "pB_grid_scope": {
            "preserved_first_failure": True,
            "failure_csv": records["coupled-pb-s4-n800-first-grid-failure"]["csv"]["file"],
            "failure_log": records["coupled-pb-s4-n800-first-grid-failure"]["log_file"]["file"],
            "requested_steps": records["coupled-pb-s4-n800-first-grid-failure"]["steps"],
            "completed_steps": records["coupled-pb-s4-n800-first-grid-failure"]["completed_steps"],
            "dt_s": records["coupled-pb-s4-n800-first-grid-failure"]["dt_s"],
            "failed_step": records["coupled-pb-s4-n800-first-grid-failure"]["log"]["failed_step"],
            "return_status": records["coupled-pb-s4-n800-first-grid-failure"]["log"]["return_status"],
            "failure_reason": "charged lab source extended below the first cell center",
            "grid_lower_edge_change": {
                "first_failure_first_edge_keV": 1.0e-10,
                "successful_first_edge_keV": 1.0e-18,
                "description": (
                    "The successful pB study lowers the first edge from 1e-10 to 1e-18 keV "
                    "to contain the declared source support. This is an explicit grid change, "
                    "not clipping or silently thermalizing a spill, and is not a universal guarantee."
                ),
            },
        },
        "figures": plot_files,
        "constants": {"keV_J": KEV_J, "MeV_J": MEV_J},
    }
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="directory containing the completed coupled-thermal CSV/log files (default: /tmp)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "artifact directory for copied studies, summary, and figures "
            f"(default: {DEFAULT_OUTPUT_DIR})"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output.expanduser()
    if not output_dir.is_absolute():
        output_dir = (Path.cwd() / output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    records: dict[str, dict[str, Any]] = {}
    for spec in RUN_SPECS:
        records[spec.stem] = run_record(spec, input_dir, output_dir)

    parity_record = run_record(PARITY_SPEC, input_dir, output_dir)
    parity_comparison = compare_common_csv_values(
        records[PARITY_REFERENCE_STEM], parity_record
    )
    plot_files = make_plot(records, output_dir)
    summary = build_summary(
        records, parity_record, parity_comparison, output_dir, plot_files
    )
    summary_path = output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )

    print(f"wrote {summary_path}")
    for name in plot_files:
        print(f"wrote {output_dir / name}")
    print(
        "DT gates: "
        f"residual<{summary['dt_validation_gates']['energy_relative_residual_limit']:.0e}="
        f"{summary['dt_validation_gates']['energy_relative_residual_pass']}, "
        f"handoff projections={summary['dt_validation_gates']['actual_handoff_projections']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
