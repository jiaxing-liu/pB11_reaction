"""Plot only the completed, hash-pinned D869 temporal refinement study.

Usage: python3 plot.py PATH/TO/manifest.json NEW_OUTPUT_DIRECTORY
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    manifest_path, output = map(Path, sys.argv[1:])
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest['status'] != 'complete':
        raise ValueError('wait for all refinement resolutions to finish')
    runs = manifest['runs']
    steps = [run['steps'] for run in runs]
    if steps != manifest['parameters']['steps'] or len(runs) < 2 or any(
            b <= a for a, b in zip(steps, steps[1:])):
        raise ValueError('incomplete refinement sequence')
    data = []
    for run in runs:
        if run['status'] != 'complete' or run['returncode'] != 0:
            raise ValueError('nonterminal or failed resolution')
        path = manifest_path.parent / f"steps-{run['steps']}.csv"
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != run['csv_sha256']:
            raise ValueError('CSV hash mismatch')
        rows = [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(raw.decode().splitlines())]
        required = {'time_s', 'Te_keV', 'Ti_keV', 'Ue_J_m3', 'Ui_J_m3',
                    'fast_energy_J_m3', 'Q_J_m3', 'electron_heat_J_m3',
                    'network_ion_heat_J_m3', 'inert_heat_J_m3', 'energy_relative_residual'}
        if not rows or any(not required <= row.keys() for row in rows):
            raise ValueError('required cumulative heat or state field missing')
        if len(rows) != run['steps'] or not all(
                math.isfinite(value) for row in rows for value in row.values()):
            raise ValueError('missing or nonfinite history')
        times = [row['time_s'] for row in rows]
        if times[0] <= 0 or any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError('nonmonotone physical time')
        if not math.isclose(times[-1], manifest['parameters']['duration_s'], rel_tol=1e-12):
            raise ValueError('incomplete physical window')
        residual = max(abs(row['energy_relative_residual']) for row in rows)
        if residual >= 1e-9:
            raise ValueError('energy closure fails')
        data.append(rows)
    if manifest_path.read_bytes() != manifest_bytes:
        raise ValueError('manifest changed while reading')
    output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'figure.dpi': 150})
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    panels = [('fast_energy_J_m3', 'Fast-particle stored energy', 1e-3, 'kJ / m³'),
              ('Q_J_m3', 'Cumulative nuclear Q', 1e-3, 'kJ / m³'),
              ('electron_heat_J_m3', 'Cumulative electron heat', 1e-3, 'kJ / m³'),
              ('network_ion_heat_J_m3', 'Cumulative network-ion heat', 1, 'J / m³'),
              ('inert_heat_J_m3', 'Cumulative carbon-ion heat', 1, 'J / m³')]
    colors = ['#0072B2', '#D55E00', '#009E73']
    for index, (run, rows) in enumerate(zip(runs, data)):
        color = colors[index % len(colors)]
        time_ms = [1000 * row['time_s'] for row in rows]
        label = f"{run['steps']} steps"
        axes.flat[0].plot(time_ms, [row['Te_keV'] for row in rows], color=color, label=label)
        axes.flat[0].plot(time_ms, [row['Ti_keV'] for row in rows], color=color, linestyle='--')
        for ax, (key, title, scale, unit) in zip(list(axes.flat)[1:], panels):
            ax.plot(time_ms, [scale * row[key] for row in rows], color=color, label=label)
            ax.set(title=title, ylabel=unit)
    axes.flat[0].set(title='Bath temperatures (Te solid; Ti dashed)', ylabel='keV')
    axes.flat[0].legend(frameon=False)
    for ax in axes.flat:
        ax.set_xlabel('Time / ms')
        ax.grid(alpha=.2)
    fig.suptitle('p–B11 local time refinement — fixed energy/angular grids\n'
                 'ne=10²² m⁻³, initial Te=5 keV / Ti=100 keV; not an EXL-50U case', fontsize=12)
    fig.tight_layout(rect=(0, .02, 1, .91))
    fig.savefig(output / 'histories.png')
    fig.savefig(output / 'histories.pdf')
    plt.close(fig)

    # Compute directly from the pinned CSV endpoints, without a denominator floor.
    keys = [panel[0] for panel in panels] + ['Ue_J_m3', 'Ui_J_m3']
    comparisons = []
    for coarse, fine, x, y in zip(runs, runs[1:], data, data[1:]):
        values = {}
        for key in keys:
            delta = abs(y[-1][key] - x[-1][key])
            values[key] = {'absolute_difference_J_m3': delta,
                           'relative_to_finer': delta / abs(y[-1][key]) if y[-1][key] else None}
        comparisons.append({'coarse_steps': coarse['steps'], 'fine_steps': fine['steps'],
                            'endpoint_differences': values})
    fig, ax = plt.subplots(figsize=(9, 4.5))
    labels = ['Fast energy', 'Nuclear Q', 'Electron heat', 'Network-ion heat', 'Carbon heat']
    for index, comparison in enumerate(comparisons):
        values = [comparison['endpoint_differences'][key]['relative_to_finer'] for key in keys[:5]]
        if any(value is None for value in values):
            raise ValueError('undefined relative difference; inspect the numerical summary')
        width = .8 / len(comparisons)
        positions = [x - .4 + width * (index + .5) for x in range(5)]
        ax.bar(positions, [100 * value for value in values], width=width, color=colors[index],
               label=f"{comparison['coarse_steps']} → {comparison['fine_steps']} steps")
    ax.set_xticks(range(5), labels)
    ax.set(ylabel='Endpoint difference / finer endpoint (%)',
           title='Temporal refinement at 1 ms — numerical differences, not error bounds')
    ax.grid(alpha=.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / 'endpoint-differences.png')
    fig.savefig(output / 'endpoint-differences.pdf')
    plt.close(fig)
    summary = {'scope': manifest['scope'], 'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
               'comparisons': comparisons,
               'max_energy_residual': max(abs(row['energy_relative_residual']) for rows in data for row in rows),
               'm3_accepted': False,
               'limitations': ['No time-zero rows are invented.',
                               'Fixed energy/angular grids and source model; not a full convergence proof.',
                               'No EXL-50U or four-fuel comparison is represented.']}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(output)


if __name__ == '__main__':
    main()
