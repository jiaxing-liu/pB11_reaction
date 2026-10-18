#!/usr/bin/env python3
"""Current-source pB time refinement; hot local benchmark, not an EXL case."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

base = Path(__file__).resolve().parent
binary = base / 'study_coupled_thermal'
manifest = {
    'id': 'D869', 'status': 'running', 'pid': os.getpid(),
    'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    'validated_kernel_source_commit': 'ff6c200d6b61bfd13e004ed143e7c614dad810fb',
    'external_evidence_commit': 'f8524cc32af8686a6c2749c9c5c9e60a10af7aae',
    'kernel_identity': '54b35a3cf27a6224908c5fcc9a0bc61ade2ba46a237074e68f087f659af928dd',
    'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
    'driver_sha256': hashlib.sha256((base/'study_coupled_thermal.cpp').read_bytes()).hexdigest(),
    'scope': 'pB local evolving-bath time refinement; not EXL-50U or a four-fuel comparison',
    'parameters': {'ne_m3': 1e22, 'nC_over_ne': 0.001, 'initial_Te_keV': 5,
                   'initial_Ti_keV': 100, 'cells': 800, 'duration_s': 0.001,
                   'handoff': 1, 'alpha_order': 8, 'steps': [64,128,256],
                   'wall_limit_s_per_resolution': 43200},
    'runs': [], 'physics_accepted': False, 'm3_accepted': False,
}

def save():
    tmp=base/'manifest.json.tmp'
    tmp.write_text(json.dumps(manifest,indent=2)+'\n')
    tmp.replace(base/'manifest.json')

save()
try:
    for steps in [64,128,256]:
        out=base/f'steps-{steps}.csv'; log=base/f'steps-{steps}.log'
        if out.exists() or log.exists():
            raise RuntimeError('refusing existing output; inspect terminal metadata before resuming')
        cmd=['/usr/bin/time','-v','-o',str(base/f'resources-{steps}.txt'),
             'timeout','--kill-after=60s','43200s','stdbuf','-oL',str(binary),
             str(steps),'800','pb','0.001','1','8']
        entry={'steps':steps,'command':cmd,'status':'running'}
        manifest['runs'].append(entry);save()
        print('START',steps,'steps; 12-hour wall guard',flush=True)
        start=time.monotonic()
        with out.open('x') as stdout, log.open('x') as stderr:
            p=subprocess.run(cmd,stdout=stdout,stderr=stderr,start_new_session=True)
        entry.update(returncode=p.returncode,elapsed_s=time.monotonic()-start)
        if p.returncode:
            entry['status']='failed_or_interrupted';save()
            raise RuntimeError(f'{steps} steps exited {p.returncode}; preserve partial output')
        with out.open() as stream: rows=list(csv.DictReader(stream))
        if len(rows)!=steps or not math.isclose(float(rows[-1]['time_s']),0.001,rel_tol=1e-12):
            raise RuntimeError('incomplete benchmark endpoint')
        endpoint={k:float(v) for k,v in rows[-1].items()}
        if not all(math.isfinite(v) for v in endpoint.values()):
            raise RuntimeError('nonfinite benchmark endpoint')
        residual=max(abs(float(row['energy_relative_residual'])) for row in rows)
        if residual>=1e-9:raise RuntimeError('independent accumulated energy closure fails')
        entry.update(status='complete',endpoint=endpoint,max_energy_residual=residual,
                     csv_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
        save();print('COMPLETE',steps,'elapsed_s',entry['elapsed_s'],flush=True)
    quantities=['Ue_J_m3','Ui_J_m3','fast_energy_J_m3','Q_J_m3',
                'electron_heat_J_m3','network_ion_heat_J_m3','inert_heat_J_m3']
    comparisons={}
    for a,b in zip(manifest['runs'],manifest['runs'][1:]):
        data={}
        for key in quantities:
            x,y=a['endpoint'][key],b['endpoint'][key]
            data[key]={'absolute_difference':abs(y-x),
                       'relative_to_finer':abs(y-x)/abs(y) if y else None}
        comparisons[f'{a["steps"]}_to_{b["steps"]}']=data
    manifest.update(status='complete',comparisons=comparisons,
                    limitations=['fixed energy/angular grids; temporal differences are not rigorous error bounds',
                                 'finite-source spectrum assumptions and EXL host coupling remain separate'])
    save();print('TIME REFINEMENT COMPLETE',flush=True)
except Exception as exc:
    manifest.update(status='stopped_unaccepted',error=str(exc));save();raise
