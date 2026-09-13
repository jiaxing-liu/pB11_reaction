from pathlib import Path
import csv,json,hashlib,shutil
base=Path('/home/cloud/research/pB11_reaction/docs/validation/coupled-thermal/pb-refinement')
base.mkdir(parents=True,exist_ok=True)
runs={}
for steps,cells,q in [(8,800,8),(16,800,8),(32,800,8),(64,800,8),(8,1600,8),(8,800,16)]:
 name=f'coupled-pb-s{steps}-n{cells}-q{q}'
 src=Path('/tmp')/name
 log=src.with_suffix('.log').read_text()
 if f'fuel=pb steps={steps} cells={cells} duration=0.001 ' not in log:raise RuntimeError(f'Incomplete: {name}')
 rows=list(csv.DictReader(src.with_suffix('.csv').open()))
 if len(rows)!=steps or float(rows[-1]['time_s'])!=.001:raise RuntimeError(name)
 for suffix in ('.csv','.log'):shutil.copyfile(src.with_suffix(suffix),base/(name+suffix))
 last={k:float(v) for k,v in rows[-1].items()}
 runs[name]={'steps':steps,'cells':cells,'nq_ncos':q,'rows':len(rows),'final':last,'sha256':{s:hashlib.sha256((base/(name+s)).read_bytes()).hexdigest() for s in ('.csv','.log')},'max_abs_energy_relative_residual':max(abs(float(r['energy_relative_residual'])) for r in rows)}
def pair(a,b,fields):
 A,B=runs[a]['final'],runs[b]['final'];return {'coarse':a,'fine':b,'fine_minus_coarse':{k:{'absolute':B[k]-A[k],'relative_to_fine':(B[k]-A[k])/B[k] if B[k] else None} for k in fields}}
fields=['Ue_J_m3','Ui_J_m3','fast_energy_J_m3','inert_heat_J_m3','Q_J_m3']
data={'scope':'Local classical high-density thermal/fast feedback stress study; not EXL or fair fuel ranking. Angular proxy accuracy not established by these observable changes.', 'units':'energies and heat J/m^3; relative quantities dimensionless','runs':runs,'time_pairs':[pair(f'coupled-pb-s{a}-n800-q8',f'coupled-pb-s{b}-n800-q8',fields) for a,b in [(8,16),(16,32),(32,64)]], 'source_heat_time_pair':pair('coupled-pb-s32-n800-q8','coupled-pb-s64-n800-q8',['electron_heat_J_m3','network_ion_heat_J_m3','removed_thermal_energy_J_m3','inert_heat_J_m3']), 'grid_pair':pair('coupled-pb-s8-n800-q8','coupled-pb-s8-n1600-q8',fields),'angular_pair':pair('coupled-pb-s8-n800-q8','coupled-pb-s8-n800-q16',fields), 'provenance':{'compiled_library':'pre inventory-roundoff-guard correction; successful numerical path unchanged; see parent parity record','study_binary_sha256_at_launch':'42b81eb82a1c0b46f71047873cb619aa6043c5303e10cb5b1d3a87ccc787a529','pre_guard_static_library_sha256':'3611811463b5e22a622dd0933ca5aace8500f28d230cf374babc851da590ec20','extra_heat_columns':'32 and64 steps only; earlier CSV layouts retained verbatim'}}
(base/'summary.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'time32to64':data['time_pairs'][-1], 'heat32to64':data['source_heat_time_pair']},indent=2))
