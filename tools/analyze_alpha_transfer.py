#!/usr/bin/env python3
"""Analyze a numerical first-passage study; no physical ash certification."""
import argparse,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
rows=[]
for path in a.input:
 for r in csv.DictReader(path.open()):
  r={k:float(v) for k,v in r.items()};r.setdefault('grid',0.);rows.append(r)
groups={}
for r in rows:groups.setdefault(tuple(r[k] for k in ('grid','cells','steps','cut_keV','penalty')),[]).append(r)
def crossing(rr,field,q):
 prev={field:0.,'time_s':0.}
 for r in rr:
  if r[field]>=q:return prev['time_s']+(r['time_s']-prev['time_s'])*(q-prev[field])/(r[field]-prev[field])
  prev=r
 return None
cases=[]
for key,rr in sorted(groups.items()):
 rr.sort(key=lambda r:r['time_s']);last=rr[-1]
 if last['time_s']!=1.:raise ValueError('Incomplete run: '+str(key))
 field='truncated_Maxwellian_overlap' if key[3]==0 else 'absorbed_N_fraction'
 case=dict(zip(('grid','cells','steps','cut_keV','penalty'),key))
 case.update(diagnostic=field,t50_s=crossing(rr,field,.5),t90_s=crossing(rr,field,.9),final_electron_heat_fraction=last['electron_heat_fraction'],final_ion_heat_after_Ti_mixing_fraction=last['ion_heat_after_Ti_mixing_fraction'],max_number_residual=max(abs(r['number_residual']) for r in rr),max_energy_residual=max(abs(r['energy_residual']) for r in rr))
 cases.append(case)
def pick(n,s,c):return next(x for x in cases if (x['grid'],x['cells'],x['steps'],x['cut_keV'])==(1,n,s,c))
coarse,fine,timefine=pick(2000,2000,0),pick(4000,2000,0),pick(4000,4000,0)
summary={'scope':'Trace alpha, fixed equal-temperature 10 keV electron/D/T baths, 3 MeV pulse, n_e=1e20 m^-3, lnLambda=15. Absorbed first-passage count and full-distribution Maxwellian overlap are different diagnostics, not interchangeable thermal-ash definitions.', 'reference_t90_grid_difference_s':abs(coarse['t90_s']-fine['t90_s']),'reference_t90_timestep_difference_s':abs(fine['t90_s']-timefine['t90_s']),'cases':cases,'max_number_residual':max(c['max_number_residual'] for c in cases),'max_energy_residual':max(c['max_energy_residual'] for c in cases)}
summary['conservation_pass']=summary['max_number_residual']<1e-9 and summary['max_energy_residual']<1e-9
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='cases'},indent=2))
