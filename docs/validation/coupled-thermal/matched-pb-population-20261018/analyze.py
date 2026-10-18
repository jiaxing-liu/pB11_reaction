"""Summarize actual selected accepted profiles; never integrate sparse snapshots in time."""
import csv,json,math,hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
profiles=list(csv.DictReader((root/'profiles.csv').open()))
populations=list(csv.DictReader((root/'populations.csv').open()))
assert len(profiles)==len(populations)==200
by_step={}
for p,r in zip(profiles,populations):
 assert r['label']==f"s{p['step']}_z{p['zone']}"
 assert math.isclose(float(r['kTa_J']),float(p['Ti_keV'])*1.6021e-16,rel_tol=1e-14)
 for k,v in r.items():
  if k!='label' and v:assert math.isfinite(float(v))
 rate=float(r['total_reactivity_m3_s'])
 components=sum(float(r[k]) for k in ['alpha0_peak_reactivity_m3_s','narrow_remainder_reactivity_m3_s','other_remainder_reactivity_m3_s'])
 assert abs(components-rate)<=1e-12*rate
 extra=float(r['continuum_extrapolated_reactivity_m3_s']);assert 0<=extra<=rate*(1+1e-12)
 by_step.setdefault(int(p['step']),[]).append((p,r))
summary=[]
fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for step,rows in sorted(by_step.items()):
 assert len(rows)==50 and [int(p['zone']) for p,r in rows]==list(range(2,52))
 weights=[float(p['nH_m3'])*float(p['nB_m3'])*float(p['volume_m3']) for p,r in rows]
 total=sum(w*float(r['total_reactivity_m3_s']) for w,(p,r) in zip(weights,rows))
 extra=sum(w*float(r['continuum_extrapolated_reactivity_m3_s']) for w,(p,r) in zip(weights,rows))
 alpha=sum(w*float(r['alpha0_peak_reactivity_m3_s']) for w,(p,r) in zip(weights,rows))
 t=float(rows[0][0]['time_s']);x=[float(p['minor_outer_m']) for p,r in rows]
 summary.append(dict(step=step,time_s=t,instantaneous_model_event_rate_s_inverse=total,
   source_weighted_continuum_extrapolated_fraction=extra/total if total>0 else None,
   source_weighted_alpha0_peak_fraction=alpha/total if total>0 else None,
   min_native_Ti_keV=min(float(p['Ti_keV']) for p,r in rows),max_native_Ti_keV=max(float(p['Ti_keV']) for p,r in rows)))
 axes[0].plot(x,[float(p['Ti_keV']) for p,r in rows],label=f'{t*1000:.3g} ms')
 axes[1].plot(x,[100*float(r['continuum_extrapolated_event_fraction']) if r['continuum_extrapolated_event_fraction'] else float('nan') for p,r in rows])
 axes[2].plot(x,[100*float(r['alpha0_peak_event_fraction']) if r['alpha0_peak_event_fraction'] else float('nan') for p,r in rows])
for ax,title,ylabel in zip(axes,['Accepted thermal ion temperature','Continuum table extrapolation','Effective alpha0 peak fraction'],['Native diagnostic / keV','Conditional event fraction / %','Conditional event fraction / %']):
 ax.set(title=title,xlabel='Zone outer minor radius / m',ylabel=ylabel);ax.grid(alpha=.25)
axes[0].legend()
fig.suptitle('Matched p–B11: accepted-state model-domain diagnostic\nEndpoint-held continuum proxy; not an uncertainty bound or exact source-entry replay')
for ext in ['png','pdf']:fig.savefig(root/f'population-domain.{ext}',dpi=150)
(root/'summary.json').write_text(json.dumps({'scope':'four accepted snapshots only; not time-integrated yield or source-entry equivalence',
 'snapshots':summary,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'profiles.csv',root/'populations.csv']},
 'limitations':['Continuum table endpoint hold is explicit model extrapolation, not data-supported deconvolution.',
 'Remainders are not declared pure alpha1. Extrapolated is an overlapping subset, not an extra population.',
 'Accepted fluid profiles are not the pre-trial library source-entry state.',
 'Sparse snapshots are not interpolated or integrated into a full-window yield.'],
 'm3_accepted':False},indent=2,allow_nan=False)+'\n')
