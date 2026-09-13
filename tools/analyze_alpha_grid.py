#!/usr/bin/env python3
"""Summarize numerical sensitivity of the declared CM amplitude model."""
import argparse,csv,json,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('csv',type=Path);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(args.csv.open())]
def select(mode,A,n,cut=.001):
 return next(r for r in rows if r['mode']==mode and abs(r['A_MeV']-A)<1e-8 and r['nq']==n and abs(r['cutoff_MeV']-cut)<1e-10)
def l1(a,b):
 return (sum(abs(a[f'N_bin{i}']-b[f'N_bin{i}']) for i in range(32))+abs(a['below_N']-b['below_N'])+abs(a['above_N']-b['above_N']))/3
cases=[]
for A in (8.829,9.3,12.):
 for mode in (1,2,3,13):
  fine=select(mode,A,256)
  refinement=[]
  for n in (32,64,128):
   old,new=select(mode,A,n),select(mode,A,2*n)
   refinement.append({'n_to_2n':[n,2*n],'normalized_L1_including_spill':l1(old,new),'raw_normalization_relative_change':abs(old['normalization_J2']/new['normalization_J2']-1)})
  cutoff=[{'cutoff_MeV':c,'normalized_L1_including_spill':l1(fine,select(mode,A,256,c)),'raw_normalization_relative_change':abs(select(mode,A,256,c)['normalization_J2']/fine['normalization_J2']-1)} for c in (.002,.004,.01)]
  cases.append({'A_MeV':A,'mode':mode,'refinement':refinement,'cutoff':cutoff})
result={'scope':'Numerical CM spectrum only; no experimental, branch, incident-energy or host validation. Histogram has 32 arithmetic-center cells on 0..8 MeV; explicit below/above spill retained. Mixture uses unit-normalized bases k=.76, delta=.67*2pi; this is a declared convention, not independently established fitted branching.', 'cases':cases,'rows':len(rows),'max_number_residual':max(abs(r['N_residual']) for r in rows),'max_energy_residual_MeV':max(abs(r['E_residual_MeV']) for r in rows),'max_128_to_256_L1':max(c['refinement'][-1]['normalized_L1_including_spill'] for c in cases),'max_cutoff_L1':max(cut['normalized_L1_including_spill'] for c in cases for cut in c['cutoff'])}
result['numerical_checks_pass']=result['max_128_to_256_L1']<.02 and result['max_cutoff_L1']<1e-9 and result['max_number_residual']<3e-12 and result['max_energy_residual_MeV']<12e-12
args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
raise SystemExit(0 if result['numerical_checks_pass'] else 1)
