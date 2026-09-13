#!/usr/bin/env python3
"""Fixed-bath alpha operator study; these are NOT BALDUR/device results."""
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--cold',type=Path,required=True);p.add_argument('--finite',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
def read(path):
 with path.open() as f:r=list(csv.DictReader(f))
 d={k:np.array([float(x[k]) for x in r]) for k in r[0]}
 assert d['time_s'][-1]==1.
 return d
c,f=read(a.cold),read(a.finite)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,'figure.dpi':130})
fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
t=f['time_s'];tc=c['time_s'];blue='#2166ac';red='#b2182b';green='#1b7837'
ax=axs[0,0];ax.plot(t,f['st_number_fraction'],label='Suprathermal component',color=blue);ax.plot(t,f['th_number_fraction'],label='Thermal-scale kinetic component',color=green);ax.plot(t,f['st_number_fraction']+f['th_number_fraction'],'k--',label='Total retained alpha');ax.set(ylabel='Number / initial alpha number',title='Finite-temperature split: particle inventory',ylim=(-.02,1.05));ax.legend(fontsize=8)
ax=axs[0,1];ax.semilogy(tc,c['normalized_half_L1_vs_FULL'],color=red,label='Cold-ion approximation');ax.semilogy(t,np.maximum(f['normalized_half_L1_vs_FULL'],1e-17),color=blue,label='Finite-temperature split');ax.set(ylabel='Distribution half-L1 / initial number',title='Difference from full collision equation',ylim=(1e-17,1));ax.legend(fontsize=8)
ax=axs[1,0];ax.plot(t,f['split_e_heat_fraction'],color=red,label='Electron deposited energy');ax.plot(t,f['split_D_heat_fraction']+f['split_T_heat_fraction'],color=blue,label='Ion deposited energy');ax.plot(t,f['st_energy_fraction']+f['th_energy_fraction'],color=green,label='Remaining kinetic energy');ax.set(ylabel='Energy / initial alpha energy',title='Finite-temperature split: energy budget',ylim=(-.02,1.05));ax.legend(fontsize=8)
ax=axs[1,1];ax.semilogy(t,np.maximum(f['th_maxwellian_half_L1'],1e-17),color=green,label='TH vs 10 keV sampled Maxwellian');ax.set(ylabel='Normalized TH distribution half-L1',title='A thermal-scale label does not imply equilibrium');ax.legend(fontsize=8)
for ax in axs.flat:ax.set_xlabel('Time (s)')
fig.suptitle('3 MeV trace alpha relaxation in fixed 10 keV e/D/T baths\n1000 energy cells, 1000 steps; internal transfer is not fluid helium ash',fontsize=13)
a.output.parent.mkdir(parents=True,exist_ok=True)
for ext in ('png','pdf'):fig.savefig(a.output.with_suffix('.'+ext),dpi=180)
