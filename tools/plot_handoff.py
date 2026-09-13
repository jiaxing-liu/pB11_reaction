#!/usr/bin/env python3
"""Plot a trace-alpha handoff reference; not an integrated device simulation."""
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
with a.input.open() as h:r=list(csv.DictReader(h))
d={k:np.array([float(x[k]) for x in r]) for k in r[0]};assert d['time_s'][-1]==1
t=d['time_s'];plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18})
fig,ax=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
for key,label,color in [('ST_N_fraction','Suprathermal kinetic','#2166ac'),('TH_N_fraction','Thermal-scale kinetic','#d6604d'),('fluid_N_fraction','Projected thermal fluid','#1b7837')]:ax[0].step(t,d[key],where='post',label=label,color=color)
ax[0].plot(t,d['ST_N_fraction']+d['TH_N_fraction']+d['fluid_N_fraction'],'k--',lw=1,label='Total alpha number');ax[0].set(ylabel='Number / initial alpha number',xlabel='Time (s)',ylim=(-.02,1.05),title='Helium classification; no new nuclear births');ax[0].legend(fontsize=8,loc='center left')
for key,label,color in [('electron_heat_fraction','Electron deposited energy','#b2182b'),('ion_heat_fraction','Ion deposited energy (incl. correction)','#2166ac'),('kinetic_U_fraction','Remaining kinetic energy','#d6604d'),('fluid_U_fraction','Thermal-fluid alpha energy','#1b7837')]:ax[1].plot(t,d[key],label=label,color=color)
ax[1].set(ylabel='Energy / initial alpha energy',xlabel='Time (s)',ylim=(-.02,1.05),title='Energy moves with the particles');ax[1].legend(fontsize=8,loc='center right')
fig.suptitle('Measured Maxwellian handoff: public API reference test\n3 MeV alpha pulse; fixed 10 keV e/D/T baths; 4000 cells, 2000 steps',fontsize=12)
fig.supxlabel('Fluid jump is a conservative representation change, not a reaction burst. Both projection tolerances = 0.001.',fontsize=9)
a.output.parent.mkdir(parents=True,exist_ok=True)
for ext in ['png','pdf']:fig.savefig(a.output.with_suffix('.'+ext),dpi=180)
