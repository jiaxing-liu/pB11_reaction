#!/usr/bin/env python3
"""Plot the documented diagnostic assumptions; this is not a data-error band."""
import csv,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open(sys.argv[1])))
fig,axes=plt.subplots(1,3,figsize=(14.5,4.7))
names=['p-B11','DD to T+p','DD to He3+n','DT','D-He3']
colors=['#8053ac','#2878b5','#41a484','#db7b2b','#ba4d69']
for ch in range(5):
 r=[x for x in rows if int(x['channel'])==ch]
 T=np.array([float(x['T_keV']) for x in r]);K=np.array([float(x['K_fit_m3_s']) for x in r])
 low=np.array([float(x['K_low_Sconst_m3_s']) for x in r]);hi=np.array([float(x['K_high_Sconst_m3_s']) for x in r])
 flat=np.array([float(x['K_high_flat_m3_s']) for x in r]);total=K+low+hi
 if ch>0:axes[0].loglog(T,low/total,'o-',ms=3,color=colors[ch],label=names[ch])
 axes[1].loglog(T,hi/total,'o-',ms=3,color=colors[ch],label=names[ch])
 axes[1].loglog(T,flat/(K+low+flat),':',lw=1,color=colors[ch])
 if ch==0:
  below=np.array([float(x['pb_K_below_SW140']) for x in r]);gap=np.array([float(x['pb_K_gap3480_5700']) for x in r])
  b22=np.array([float(x['pb_K_below_Becker22']) for x in r])
  axes[2].loglog(T,b22/K,'^-',ms=4,color='#d66c26',label='E < 22 keV (Becker floor)')
  axes[2].loglog(T,below/K,'o-',ms=4,color=colors[0],label='E < 140 keV')
  axes[2].loglog(T,gap/K,'s--',ms=3,color='#727272',label='3.48 < E < 5.70 MeV')
axes[0].set(xlim=(.01,3),ylim=(1e-12,1.1),title='Low-energy continuation',ylabel='Fraction of model reactivity')
axes[1].set(xlim=(10,500),ylim=(1e-20,1e-2),title='High-energy continuation',ylabel='Fraction of model reactivity')
axes[2].set(xlim=(.01,500),ylim=(1e-6,1.1),title='p-B11: data-coverage diagnostics',ylabel='Fraction of fit-window reactivity')
for ax in axes:
 ax.set_xlabel('Common ion temperature (keV)');ax.grid(True,which='major',alpha=.22)
 ax.legend(fontsize=8,loc='best');ax.spines[['top','right']].set_visible(False)
fig.suptitle('Thermal reactivity: fit-window and continuation diagnostics',fontsize=14,y=.98)
fig.text(.02,.053,'Continuation is an explicit sensitivity assumption: solid = constant endpoint S; dotted = constant high-energy cross section.',fontsize=9)
fig.text(.02,.018,'Right: Becker data start at 22 keV; SW at 140 keV. Coverage fractions are not error bounds. This is not a BALDUR / EXL run.',fontsize=9)
fig.tight_layout(rect=(0,.125,1,.93))
out=Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(out.with_suffix('.png'),dpi=180);fig.savefig(out.with_suffix('.pdf'))
