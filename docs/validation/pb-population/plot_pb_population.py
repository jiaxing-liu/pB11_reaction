from pathlib import Path
import csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).parent;rows=list(csv.DictReader((r/'pb-population-study.csv').open()))
T=np.array([2.,10.,50.,150.,300.,1000.]);F=np.zeros((5,6))
for j,t in enumerate(T):
 case=[v for v in rows if int(v['kind'])==0 and float(v['Ta_or_beam_keV'])==t];K=float(case[0]['K'])
 for v in case:F[int(v['component']),j]=float(v['K'])/K
fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
for i,label,color in [(1,'Ground-state peak','#c16029'),(2,'Narrow-fit remainder','#236c9a'),(3,'Other remainder','#408858')]:axs[0].semilogx(T,100*F[i],'-o',label=label,color=color)
axs[1].semilogx(T,100*F[4],'-o',color='#7b5796')
axs[0].set(ylabel='Fraction of reacting events (%)',title='Three source populations',ylim=(0,105));axs[0].legend(fontsize=8)
axs[1].set(ylabel='Fraction of reacting events (%)',title='Continuum table endpoint extension used',ylim=(0,105))
for ax in axs:ax.set_xlabel('Proton kT (keV), boron kT = 0.7 proton kT');ax.grid(alpha=.2)
fig.suptitle('Explicit effective pB population model: thermal reaction weighting\nPeak parameter 5.1%, continuum scale 1; no spectrum or device calibration',fontsize=12)
for ext in ['png','pdf']:fig.savefig(r/f'pb-population-weights.{ext}',dpi=180)
