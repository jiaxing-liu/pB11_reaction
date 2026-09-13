from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).parent
a=np.genfromtxt(r/'pb-birth-study.csv',delimiter=',',names=True)
fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
for f,label,color in [(0,'Pure alpha1 (l=2, FSCI)','#186794'),(.05,'5% alpha0 + 95% alpha1','#d46a31'),(1,'Pure alpha0 (narrow-width)','#328463')]:
 v=a[a['alpha0_fraction']==f];axs[0].plot(v['E_MeV'],v['N_per_event']/3/.06,label=label,color=color)
 if f<1:axs[1].plot(v['E_MeV'],v['N_per_event']/3/.06,label=label,color=color)
axs[0].set_yscale('log');axs[0].set_ylim(1e-4,10);axs[0].set_title('Full branches; log scale')
axs[1].set_title('Illustrative 5% branch contribution')
for ax in axs:
 ax.set(xlabel='Alpha CM kinetic energy (MeV)',ylabel='Probability density (1/MeV)',xlim=(0,6));ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.suptitle('Conditional pB CM births: A=8.84 MeV, q(gs)=91.84 keV\n60 keV grid; 5% is illustrative, not a fitted branching fraction',fontsize=12)
for ext in ['png','pdf']:fig.savefig(r/f'pb-birth-comparison.{ext}',dpi=180)
