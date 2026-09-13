import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
src=Path(__file__).parent
p=src/'m3-source-state-study.csv'
x=np.genfromtxt(p,delimiter=',',names=True)
fig,axs=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
for ax in axs:
 ax.axvline(.2,color='.45',ls=':',lw=1,label='Checkpoint/restart')
 ax.set_xlabel('Time (s)');ax.grid(alpha=.2)
axs[0].plot(x['time_s'],x['born_N_m3']/1e12,label='Cumulative injection',color='black',ls='--')
axs[0].plot(x['time_s'],x['fast_N_m3']/1e12,label='Fast alpha inventory',color='#126eab')
axs[0].plot(x['time_s'],x['escaped_N_m3']/1e12,label='Cumulative escape',color='#c36123')
axs[0].set_ylabel(r'Number ($10^{12}$ m$^{-3}$)')
axs[0].set_title('Particle accounting')
axs[1].plot(x['time_s'],x['born_U_J_m3'],label='Cumulative injected energy',color='black',ls='--')
for key,label,color in [('fast_U_J_m3','Fast alpha energy','#126eab'),('escaped_U_J_m3','Escape energy','#c36123'),('heat_e_J_m3','Heat to electrons','#b02859'),('heat_p_J_m3','Heat to protons','#35804c')]:
 axs[1].plot(x['time_s'],x[key],label=label,color=color)
axs[1].set_ylabel(r'Energy (J m$^{-3}$)');axs[1].set_title('Energy accounting')
for ax in axs:ax.legend(fontsize=8,loc='upper left')
fig.suptitle('Library validation: continuous external alpha injection\nPrescribed changing baths; restarted and uninterrupted states are identical',fontsize=12)
for extension in ['png','pdf']:
 fig.savefig(src/f'source-state-validation.{extension}',dpi=180)
