import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).parent

def read(stem):
 data=np.genfromtxt(root/(stem+'.csv'),delimiter=',',names=True)
 ledger=np.loadtxt(root/(stem+'.log'),ndmin=2)
 cases={}
 for row in ledger:
  mode,policy=map(int,row[:2]);grid=data[(data['mode']==mode)&(data['policy']==policy)]['number_per_event']
  cases[(mode,policy)]={'grid':grid,'all':np.r_[grid,row[4:6]],'ledger':row}
 return cases

def distance(a,b):return float(np.sum(np.abs(a['all']-b['all']))/6)
levels={n:read(f'fsci-spectrum-n{n}') for n in [128,256,512,1024]}
cut=read('fsci-spectrum-cut10-n1024');mix=read('fsci-spectrum-mix-n1024');fixed=read('fsci-spectrum-fixed-raw-n1024')
results={'metric':'half-L1 on 100 energy bins plus below/above categories, divided by 3 alphas/event','convergence':{},'cutoff_1_to_10_keV':{},'policy_difference':{},'number_energy_residuals':{}}
for case in levels[512]:
 label=f'l{case[0]}_policy{case[1]}'
 results['convergence'][label]={'128_to_256':distance(levels[128][case],levels[256][case]),'256_to_512':distance(levels[256][case],levels[512][case]),'512_to_1024':distance(levels[512][case],levels[1024][case])}
 results['cutoff_1_to_10_keV'][label]=distance(levels[1024][case],cut[case])
 row=levels[1024][case]['ledger']
 results['number_energy_residuals'][label]={'number_per_event':float(sum(row[3:6])-3),'energy_relative':float(sum(row[6:9])/row[2]-1),'below_N':float(row[4]),'above_N':float(row[5])}
results['policy_difference']['l2']=distance(levels[1024][(2,0)],levels[1024][(2,1)])
results['policy_difference']['l13_fixed_unit_k']=distance(mix[(13,0)],mix[(13,1)])
results['policy_difference']['l13_fixed_raw_coefficients']=distance(mix[(13,0)],fixed[(13,1)])
low_delta_512=levels[512][(2,1)]['all']-levels[512][(2,0)]['all']
low_delta_1024=levels[1024][(2,1)]['all']-levels[1024][(2,0)]['all']
mix512=read('fsci-spectrum-mix-full')
broad_delta_512=mix512[(13,1)]['all']-mix512[(13,0)]['all']
broad_delta_1024=mix[(13,1)]['all']-mix[(13,0)]['all']
results['paired_policy_difference_refinement_512_to_1024']={
 'l2':float(np.sum(np.abs(low_delta_512-low_delta_1024))/6),
 'l13_fixed_unit_k':float(np.sum(np.abs(broad_delta_512-broad_delta_1024))/6)}
(root/'fsci-spectrum-comparison.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained');E=(np.arange(100)+.5)*.06
axs[0].plot(E,levels[1024][(2,0)]['grid']/.18,label='No added FSCI',color='#16689f')
axs[0].plot(E,levels[1024][(2,1)]['grid']/.18,label='Refsgaard Model II, 16 fm',color='#c45c25')
axs[0].set_title(r'Primary $l=2$, $A=8.84$ MeV')
axs[1].plot(E,mix[(13,0)]['grid']/.18,label='No added FSCI, k=0.76',color='#16689f')
axs[1].plot(E,mix[(13,1)]['grid']/.18,label='FSCI, same unit-basis k',color='#c45c25')
axs[1].plot(E,fixed[(13,1)]['grid']/.18,label='FSCI, same raw coefficients',color='#438448',ls='--')
axs[1].set_title(r'Coherent $l=1,3$, $A=9.30$ MeV')
for ax in axs:
 ax.set_xlabel('Alpha kinetic energy in parent CM (MeV)');ax.set_ylabel(r'Probability density (MeV$^{-1}$)')
 ax.grid(alpha=.2);ax.legend(fontsize=8);ax.set_xlim(0,6)
fig.suptitle('Conditional CM source models: final-state Coulomb sensitivity\n1024-point quadrature in each variable; no laboratory or host feedback',fontsize=12)
for ext in ['png','pdf']:fig.savefig(root/f'fsci-source-comparison.{ext}',dpi=180)
