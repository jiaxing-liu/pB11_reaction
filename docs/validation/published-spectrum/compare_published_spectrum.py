import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=list(csv.DictReader(open(sys.argv[3])))
s=list(csv.DictReader(open(sys.argv[1])))
fig,axes=plt.subplots(1,3,figsize=(12,3.9));metrics=[]
for panel,mode in enumerate((1,3,13)):
 a=[x for x in s if int(x['mode'])==mode];E=np.array([float(x['energy_keV']) for x in a]);v=np.array([float(x['number_per_event']) for x in a]);v/=np.trapezoid(v,E)
 ax=axes[panel];ax.plot(E/1000,v*1000,label='Library',c='#245f9b')
 for c,label,color in [('data','Published data','#202020'),('model','Published model','#d17b24')]:
  q=[x for x in r if int(x['panel'])==panel and x['curve']==c];x=np.array([float(x['energy_keV']) for x in q]);y=np.array([float(x['relative_height']) for x in q]);y/=np.trapezoid(y,x)
  interp=np.interp(E,x,y);interp/=np.trapezoid(interp,E)
  metrics.append({'mode':mode,'reference':c,'half_L1':float(.5*np.trapezoid(abs(v-interp),E)),'mean_keV_library':float(np.trapezoid(E*v,E)),'mean_keV_reference':float(np.trapezoid(E*interp,E))})
  ax.plot(x/1000,y*1000,label=label,c=color,ls='--' if c=='model' else '-',lw=1)
 if mode==13:
  raw=list(csv.DictReader(open(sys.argv[4])))
  raw_y=np.array([float(x['number_per_event']) for x in raw]);raw_y/=np.trapezoid(raw_y,E)
  ax.plot(E/1000,raw_y*1000,label='Raw-coefficient hypothesis',c='#43804e',ls=':',lw=1.5)
 ax.set(xlabel='CM alpha energy (MeV)',ylabel='Area-normalized spectrum (1/MeV)',title=f'l = {mode}' if mode!=13 else 'Mixture: k=0.76, phase=0.67 x 2pi');ax.legend(fontsize=8);ax.spines[['top','right']].set_visible(False)
fig.suptitle('Kuhlwein 2021 Fig. 3: published spectra and declared library model',fontsize=12)
fig.text(.02,.02,'Curves extracted from published vector paths. Data include event selection; no detector correction or fit-parameter inference is claimed.',fontsize=8)
fig.tight_layout(rect=(0,.065,1,.94));fig.savefig(sys.argv[2]+'.png',dpi=170);fig.savefig(sys.argv[2]+'.pdf')
print(json.dumps(metrics,indent=2))
