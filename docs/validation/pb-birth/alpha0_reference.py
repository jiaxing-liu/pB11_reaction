import mpmath as mp, random, subprocess, json
from pathlib import Path
mp.mp.dps=80
mev=1.602176634e-13
mass=mp.mpf(float(mp.mpf('4.001506179129')*mp.mpf('1.66053906892e-27')))
rest=mass*299792458**2
rng=random.Random(141926)
cases=[]
for i in range(120):
 A=(8.84 if i<30 else rng.uniform(.1,12))*mev
 q=([0,A,.09184*mev,A/4,A*(1-1e-12),A*1e-12][i%6] if i<60 else rng.uniform(0,1)*A)
 n=[1,7,31,100][i%4]
 lo=[0,.2,.8][i%3]*A;hi=[1.1,.65,1.2][i%3]*A
 edges=[float(lo+(hi-lo)*(j/n)**1.3) for j in range(n+1)]
 cases.append((A,q,edges))
text=''.join(f'{len(e)-1} {A:.17g} {q:.17g} '+' '.join(f'{v:.17g}' for v in e)+'\n' for A,q,e in cases)
p=subprocess.run([str(Path(__file__).with_name('alpha0_reference_driver'))],input=text,text=True,capture_output=True,check=True)
rows=p.stdout.splitlines();assert len(rows)==len(cases)
maxn=maxe=maxend=0.;results=[]
for idx,((A0,q0,edges),line) in enumerate(zip(cases,rows)):
 vals=list(map(float,line.split()));assert vals[0]==0,(idx,vals[0])
 A,q=mp.mpf(A0),mp.mpf(q0);M=3*rest+A;B=2*rest+q
 p1=mp.sqrt((M*M-(rest+B)**2)*(M*M-(rest-B)**2))/(2*M)
 K=mp.sqrt(rest*rest+p1*p1)-rest;gb=mp.sqrt(1+(p1/B)**2)
 h=p1/B*mp.sqrt((B/2)**2-rest*rest)
 L=gb*(B/2)-rest-h;H=gb*(B/2)-rest+h
 c=[(mp.mpf(x)+mp.mpf(y))/2 for x,y in zip(edges[:-1],edges[1:])]
 b=[mp.mpf(0)]*len(c);bn=be=an=ae=mp.mpf(0)
 def delta(E,w):
  global bn,be,an,ae
  if E<c[0]:bn+=w;be+=w*E
  elif E>c[-1]:an+=w;ae+=w*E
  elif E==c[-1]:b[-1]+=w
  else:
   j=next(j for j in range(len(c)-1) if c[j]<=E<c[j+1]);b[j]+=w*(c[j+1]-E)/(c[j+1]-c[j]);b[j+1]+=w*(E-c[j])/(c[j+1]-c[j])
 delta(K,mp.mpf(1))
 if H==L:delta(L,mp.mpf(2))
 else:
  density=2/(H-L)
  def integral(l,u):return (density*(u-l),density*(u*u-l*l)/2) if u>l else (mp.mpf(0),mp.mpf(0))
  w,en=integral(L,min(H,c[0]));bn+=w;be+=en
  w,en=integral(max(L,c[-1]),H);an+=w;ae+=en
  for j in range(len(c)-1):
   l=max(L,c[j]);u=min(H,c[j+1])
   if u>l:
    # Independent antiderivatives of the two piecewise linear basis functions.
    den=(c[j+1]-c[j]);b[j]+=density*(c[j+1]*(u-l)-(u*u-l*l)/2)/den
    b[j+1]+=density*((u*u-l*l)/2-c[j]*(u-l))/den
 expected=[K,L,H,bn,be,an,ae]+b
 errs=[abs(mp.mpf(v)-r) for v,r in zip(vals[1:],expected)]
 ne=float(max(errs[3],errs[5],*errs[7:]));ee=float(max(errs[4],errs[6])/A);endpoint=float(max(errs[:3])/A)
 assert ne<1e-11 and ee<1e-12 and endpoint<1e-14,(idx,ne,ee,endpoint)
 maxn=max(maxn,ne);maxe=max(maxe,ee);maxend=max(maxend,endpoint)
 results.append({'case':idx,'A_J':A0,'q_J':q0,'cells':len(c),'max_number_error':ne,'max_spill_energy_error_over_A':ee,'max_endpoint_error_over_A':endpoint})
summary={'cases':len(cases),'digits':mp.mp.dps,'canonical_alpha_mass_kg':str(mass),'max_number_error':maxn,'max_spill_energy_error_over_A':maxe,'max_endpoint_error_over_A':maxend,'results':results}
Path(__file__).with_name('alpha0-reference.json').write_text(json.dumps(summary,indent=2)+'\n')
print({k:v for k,v in summary.items() if k!='results'})
