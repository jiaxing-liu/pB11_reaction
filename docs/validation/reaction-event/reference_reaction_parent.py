from pathlib import Path
import subprocess,json,random,mpmath as m
m.mp.dps=80;c=m.mpf(299792458);rng=random.Random(314159)
cases=[]
for ch in range(5):
 for mode in range(2):
  for scale in [0,1e-25,1e-23,1e-21,1e-19]:cases.append((ch,mode,[rng.uniform(-1,1)*scale for _ in range(6)]))
text=''.join(f'{ch} {mode} '+' '.join(f'{p:.17g}' for p in ps)+'\n' for ch,mode,ps in cases)
r=Path(__file__).parent
p=subprocess.run([str(r/'reaction_parent_reference_driver')],input=text,text=True,capture_output=True,check=True)
lines=p.stdout.splitlines();assert len(lines)==len(cases)
records=[]
for (ch,mode,ps),line in zip(cases,lines):
 v=list(map(float,line.split()));assert v[0]==0,(ch,mode,ps,v[0])
 Q=m.mpf(v[1]);ma,mb=map(m.mpf,v[2:4]);rest=sum(map(m.mpf,v[4:7]))*c*c
 pa=list(map(m.mpf,ps[:3]));pb=list(map(m.mpf,ps[3:]));P=[a+b for a,b in zip(pa,pb)]
 nr=[sum(x*x for x in p)/(2*mass) for p,mass in [(pa,ma),(pb,mb)]]
 on=[m.sqrt((mass*c*c)**2+sum(x*x for x in p)*c*c)-mass*c*c for p,mass in [(pa,ma),(pb,mb)]]
 T=on if mode else nr;lab=Q+sum(T);Etot=rest+lab
 # Deliberately direct invariant subtraction at80 digits, unlike stable production excess formula.
 A=m.sqrt(Etot*Etot-sum(x*x for x in P)*c*c)-rest
 mu=ma*mb/(ma+mb);Erel=mu/2*sum((a/ma-b/mb)**2 for a,b in zip(pa,pb));diff=sum(nr)-sum(on)
 expected=[A,lab,Erel,diff]+T+[c*c*x/Etot for x in P]
 actual=v[7:]
 errors=[]
 for i,(a,b) in enumerate(zip(actual,expected)):
  den=max(abs(b),abs(Q)*m.mpf('1e-40')) if i<6 else max(abs(b),m.mpf('1e-40'))
  err=float(abs(m.mpf(a)-b)/den);errors.append(err)
 assert max(errors)<1e-10,(ch,mode,ps,errors)
 records.append({'channel':ch,'mode':mode,'momenta':ps,'max_relative_error':max(errors),'errors':errors})
summary={'cases':len(records),'digits':m.mp.dps,'max_relative_error':max(x['max_relative_error'] for x in records),'records':records}
(r/'reaction-parent-reference.json').write_text(json.dumps(summary,indent=2)+'\n');print('checked',summary['cases'],'max relative',summary['max_relative_error'])
