import mpmath as m,csv,json
from pathlib import Path
m.mp.dps=40
r=Path(__file__).parent
rows=list(csv.DictReader((r/'pb-population-study.csv').open()))
# Canonical rounded masses, independently reconstructed from source constants.
u=m.mpf('1.602176634e-16');barn=m.mpf('1e-28')
unit=m.mpf('1.66053906892e-27');c=m.mpf(299792458)
mp=m.mpf(float(m.mpf('1.0072764665789')*unit))
mb=m.mpf(float((m.mpf('11.009305166')-5*m.mpf('5.485799090441e-4'))*unit+m.mpf('670.9838405')*m.mpf('1.602176634e-19')/c**2))
M=mp+mb;mu=mp*mb/M;ratio=mb/M
energies=list(map(m.mpf,['247','355','461','565','668','771','873','975','1077','1178','1279','1380','1481','1582','1683','1783','1884','1985','2085','2186']))
a0=list(map(m.mpf,['.64','1.40','2.02','2.53','3.31','3.45','3.04','2.90','2.61','2.27','2.07','1.66','1.27','.97','1.35','1.96','2.43','8.08','9','4.09']))
a1=list(map(m.mpf,['40','148','357','598','668','386','234','171','152','156','161','160','146','129','117','110','102','102','84','84']))
values=[a/(a+b) for a,b in zip(a0,a1)]
def weights(E):
 if E<=400:
  pole=m.mpf('1.82e4')/((E-148)**2+m.mpf('2.35')**2)
  narrow=pole/(197+m.mpf('.269')*E+m.mpf('2.54e-4')*E*E+pole)
 else:narrow=m.mpf(0)
 lab=E/ratio
 if lab<=energies[0]:z=values[0]
 elif lab>=energies[-1]:z=values[-1]
 else:
  j=next(j for j in range(1,20) if lab<energies[j]);f=(lab-energies[j-1])/(energies[j]-energies[j-1]);z=values[j-1]*(1-f)+values[j]*f
 return [m.mpf(1),narrow*m.mpf('.051')+(1-narrow)*z,narrow*m.mpf('.949'),(1-narrow)*(1-z),(1-narrow) if lab<energies[0] or lab>energies[-1] else m.mpf(0)]
def S(E):
 if E<=400:return 197+m.mpf('.269')*E+m.mpf('2.54e-4')*E*E+m.mpf('1.82e4')/((E-148)**2+m.mpf('2.35')**2)
 if E<=668:
  x=(E-400)/100;return 346+150*x-m.mpf('59.9')*x*x-m.mpf('.460')*x**5
 return m.mpf('.381')+sum(A/((E-e)**2+g*g) for A,e,g in zip(map(m.mpf,['1.98e6','3.89e6','1.36e6','3.71e6']),map(m.mpf,['640.9','1211','2340','3294']),map(m.mpf,['85.5','414','221','351'])))
def sigma_barn(E):
 if E==0:return m.mpf(0)
 return 1000*S(min(E,m.mpf(9760)))/E*m.exp(-m.sqrt(22589/E))
basecuts=[m.mpf(x) for x in [0,22,101,148,195,400,668,1211,2340,3294,5700,9760]]+[e*ratio for e in energies]
results=[]
# Six unequal-temperature thermal cases; two beam cases cover narrow and broad ranges.
for kind,T in [(0,x) for x in [2,10,50,150,300,1000]]+[(1,150),(1,1000)]:
 ta=m.mpf(T);tb=ta*m.mpf('.7') if kind==0 else m.mpf(3)
 cuts=sorted(set(basecuts+[ta*x for x in [1,4,16,64,1600]]))
 if kind==1:
  v=m.sqrt(2*ta*u/mp);speed=m.sqrt(2*tb*u/mb)
  cuts=sorted(set(cuts+[mu*(max(m.mpf(0),v+speed*d))**2/(2*u) for d in [-40,-16,-8,-4,-2,-1,0,1,2,4,8,16,40]]))
 else:tr=(mb*ta+mp*tb)/M
 # 1D energy integrals are independent of the production speed-space GK61 kernel.
 for part in range(5):
  if kind==0:
   def fn(E,power):return weights(E)[part]*sigma_barn(E)*E**power*m.exp(-E/tr)
   I=[m.quad(lambda e:fn(e,power),cuts+[m.inf]) for power in [1,2]]
   K=m.sqrt(8*u/(m.pi*mu))/tr**m.mpf('1.5')*barn*I[0]
   Er=m.sqrt(8*u/(m.pi*mu))/tr**m.mpf('1.5')*barn*u*I[1]
   cmrandom=m.mpf('1.5')*ta*tb/tr*u*K
   Ea=mp/M*cmrandom+mb/M*(ta/tr)**2*Er
   Eb=mb/M*cmrandom+mp/M*(tb/tr)**2*Er
   Ecm=cmrandom+mu/M*((ta-tb)/tr)**2*Er
  else:
   shift=((m.sqrt(2*400*u/mu)-v)/speed)**2 if part==2 and mu*v*v/(2*u)>400 else m.mpf(0)
   def fn(E,moment):
    if E==0:return m.mpf(0)
    w=m.sqrt(2*E*u/mu);a=2*v*w/speed**2
    diff=m.exp(-((w-v)/speed)**2)*(-m.expm1(-2*a))
    prob=diff*u/(m.sqrt(m.pi)*speed*v*mu)
    angular=(1/m.tanh(a)-1/a) if a else m.mpf(0)
    Ebcond=mb*(v*v+w*w-2*v*w*angular)/(2*u)
    return m.exp(shift)*weights(E)[part]*sigma_barn(E)*w*prob*(1 if moment==0 else E if moment==1 else Ebcond)
   I=[m.quad(lambda e:fn(e,k),cuts+[m.inf]) for k in [0,1,2]]
   K=barn*I[0]*m.exp(-shift);Er=barn*u*I[1]*m.exp(-shift);Eb=barn*u*I[2]*m.exp(-shift);Ea=ta*u*K;Ecm=Ea+Eb-Er
  row=next(x for x in rows if int(x['kind'])==kind and float(x['Ta_or_beam_keV'])==T and int(x['component'])==part)
  errs={name:float(abs(m.mpf(row[name])-expected)/max(abs(expected),m.mpf('1e-300'))) for name,expected in zip(['K','Ea',' Eb','Erel','Ecm'],[K,Ea,Eb,Er,Ecm])}
  assert max(errs.values())<2e-7,(kind,T,part,errs)
  results.append({'kind':kind,'T_keV':T,'component':part,'errors':errs})
 print('checked',kind,T,flush=True)
summary={'digits':m.mp.dps,'component_cases':len(results),'max_relative_error':max(e for z in results for e in z['errors'].values()),'results':results}
(r/'pb-population-reference.json').write_text(json.dumps(summary,indent=2)+'\n');print('max error',summary['max_relative_error'])
