import ctypes as c, random, math, json, sys
import mpmath as mp
mp.mp.dps=100
lib=c.CDLL(sys.argv[1])
A6=c.c_double*6; A5=c.c_double*5
class Ledger(c.Structure):
 _fields_=[('events',A5),('loss',A6),('debit',A6),('birth',A6),('neutron',c.c_double),('nr',A6),('er',A6)]
f=lib.fusion_c_thermal_burn_trial
f.argtypes=[c.c_double,c.POINTER(c.c_double),c.POINTER(c.c_double),c.POINTER(c.c_double),c.POINTER(c.c_double),c.POINTER(c.c_double),c.POINTER(c.c_double),c.POINTER(c.c_double),c.POINTER(Ledger)]
a=[0,1,1,1,1];b=[5,1,1,2,3];mass=[1,2,3,3,4,11];charge=[1,1,1,2,2,5]
rng=random.Random(970314); worst_rate=mp.mpf(0);worst_balance=mp.mpf(0);worst_conservation=mp.mpf(0)
for case in range(500):
 # Independent raw returned-state BE residual audit over broad finite scales.
 old=[10**rng.uniform(-20,80) for _ in range(6)]
 k=[10**rng.uniform(-70,30) for _ in range(5)];dt=10**rng.uniform(-30,30)
 new=A6();en=A6();r=Ledger();code=f(dt,A6(*old),A6(),A5(*k),A5(),A5(),new,en,c.byref(r))
 if code: raise RuntimeError((case,code,old,k,dt))
 for j in range(5):
  expected=mp.mpf(dt)*k[j]*new[a[j]]*new[b[j]]/(2 if a[j]==b[j] else 1)
  err=abs(mp.mpf(r.events[j])-expected)/max(abs(expected),mp.mpf('1e-300'))
  worst_rate=max(worst_rate,err)
  assert err<mp.mpf('1e-11'),(case,j,float(err))
 for s in range(6):
  loss=sum(mp.mpf(r.events[j])*((a[j]==s)+(b[j]==s)) for j in range(5))
  err=abs(mp.mpf(old[s])-new[s]-loss)/max(mp.mpf(old[s]),loss)
  worst_balance=max(worst_balance,err);assert err<mp.mpf('2e-12'),(case,s,float(err))
 for weights,neutron_weight in [(mass,1),(charge,0)]:
  initial=sum(mp.mpf(old[s])*weights[s] for s in range(6))
  final=sum((mp.mpf(new[s])+r.birth[s])*weights[s] for s in range(6))+mp.mpf(r.neutron)*neutron_weight
  err=abs(initial-final)/initial
  worst_conservation=max(worst_conservation,err);assert err<mp.mpf('2e-12')
print(json.dumps(dict(cases=500,seed=970314,mpmath_digits=100,max_relative_BE_rate_residual=float(worst_rate),max_relative_species_balance=float(worst_balance),max_relative_baryon_charge_residual=float(worst_conservation)),indent=2))
