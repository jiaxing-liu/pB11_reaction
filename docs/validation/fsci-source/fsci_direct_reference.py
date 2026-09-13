import csv
import mpmath as mp
mp.mp.dps=50
m=mp.mpf('3727.3794118');alpha=1/mp.mpf('137.035999084');hbarc=mp.mpf('197.3269804')
def log_inverse_penetration(l,E,primary):
 mu=2*m/3 if primary else m/2
 z=8 if primary else 4
 eta=z*alpha*mp.sqrt(mu/(2*E));rho=16*mp.sqrt(2*mu*E)/hbarc
 f=mp.coulombf(l,eta,rho);g=mp.coulombg(l,eta,rho)
 return mp.log(f*f+g*g)
cases=[(1,'8.84','3.129','0'),(2,'8.84','3.129','.75'),(3,'9.3','3.129','-.3'),(2,'8.84','.1','.8'),(2,'8.84','8.6','-.4'),(2,'8.84','4.42','0')]
with open('/tmp/fsci-direct-reference.csv','w',newline='') as stream:
 writer=csv.writer(stream,lineterminator='\n');writer.writerow(['l','A_MeV','q_MeV','cosine','E12_MeV','E13_MeV','log_correction','amplitude_factor'])
 for l,aa,qq,cc in cases:
  A,q,c=map(mp.mpf,(aa,qq,cc));shift=mp.sqrt(3*q*(A-q))*c/2
  e12=3*A/4-q/2+shift;e13=3*A/4-q/2-shift
  logC=log_inverse_penetration(l,A-q,True)-log_inverse_penetration(2,e12,False)-log_inverse_penetration(2,e13,False)
  writer.writerow([l,aa,qq,cc,*[mp.nstr(v,40) for v in [e12,e13,logC,mp.exp(logC/2)]]])
  stream.flush()
