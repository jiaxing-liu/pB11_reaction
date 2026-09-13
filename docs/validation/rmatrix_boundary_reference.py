import json
import mpmath as mp
mp.mp.dps=40
# Explicit numerical conventions; source-independent reproducibility inputs.
a=mp.mpf('4.5'); mu=mp.mpf('3727.3794118')/2
alpha=1/mp.mpf('137.035999084'); hbarc=mp.mpf('197.3269804')
def values(E):
 eta=4*alpha*mp.sqrt(mu/(2*E)); rho=a*mp.sqrt(2*mu*E)/hbarc
 F=mp.coulombf(2,eta,rho); G=mp.coulombg(2,eta,rho)
 Fp=mp.diff(lambda r:mp.coulombf(2,eta,r),rho)
 Gp=mp.diff(lambda r:mp.coulombg(2,eta,r),rho)
 den=F*F+G*G
 return dict(eta=eta,rho=rho,P=rho/den,S=rho*(F*Fp+G*Gp)/den,wronskian=Fp*G-F*Gp)
r=values(mp.mpf('3.129'))
assert abs(r['wronskian']-1)<mp.mpf('1e-30')
r['formal_width_MeV']=2*mp.mpf('1.075')*r['P']
print(json.dumps({'model':'Laursen2016 Be8 2+ boundary S(E0), not complete spectrum','E0_MeV':'3.129','radius_fm':str(a),'mu_MeV_c2':str(mu),'alpha_inverse':'137.035999084','hbarc_MeV_fm':str(hbarc),'dps':mp.mp.dps,'values':{k:mp.nstr(v,30) for k,v in r.items()}},indent=2))
E=mp.mpf('3.129'); g2=mp.mpf('1.075')
for h in [mp.mpf('.0001'),mp.mpf('.00001')]:
 derivative=(values(E+h)['S']-values(E-h)['S'])/(2*h)
 print('step_MeV',h,'dS_dE',mp.nstr(derivative,24),'linearized_width_MeV',mp.nstr(2*g2*r['P']/(1+g2*derivative),24))
