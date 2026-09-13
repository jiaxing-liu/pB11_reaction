#!/usr/bin/env python3
"""Compare a new-API grid to the frozen independent energy-space study."""
import csv,json,sys
A=list(csv.DictReader(open(sys.argv[1])))
B=list(csv.DictReader(open(sys.argv[2])))
assert len(A)==len(B)==65
pairs={'K_fit':'K_fit_m3_s','K_low':'K_low_Sconst_m3_s','K_high':'K_high_Sconst_m3_s','E_fit':'relative_E_fit_keV','E_low':'relative_E_low_keV','E_high':'relative_E_high_keV'}
errors={}
for k,v in pairs.items():
    worst=(0,None)
    for a,b in zip(A,B):
        assert a['channel']==b['channel'] and float(a['T_keV'])==float(b['T_keV'])
        x,y=float(a[k]),float(b[v])
        error=abs(x-y)/max(abs(x),abs(y),1e-300)
        if error>worst[0]: worst=(error,(a['channel'],a['T_keV'],x,y))
    errors[k]=worst
print(json.dumps(errors,indent=2))
assert max(v[0] for v in errors.values())<2e-7
