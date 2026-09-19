#!/usr/bin/env python3
"""Independent hashlib oracle; run the C++ fixture in an isolated directory."""
import hashlib,json,pathlib,random,subprocess,sys
exe=pathlib.Path(sys.argv[1]).resolve();out=pathlib.Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
p=subprocess.run([str(exe)],cwd=out,text=True,capture_output=True,check=True)
(out/'run.log').write_text(p.stdout+p.stderr)
results=[]
for line in p.stdout.splitlines():
 if line.startswith(('thermal ','beam ')):
  kind,n,digest=line.split();b=(out/f'{kind}-digest.bin').read_bytes()
  assert len(b)==int(n) and hashlib.sha256(b).hexdigest()==digest
  results.append({'table':kind,'bytes':len(b),'sha256':digest})
rng=random.Random(188)
for n in [0,1,55,56,63,64,65,119,120,127,128,129,4096,65537,1000000]:
 b=rng.randbytes(n);f=out/'hash-input.bin';f.write_bytes(b)
 got=subprocess.check_output([str(exe),'--hash',str(f)],text=True).strip()
 assert got==hashlib.sha256(b).hexdigest(),n
 results.append({'input_bytes':n,'sha256':got})
(out/'hashlib-results.json').write_text(json.dumps(results,indent=2)+'\n')
print('HASHLIB_ORACLE_PASS',len(results))
