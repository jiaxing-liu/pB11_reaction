from pathlib import Path
import subprocess,json,hashlib,statistics,time
p=Path('/tmp/pb11-transfer-precompute-d885/coupled_benchmark')
def digest(q):return hashlib.sha256(q.read_bytes()).hexdigest()
sha=json.loads((p/'initial-sha.json').read_text())
for name in ['fp_bench_old','fp_bench_prepared']:sha[str(p/name)]=digest(p/name)
(p/'sha.json').write_text(json.dumps(sha,indent=2)+'\n')
runs=[]
def run(kind,reps,name):
 result=subprocess.run([str(p/('fp_bench_'+kind)),str(p/'edges.txt'),str(reps),str(p/(name+'.bin'))],check=True,capture_output=True,text=True)
 (p/(name+'.stdout')).write_text(result.stdout);(p/(name+'.stderr')).write_text(result.stderr)
 data=json.loads(result.stdout);data.update(kind=kind,output_sha256=digest(p/(name+'.bin')),output_bytes=(p/(name+'.bin')).stat().st_size)
 return data
base=run('old',1,'gate_old');new=run('prepared',1,'gate_prepared')
if (p/'gate_old.bin').read_bytes()!=(p/'gate_prepared.bin').read_bytes():raise RuntimeError('byte gate failed; no timing')
reference=(p/'gate_old.bin').read_bytes();print('1-trial byte gate PASS '+base['output_sha256'],flush=True)
for pair in range(6):
 for kind in (['old','prepared'] if pair%2==0 else ['prepared','old']):
  name=f'pair{pair}_{kind}';data=run(kind,1000,name)
  if (p/(name+'.bin')).read_bytes()!=reference:raise RuntimeError('timed output byte mismatch '+name)
  data['pair']=pair;runs.append(data);(p/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(json.dumps(data),flush=True)
pairs=[]
for pair in range(6):
 d={x['kind']:x['seconds'] for x in runs if x['pair']==pair};pairs.append({'pair':pair,**d,'speedup':d['old']/d['prepared'],'reduction':1-d['prepared']/d['old']})
summary={'scope':'same frozen coupled collision/handoff microbenchmark; not full BALDUR host speedup','all_outputs_byte_identical':True,'output_sha256':base['output_sha256'],'output_bytes':len(reference),'gate_runs':[base,new],'runs':runs,'pairs':pairs,'old_median_s':statistics.median(x['old'] for x in pairs),'prepared_median_s':statistics.median(x['prepared'] for x in pairs),'paired_speedup_median':statistics.median(x['speedup'] for x in pairs),'paired_reduction_median':statistics.median(x['reduction'] for x in pairs),'paired_speedup_range':[min(x['speedup'] for x in pairs),max(x['speedup'] for x in pairs)]}
(p/'benchmark.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY '+json.dumps({k:v for k,v in summary.items() if k not in ['runs','pairs','gate_runs']}),flush=True)
