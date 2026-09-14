from pathlib import Path
import csv,json,hashlib,shutil
root=Path('/home/cloud/research/pB11_reaction')
out=root/'docs/validation/birth-table/companion-fuels';out.mkdir(parents=True,exist_ok=True)
summary={}
for fuel in ('dd','dhe3'):
    names={mode:f'coupled-{fuel}-s16-n800{suffix}' for mode,suffix in [('direct','-table-era-direct'),('table','-table001'),('previous_direct','')]}
    data={}
    for mode,name in names.items():
        for ext in ('csv','log'): shutil.copyfile('/tmp/'+name+'.'+ext,out/(name+'.'+ext))
        with (out/(name+'.csv')).open() as h: data[mode]=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(h)]
    a,b,c=(data[k] for k in ('direct','table','previous_direct'))
    assert len(a)==len(b)==len(c)==16 and a[-1]['time_s']==b[-1]['time_s']==.01
    assert all(x['time_s']==y['time_s'] for x,y in zip(a,b))
    parity=all(x[k]==z[k] for x,z in zip(a,c) for k in x.keys()&z.keys());assert parity
    fields={}
    for k in a[0]:
        if k in ('step','time_s'):continue
        diffs=[y[k]-x[k] for x,y in zip(a,b)]
        rel=[abs(d/x[k]) for x,d in zip(a,diffs) if x[k]!=0]
        fields[k]={'direct_final':a[-1][k],'table_final':b[-1][k],'final_difference':diffs[-1],'final_relative_percent':100*diffs[-1]/a[-1][k] if a[-1][k] else None,'max_history_absolute_difference':max(map(abs,diffs)),'max_history_relative_percent_nonzero_reference':100*max(rel) if rel else None,'zero_reference_nonzero_difference_count':sum(x[k]==0 and d!=0 for x,d in zip(a,diffs))}
    summary[fuel]={'rows':16,'duration_s':.01,'old_new_direct_common_columns_exact':parity,'max_table_energy_relative_residual':max(abs(x['energy_relative_residual']) for x in b),'fields':fields}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
shutil.copyfile(__file__,out/'archive_and_compare.py')
manifest={'binary':'/tmp/study_coupled_thermal_table','binary_sha256':hashlib.sha256(Path('/tmp/study_coupled_thermal_table').read_bytes()).hexdigest(),'commands':[f'/tmp/study_coupled_thermal_table 16 800 {f} 0.01 1 8'+s for f in ('dd','dhe3') for s in ('',' 0.001 19.5 23')],'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.name!='manifest.json'}}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(out/'README.md').write_text('''# DD and D–He3 companion source-table comparisons

These are local numerical stress tests, not EXL predictions or a fair fuel comparison. Both use 800 energy cells, 16 steps over 10 ms, the existing common-temperature birth model at nq=ncos8, and table bounds 19.5–23 keV with 0.1% sampled interpolation gates. The direct and table paths have otherwise identical inputs. Exact commands, binary hash, full CSV/logs and per-column full-history differences are retained. No thermal handoff projections occur.

DD evaluates both primary branches (channels1/2); the driver also constructs secondary DT/DHe3 tables for its enabled network. Final table/direct differences are fast energy0.01159%, electron heat0.01051%, network-ion heat0.01144%, carbon heat0.01833% (0.0109775 J/m3), and Q0.01050%. Maximum local energy residual is1.9423e-16.

D–He3 final differences are fast energy0.01744%, electron heat0.01613%, network-ion heat0.01726%, carbon heat0.01829% (0.0129664 J/m3), and Q0.01624%. Maximum local energy residual is1.3629e-16. Full-history electron-heat discrepancy peaks at0.01622%; terminal comparisons do not replace history checks.

The previous direct files predate the additional separate heat columns. All16 rows and every shared column agree exactly with the newly generated direct files. This checks the source-provider refactor independently of interpolated-mode sensitivity. Zero reference values are reported explicitly rather than assigned an arbitrary relative-error floor. These runs do not establish timestep/grid/source-model convergence; those remain separate checks.
''')
