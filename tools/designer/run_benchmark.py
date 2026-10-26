#!/usr/bin/env python3
"""Run every frozen case in an isolated measured R process; fail on a lost gate."""
from pathlib import Path
import csv,datetime,hashlib,json,os,re,signal,subprocess,time
root=Path(__file__).resolve().parents[2]
dest=root/'docs/evidence/designer-v1';dest.mkdir(parents=True,exist_ok=True)
freeze=json.loads((dest/'freeze-manifest.json').read_text())
assert all(hashlib.sha256((root/f).read_bytes()).hexdigest()==h for f,h in freeze['files'].items())
policy={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'rss_limit_bytes':2*1024**3,'process_wall_limit_seconds':60,'sample_seconds':.25,
        'note':'Supplemental process cap declared before comparative benchmark. Frozen fixtures/objectives/tolerances unchanged. RSS includes R startup, optimizer, ranked comparator and exhaustive oracle.'}
(dest/'run-policy.json').write_text(json.dumps(policy,indent=2)+'\n')
cases=['knapsack','set_knapsack','ties','small_set','per_drawing_infeasible','smaller_budget','exhausted']
env=os.environ.copy();env.update(LC_ALL='en_US.UTF-8',LANG='en_US.UTF-8')
rows=[]
for name in cases:
    log=dest/(name+'.log');started=time.monotonic();peak=0;reason='exited'
    args=['/usr/bin/time','-l','Rscript','--vanilla','tools/designer/run_benchmark.R',name,str(dest)]
    with log.open('wb') as stream:
        proc=subprocess.Popen(args,cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            try:
                table=[list(map(int,line.split())) for line in subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True).splitlines()]
            except Exception:
                if proc.poll() is None:os.killpg(proc.pid,signal.SIGTERM)
                proc.wait(timeout=5)
                raise
            owned={proc.pid}
            for _ in range(8):owned.update(pid for pid,parent,rss in table if parent in owned)
            usage=sum(rss*1024 for pid,parent,rss in table if pid in owned);peak=max(peak,usage)
            if usage>policy['rss_limit_bytes'] or time.monotonic()-started>policy['process_wall_limit_seconds']:
                reason='rss_limit' if usage>policy['rss_limit_bytes'] else 'wall_limit'
                os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=5);break
            time.sleep(.25)
    match=re.search(r'(\d+)\s+maximum resident set size',log.read_text())
    receipt={'case':name,'args':args,'exit_code':proc.returncode,'reason':reason,
             'process_seconds':time.monotonic()-started,'sampled_tree_rss_bytes':peak,
             'maximum_resident_set_size_bytes':int(match.group(1)) if match else None}
    (dest/(name+'-process.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    if proc.returncode or reason!='exited':raise SystemExit(f'{name}: {receipt}')
    with (dest/(name+'.csv')).open() as inp:row=next(csv.DictReader(inp))
    row.update({k:receipt[k] for k in ['process_seconds','sampled_tree_rss_bytes','maximum_resident_set_size_bytes']})
    rows.append(row);print(name,row['status'],row['objective'],flush=True)
with (dest/'results.csv').open('w') as out:
    writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
