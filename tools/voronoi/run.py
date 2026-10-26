#!/usr/bin/env python3
"""Bounded research driver. Owns only its prototype process groups and paths."""
from pathlib import Path
import datetime,json,os,signal,subprocess,time
root=Path(__file__).resolve().parents[2];dest=root/'docs/evidence/voronoi-spike'
work=root/'.agent-work/voronoi-spike';work.mkdir(exist_ok=True)
preflight=json.loads((dest/'preflight.json').read_text());box=json.loads((dest/'timebox.json').read_text())
env=os.environ.copy();env.update(LC_ALL='en_US.UTF-8',LANG='en_US.UTF-8',R_LIBS=str(root/'.agent-work/graph-build/library'))
def run(name,args,limit):
    remaining=box['start_epoch']+box['limit_seconds']-time.time()
    if remaining<=0:raise SystemExit('four-hour prototype timebox exhausted')
    limit=min(limit,remaining);started=time.monotonic();peak=0;reason='exited'
    log=dest/(name+'.log')
    with log.open('wb') as stream:
        proc=subprocess.Popen(['/usr/bin/time','-l',*args],cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            while proc.poll() is None:
                rows=[list(map(int,line.split())) for line in subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True).splitlines()]
                owned={proc.pid}
                for _ in range(8):owned.update(pid for pid,parent,rss in rows if parent in owned)
                memory=sum(rss*1024 for pid,parent,rss in rows if pid in owned);peak=max(peak,memory)
                usage=sum(p.stat().st_size for base in [work,dest] for p in base.rglob('*') if p.is_file())
                if memory>2*1024**3 or usage>200*1024**2 or time.monotonic()-started>limit:
                    reason='resource_limit';os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=5);break
                time.sleep(.25)
        except Exception:
            if proc.poll() is None:os.killpg(proc.pid,signal.SIGTERM)
            proc.wait(timeout=5);raise
    receipt={'phase':name,'args':args,'exit_code':proc.returncode,'reason':reason,
             'seconds':time.monotonic()-started,'sampled_tree_peak_rss_bytes':peak,
             'rss_cap_bytes':2*1024**3,'workspace_cap_bytes':200*1024**2}
    (dest/(name+'.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
    if proc.returncode or reason!='exited':raise SystemExit(proc.returncode or 125)
executable=work/'prototype'
if executable.exists():executable.unlink()
run('compile',['clang++','-std=c++17','-O2','-I'+preflight['BH_include'],'tools/voronoi/prototype.cpp','-o',str(executable)],box['compiler_seconds_limit'])
run('experiment',['Rscript','--vanilla','tools/voronoi/experiment.R'],box['experiment_seconds_limit'])
box['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();box['elapsed_seconds']=time.time()-box['start_epoch']
(dest/'timebox-result.json').write_text(json.dumps(box,indent=2)+'\n')
