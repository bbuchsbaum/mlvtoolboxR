#!/usr/bin/env python3
"""Build/install/check this package with an owned 200 MiB work-directory cap."""
from pathlib import Path
import os,subprocess,signal,time,json,sys,tarfile
root=Path(__file__).resolve().parents[1]
work=root/'build';work.mkdir(exist_ok=True)
tmp=work/'tmp';tmp.mkdir(exist_ok=True)
lib=work/'library';lib.mkdir(exist_ok=True)
logs=root/'.agent-work';logs.mkdir(exist_ok=True)
limit=200*1024*1024
env=os.environ.copy()
env.update({'LC_ALL':'en_US.UTF-8','LANG':'en_US.UTF-8','TMPDIR':str(tmp),
            'R_PROFILE_USER':'/dev/null','R_ENVIRON_USER':'/dev/null'})
commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
def usage():
    return sum(p.stat().st_size for p in work.rglob('*') if p.is_file())
def run(stage,args,timeout=300):
    log=logs/f'{stage}.log'
    if log.exists():raise RuntimeError(f'Refusing to overwrite {log}')
    started=time.monotonic();peak=usage();reason='exited'
    with log.open('wb') as stream:
        proc=subprocess.Popen(args,cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            peak=max(peak,usage())
            if peak>limit or time.monotonic()-started>timeout:
                reason='work_directory_limit' if peak>limit else 'timeout'
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                break
            time.sleep(.25)
    meta={'stage':stage,'args':args,'exit_code':proc.returncode,'reason':reason,
          'duration_seconds':time.monotonic()-started,'peak_work_bytes':peak,
          'limit_bytes':limit,'source_commit':commit,'log':str(log)}
    (logs/f'{stage}.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta),flush=True)
    if proc.returncode or reason!='exited':raise SystemExit(proc.returncode or 125)
run('build-01',['R','CMD','build','--no-manual',str(root)])
tarball=work/'mlvtoolbox_0.1.0.tar.gz'
with tarfile.open(tarball) as archive:
    names=archive.getnames()
    assert not any('/vendor/' in x or '/reference-data/' in x or '/.git/' in x or '/.mote/' in x for x in names)
    (logs/'tarball-contents.txt').write_text('\n'.join(names)+'\n')
run('install-01',['R','CMD','INSTALL','--preclean','--library='+str(lib),str(tarball)])
run('check-01',['R','CMD','check','--as-cran','--no-manual',str(tarball)],timeout=420)
print('Bounded package check completed; log: build/mlvtoolbox.Rcheck/00check.log')
