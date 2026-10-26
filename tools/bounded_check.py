#!/usr/bin/env python3
"""Build/install/check this package with an owned 200 MiB work-directory cap."""
from pathlib import Path
import os,subprocess,signal,time,json,sys,tarfile,shutil,hashlib
root=Path(__file__).resolve().parents[1]
work=root/'.agent-work/package-build';work.mkdir(parents=True,exist_ok=True)
tmp=work/'tmp';tmp.mkdir(exist_ok=True)
lib=work/'library';lib.mkdir(exist_ok=True)
logs=root/'.agent-work';logs.mkdir(exist_ok=True)
limit=200*1024*1024
env=os.environ.copy()
env.update({'LC_ALL':'en_US.UTF-8','LANG':'en_US.UTF-8','TMPDIR':str(tmp),
            'R_PROFILE_USER':'/dev/null','R_ENVIRON_USER':'/dev/null'})
commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
staged=work/'package-source';staged.mkdir(exist_ok=True)
tracked=subprocess.check_output(['git','ls-files','-z']).decode().split(chr(0))
manifest={}
for rel in tracked:
    if not rel: continue
    path=Path(rel)
    if path.parts[0] not in {'R','src','man','inst','tests','vignettes'} and rel not in {'DESCRIPTION','NAMESPACE','README.md','.Rbuildignore','LICENSE','NEWS.md'}:continue
    original=root/path
    if not original.is_file():continue
    dest=staged/path;dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(original,dest)
    manifest[rel]=hashlib.sha256(original.read_bytes()).hexdigest()
(logs/'package-source-manifest.json').write_text(json.dumps(manifest,indent=2)+chr(10))
def usage():
    return sum(p.stat().st_size for p in work.rglob('*') if p.is_file())
def run(stage,args,timeout=300):
    stem=stage.rsplit('-',1)[0];attempt=1
    while (logs/f'{stem}-{attempt:02d}.log').exists():attempt+=1
    stage=f'{stem}-{attempt:02d}';log=logs/f'{stage}.log'
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
if '--coverage-only' in sys.argv:
    run('coverage-01',['Rscript','--vanilla',str(root/'tools/coverage.R'),str(staged),str(root/'docs/evidence')],timeout=240)
    raise SystemExit(0)
run('build-02',['R','CMD','build','--no-manual',str(staged)])
tarball=work/'mlvtoolbox_0.1.0.tar.gz'
with tarfile.open(tarball) as archive:
    names=archive.getnames()
    assert not any('/vendor/' in x or '/reference-data/' in x or '/.git/' in x or '/.mote/' in x for x in names)
    (logs/'tarball-contents.txt').write_text('\n'.join(names)+'\n')
run('install-01',['R','CMD','INSTALL','--preclean','--library='+str(lib),str(tarball)])
run('check-01',['R','CMD','check','--as-cran','--no-manual',str(tarball)],timeout=420)
print('Bounded package check completed; log: '+str(work/'mlvtoolbox.Rcheck/00check.log'))
