#!/usr/bin/env python3
"""Bounded benchmark of unmodified upstream AOF with explicit synthetic inputs."""
from pathlib import Path
import subprocess,tempfile,csv,json,platform
root=Path(__file__).resolve().parents[1]
out=root/'docs/evidence/benchmark-octave.csv'
rows=[]
for n in [32,64]:
 with tempfile.TemporaryDirectory(prefix='mlv-bench-') as td:
  dest=Path(td)/'times.csv'
  expr=f"""addpath('{root/'vendor/MLV_toolbox/MLVcode'}');
[rr,cc]=ndgrid(1:{n},1:{n});r=round({n}/2);c=round({n}/2);
D=sqrt((rr-r).^2+(cc-c).^2);IDX=repmat(sub2ind(size(D),r,c),size(D));
alpha=2*pi/60;s=transpose(0:59);circle=[cos(alpha*s),sin(alpha*s)];
tic;F=computeAOF(D,IDX,circle,1);cold=toc;
t=zeros(1,3);for k=1:3,tic;F=computeAOF(D,IDX,circle,1);t(k)=toc;end
dlmwrite('{dest}',[cold,t],',','precision',17);
assert(all(isfinite(F(:))));
dlmwrite('{root/'docs/evidence'}/aof-octave-{n}.csv',F,',','precision',17);
img=ones({n},{n});img(3:end-2,round({n}*.25))=0;img(3:end-2,round({n}*.75))=0;

"""
  try:
   p=subprocess.run(['octave-cli','--quiet','--no-history','--eval',expr],capture_output=True,text=True,timeout=60)
   if p.returncode:raise RuntimeError(p.stderr)
   times=[float(s) for s in dest.read_text().strip().split(',')]
   rows.append({'stage':'aof','size':n,'cold_seconds':times[0],'median_seconds':sorted(times[1:])[1],
                'min_seconds':min(times[1:]),'max_seconds':max(times[1:]),'status':'completed',
                'runtime':'Octave 11.1.0','source':'f12418a5e03c2f38e222fed1dffccf42b227cd29'})
  except subprocess.TimeoutExpired:
   rows.append({'stage':'aof','size':n,'status':'timeout_60s','runtime':'Octave 11.1.0',
                'source':'f12418a5e03c2f38e222fed1dffccf42b227cd29'})
out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w') as f:
 writer=csv.DictWriter(f,fieldnames=['stage','size','cold_seconds','median_seconds','min_seconds','max_seconds','status','runtime','source'])
 writer.writeheader();writer.writerows(rows)
print(json.dumps(rows,indent=2))
