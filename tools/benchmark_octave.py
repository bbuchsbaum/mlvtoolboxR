#!/usr/bin/env python3
"""Bounded upstream timing and output witnesses for common AOF/MAT work."""
from pathlib import Path
import subprocess,tempfile,csv,json
root=Path(__file__).resolve().parents[1]
out=root/'docs/evidence/benchmark-octave.csv'
rows=[]
source=root/'vendor/MLV_toolbox/MLVcode'
for n in [32,64]:
 for stage in ['aof','mat_core']:
  with tempfile.TemporaryDirectory(prefix='mlv-bench-') as td:
   dest=Path(td)/'times.csv'
   text=(source/'computeGradientVectorField.m').read_text()
   assert text.count('jj++m_Neighbors8')==2
   (Path(td)/'computeGradientVectorField.m').write_text(text.replace('jj++m_Neighbors8','jj + m_Neighbors8'))
   setup=f"addpath('{source}');addpath('{td}','-begin');"
   if stage=='aof':
    setup+=f"""[rr,cc]=ndgrid(1:{n},1:{n});r=round({n}/2);c=round({n}/2);
D=sqrt((rr-r).^2+(cc-c).^2);IDX=repmat(sub2ind(size(D),r,c),size(D));
alpha=2*pi/60;s=transpose(0:59);circle=[cos(alpha*s),sin(alpha*s)];"""
    call="F=computeAOF(D,IDX,circle,1);"
    save=f"dlmwrite('{root/'docs/evidence'}/aof-octave-{n}.csv',F,',','precision',17);"
   else:
    setup+=f"img=ones({n},{n});img(3:end-2,round({n}*.25))=0;img(3:end-2,round({n}*.75))=0;"
    call="MAT=computeMAT(img,28);"
    save=f"""dlmwrite('{root/'docs/evidence'}/mat-aof-octave-{n}.csv',MAT.AOF,',','precision',17);
dlmwrite('{root/'docs/evidence'}/mat-distance-octave-{n}.csv',MAT.distance_map,',','precision',17);
dlmwrite('{root/'docs/evidence'}/mat-skeleton-octave-{n}.csv',MAT.skeleton,',','precision',17);"""
   expr=setup+f"""tic;{call}cold=toc;
t=zeros(1,3);for k=1:3,tic;{call}t(k)=toc;end
dlmwrite('{dest}',[cold,t],',','precision',17);"""+save
   try:
    p=subprocess.run(['octave-cli','--quiet','--no-history','--eval',expr],capture_output=True,text=True,timeout=60)
    if p.returncode:raise RuntimeError(p.stderr)
    times=[float(s) for s in dest.read_text().strip().split(',')]
    rows.append({'stage':stage,'size':n,'cold_seconds':times[0],'median_seconds':sorted(times[1:])[1],
                 'min_seconds':min(times[1:]),'max_seconds':max(times[1:]),'status':'completed',
                 'runtime':'Octave 11.1.0; image 2.20.0','source':'f12418a5e03c2f38e222fed1dffccf42b227cd29'})
   except subprocess.TimeoutExpired:
    rows.append({'stage':stage,'size':n,'status':'timeout_60s','runtime':'Octave 11.1.0; image 2.20.0',
                 'source':'f12418a5e03c2f38e222fed1dffccf42b227cd29'})
   print(json.dumps(rows[-1]),flush=True)
out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w') as f:
 writer=csv.DictWriter(f,fieldnames=['stage','size','cold_seconds','median_seconds','min_seconds','max_seconds','status','runtime','source'])
 writer.writeheader();writer.writerows(rows)
