#!/usr/bin/env python3
"""Optional local parity/timing witness using the retained upstream forest demo."""
from pathlib import Path
import json,hashlib,subprocess,sys,time,statistics,copy,csv
import numpy as np
root=Path(__file__).resolve().parents[2];out=root/'.agent-work/parity-perf/natural'
out.mkdir(parents=True,exist_ok=True)
manifest=json.loads((root/'vendor/SOURCES.json').read_text())
entry=next(x for x in manifest['shared_data']['files'] if x['source_path']=='demos/Forest_DBW.mat')
asset=Path(entry['destination'])
assert hashlib.sha256(asset.read_bytes()).hexdigest()==entry['sha256']
source=root/'vendor/MLV_toolbox/MLVcode'
expr=f"""addpath('{source}');load('{asset}','vecLD');
originalContours=vecLD.numContours;removed=0;cells={{}};ids=[];
for c=1:vecLD.numContours
 z=vecLD.contours{{c}};keep=any(z(:,1:2)~=z(:,3:4),2);removed=removed+sum(~keep);
 if any(keep),cells{{end+1}}=z(keep,:);ids(end+1)=c;end
end
v=struct('contours',{{cells}},'numContours',numel(cells),'imsize',vecLD.imsize);
Z=[];for c=1:v.numContours,Z=[Z;repmat(c,size(cells{{c}},1),1),cells{{c}}];end
dlmwrite('{out/'input.csv'}',Z,',','precision',17);
dlmwrite('{out/'dimensions.csv'}',v.imsize,',','precision',17);
dlmwrite('{out/'cleaning.csv'}',[originalContours,v.numContours,size(Z,1),removed],',','precision',17);
t=zeros(1,4);for k=1:4,tic;g=computeContourProperties(v,{{'length','orientation','curvature'}});t(k)=toc;end
dlmwrite('{out/'octave-geometry-time.csv'}',t,',','precision',17);
G=[];for c=1:g.numContours,G=[G;g.lengths{{c}}(:),g.orientations{{c}}(:),g.curvatures{{c}}(:)];end
dlmwrite('{out/'octave-geometry.csv'}',G,',','precision',17);
t=zeros(1,4);for k=1:4,tic;h=getOrientationStats(g,8);h=getLengthStats(h,8);h=getCurvatureStats(h,8);t(k)=toc;end
dlmwrite('{out/'octave-histogram-time.csv'}',t,',','precision',17);
dlmwrite('{out/'octave-histogram.csv'}',[h.normOrientationHistograms,h.normLengthHistograms,h.normCurvatureHistograms],',','precision',17);
"""
subprocess.run(['octave-cli','--quiet','--no-history','--eval',expr],check=True,timeout=60,capture_output=True,text=True)
sys.path.insert(0,str(root/'vendor/pyMLV'))
from MLVcode.computeContourProperties import computeContourProperties
from MLVcode.getOrientationStats import getOrientationStats
from MLVcode.getLengthStats import getLengthStats
from MLVcode.getCurvatureStats import getCurvatureStats
z=np.loadtxt(out/'input.csv',delimiter=',');n=int(z[:,0].max())
cells=np.empty((1,n),dtype=object)
for i in range(n):cells[0,i]=z[z[:,0]==i+1,1:]
v={'contours':cells,'numContours':np.array([[n]]),'imsize':np.loadtxt(out/'dimensions.csv',delimiter=',').reshape(1,2)}
def timed(f):
 times=[]
 for _ in range(4):
  start=time.perf_counter();value=f();times.append(time.perf_counter()-start)
 return value,times
g,tg=timed(lambda:computeContourProperties(v,['length','orientation','curvature']))
G=np.vstack([np.column_stack((np.ravel(g['lengths'][0][i]),np.ravel(g['orientations'][i]),np.ravel(g['curvatures'][i]))) for i in range(n)])
np.savetxt(out/'python-geometry.csv',G,delimiter=',',fmt='%.17g')
# Explicit schema normalization; raw getLengthStats fails on the returned (n,1) field.
g['contourLengths']=g['contourLengths'].reshape(-1)
def histogram():
 h,_,_,_=getOrientationStats(g,8);h,_,_,_=getLengthStats(h,8);h,_,_,_=getCurvatureStats(h,8)
 return np.column_stack((h['normOrientstionHistograms'],h['normLengthHistogram'],h['normCurvatureHistograms']))
H,th=timed(histogram);np.savetxt(out/'python-histogram.csv',H,delimiter=',',fmt='%.17g')
rows=[]
for backend in ['octave','python']:
 for stage,t in [('geometry',tg),('histogram',th)]:
  if backend=='octave':t=np.loadtxt(out/f'octave-{stage}-time.csv',delimiter=',').tolist()
  rows.append(dict(backend=backend,stage=stage,first_seconds=t[0],median_seconds=statistics.median(t[1:]),min_seconds=min(t[1:]),max_seconds=max(t[1:]),schema_adapter=backend=='python' and stage=='histogram'))
with (out/'timings-reference.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
receipt={'source':entry,'cleaning':np.loadtxt(out/'cleaning.csv',delimiter=',').tolist(),
 'cleaning_fields':['original_contours','retained_contours','retained_segments','zero_segments_removed'],
 'python_histogram_adapter':'reshape contourLengths from (n,1) to (n,); numerical functions unchanged',
 'reference_commits':{r['name']:r['commit'] for r in manifest['repositories']}}
(out/'provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'cleaning':receipt['cleaning'],'timings':rows},indent=2))
