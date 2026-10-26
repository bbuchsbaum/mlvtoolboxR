#!/usr/bin/env python3
"""Benchmarks only runnable Python paths with correct synthetic outputs."""
from pathlib import Path
import sys,time,statistics,json,csv
import numpy as np,scipy
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'vendor/pyMLV'))
from MLVcode.computeContourProperties import computeContourProperties
from MLVcode.bwdist import bwdist
def timed(f):
 start=time.perf_counter();f();cold=time.perf_counter()-start
 number=1
 while True:
  start=time.perf_counter()
  for _ in range(number):f()
  elapsed=time.perf_counter()-start
  if elapsed>=.12 or number>=1024:break
  number*=2
 times=[]
 for _ in range(5):
  start=time.perf_counter()
  for _ in range(number):f()
  times.append((time.perf_counter()-start)/number)
 return dict(cold_seconds=cold,median_seconds=statistics.median(times),min_seconds=min(times),max_seconds=max(times),iterations=number)
rows=[]
for contours,segments in [(100,1),(1000,1),(1,10000)]:
 cells=np.empty((1,contours),dtype=object)
 for i in range(contours):
  x=np.arange(segments,dtype=float)+1
  cells[0,i]=np.column_stack([x,np.ones(segments)*(i+1),x+1,np.ones(segments)*(i+1)])
 drawing={'numContours':np.array([[contours]]),'contours':cells}
 def geometry():return computeContourProperties(drawing,['orientation','length','curvature'])
 result=geometry()
 assert all(np.all(z==0) for z in result['orientations'])
 assert all(np.all(np.asarray(z)==0) for z in result['curvatures'])
 assert np.all(np.concatenate(result['lengths'][0])==1)
 rows.append(dict(stage='contour_properties',size=contours,segments_per_contour=segments,**timed(geometry)))
for n in [32,64,128,512]:
 w=n+7;feature=np.zeros((n,w),dtype=bool);r=round(n/2)-1;c=round(w/2)-1;feature[r,c]=True
 D,I=bwdist(feature);rr,cc=np.indices(feature.shape)
 assert np.max(np.abs(D-np.hypot(rr-r,cc-c)))<1e-12
 assert np.all(I==r*w+c)
 rows.append(dict(stage='distance_transform',size=n,segments_per_contour=0,**timed(lambda:bwdist(feature))))
for row in rows:
 row.update(runtime=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,source='14b3183b91146b1021eb7ef97db75ad703b06832')
out=root/'docs/evidence/benchmark-python.csv'
with out.open('w') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print(json.dumps(rows,indent=2))
