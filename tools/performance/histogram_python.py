"""Fresh adapted Python comparison for the isolated forest histogram experiment."""
from pathlib import Path
import sys,time,csv,json,statistics
import numpy as np
root=Path(__file__).resolve().parents[2]
work=Path(sys.argv[1]);dest=Path(sys.argv[2])
sys.path.insert(0,str(root/'vendor/pyMLV'))
from MLVcode.computeContourProperties import computeContourProperties
from MLVcode.getOrientationStats import getOrientationStats
from MLVcode.getLengthStats import getLengthStats
from MLVcode.getCurvatureStats import getCurvatureStats
z=np.loadtxt(work/'forest.csv',delimiter=',');n=int(z[:,0].max())
cells=np.empty((1,n),dtype=object)
for i in range(n):cells[0,i]=z[z[:,0]==i+1,1:]
v={'contours':cells,'numContours':np.array([[n]]),
   'imsize':np.loadtxt(work/'dimensions.csv',delimiter=',').reshape(1,2)}
g=computeContourProperties(v,['length','orientation','curvature'])
raw_failure=None
try:getLengthStats(g,8)
except ValueError as e:raw_failure=str(e)
assert raw_failure is not None
g['contourLengths']=g['contourLengths'].reshape(-1)
def histogram():
 h,_,_,_=getOrientationStats(g,8)
 h,_,_,_=getLengthStats(h,8)
 h,_,_,_=getCurvatureStats(h,8)
 return np.column_stack((h['normOrientstionHistograms'],h['normLengthHistogram'],h['normCurvatureHistograms']))
actual=histogram();ref=np.loadtxt(work/'R-histogram.csv',delimiter=',')
error=float(np.max(np.abs(actual-ref)))
assert actual.shape==ref.shape and np.all(np.isfinite(actual)) and error<1e-8
rows=[]
for round in range(1,4):
 number=1
 while True:
  start=time.perf_counter()
  for _ in range(number):histogram()
  elapsed=time.perf_counter()-start
  if elapsed>=.12 or number>=512:break
  number*=2
 times=[]
 for _ in range(5):
  start=time.perf_counter()
  for _ in range(number):histogram()
  times.append((time.perf_counter()-start)/number)
 rows.append(dict(round=round,workload='forest',method='adapted_python',
   median_seconds=statistics.median(times),min_seconds=min(times),
   max_seconds=max(times),iterations=number))
with (dest/'python-timings.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),lineterminator='\n');w.writeheader();w.writerows(rows)
(dest/'python-validation.json').write_text(json.dumps(dict(source='14b3183b91146b1021eb7ef97db75ad703b06832',
 python=sys.version,numpy=np.__version__,max_abs_error=error,tolerance=1e-8,
 raw_length_failure=raw_failure,adapter='reshape contourLengths from (n,1) to (n,) before histogram timing'),indent=2)+'\n')
print(json.dumps(rows,indent=2))
