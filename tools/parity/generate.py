#!/usr/bin/env python3
"""Generate bounded synthetic differential fixtures from untouched pinned sources."""
from pathlib import Path
import csv,json,hashlib,subprocess,tempfile,sys,platform,copy,traceback
import numpy as np
root=Path(__file__).resolve().parents[2]
out=root/'inst/fixtures/parity';out.mkdir(parents=True,exist_ok=True)
source=root/'vendor/MLV_toolbox/MLVcode'
sys.path.insert(0,str(root/'vendor/pyMLV'))
from MLVcode.computeContourProperties import computeContourProperties
from MLVcode.getOrientationStats import getOrientationStats
from MLVcode.getLengthStats import getLengthStats
from MLVcode.getCurvatureStats import getCurvatureStats
from MLVcode.rotateLinedrawing import rotateLinedrawing
from MLVcode.bwdist import bwdist
from MLVcode.computeAOF import computeAOF
cases=[];commands=[];python_jobs=[];failures=[]
def save(name,value):
 np.savetxt(out/name,np.atleast_2d(value),delimiter=',',fmt='%.17g')
def ocsv(name,expression):
 return f"dlmwrite('{out/name}',{expression},',','precision',17);"
def geometry(name,contours,domain):
 rows=np.vstack([np.column_stack((np.full(len(z),i+1),z)) for i,z in enumerate(contours)])
 save(name+'.input.csv',rows)
 cases.append(dict(family='geometry',name=name,domain=domain))
 cmd=f"Z=dlmread('{out/(name+'.input.csv')}',',');v=struct('imsize',[256,192],'numContours',max(Z(:,1)),'originalImage','synthetic','lineMethod','fixture');v.contours=cell(1,v.numContours);for c=1:v.numContours,v.contours{{c}}=Z(Z(:,1)==c,2:5);end;v=computeContourProperties(v,{{'length','orientation','curvature'}});"
 cmd+="G=[];for c=1:v.numContours,G=[G;v.lengths{c}(:),v.orientations{c}(:),v.curvatures{c}(:)];end;"
 cmd+=ocsv(name+'.octave.geometry.csv','G')
 for bins in [5,8,13]:
  cmd+=f"[v,~,~,~]=getOrientationStats(v,{bins});[v,~,~,~]=getLengthStats(v,{bins});[v,~,~,~]=getCurvatureStats(v,{bins});"
  for prop in ['orientation','length','curvature']:
   cmd+=ocsv(f'{name}.octave.{prop}{bins}.csv',f'v.{prop}Histograms')
   field='norm'+prop[0].upper()+prop[1:]+'Histograms'
   cmd+=ocsv(f'{name}.octave.norm_{prop}{bins}.csv','v.'+field)
 cmd+="rot=rotateLinedrawing(v,-37);G=vertcat(rot.contours{:});"+ocsv(name+'.octave.rotate37.csv','G')
 commands.append(cmd)
 def run_python():
  cells=np.empty((1,len(contours)),dtype=object)
  for i,z in enumerate(contours):cells[0,i]=z.copy()
  v={'contours':cells,'numContours':np.array([[len(contours)]]),'imsize':np.array([[256,192]]),
     'originalImage':np.array(['synthetic']),'lineMethod':np.array(['fixture'])}
  v=computeContourProperties(v,['length','orientation','curvature'])
  G=np.vstack([np.column_stack((np.ravel(v['lengths'][0][i]),np.ravel(v['orientations'][i]),np.ravel(v['curvatures'][i]))) for i in range(len(contours))])
  save(name+'.python.geometry.csv',G)
  for bins in [5,8,13]:
   for prop,fn,field,norm in [('orientation',getOrientationStats,'orientationHistograms','normOrientstionHistograms'),('length',getLengthStats,'lengthHistogram','normLengthHistogram'),('curvature',getCurvatureStats,'curvatureHistograms','normCurvatureHistograms')]:
    try:
     result,_,_,_=fn(copy.deepcopy(v),bins)
     save(f'{name}.python.{prop}{bins}.csv',result[field]);save(f'{name}.python.norm_{prop}{bins}.csv',result[norm])
    except Exception as error:
     failures.append({'case':name,'stage':prop+str(bins),'exception':type(error).__name__,'message':str(error),'location':traceback.extract_tb(error.__traceback__)[-1].line})
     if prop=='length':
      canonical=copy.deepcopy(v);canonical['contourLengths']=canonical['contourLengths'].reshape(-1)
      result,_,_,_=fn(canonical,bins)
      save(f'{name}.python.adapted_length{bins}.csv',result[field]);save(f'{name}.python.adapted_norm_length{bins}.csv',result[norm])
  rot=rotateLinedrawing(v,-37);save(name+'.python.rotate37.csv',np.vstack(rot['contours'][0]))
 python_jobs.append((name,run_python))
rng=np.random.default_rng(9201)
for seed in range(3):
 contours=[]
 for i in range(8):
  steps=rng.integers(-8,9,size=(13+i,2)).astype(float)
  steps[np.all(steps==0,axis=1),0]=1
  xy=np.vstack((np.array([[40+i*8,30+i*6]]),np.array([40+i*8,30+i*6])+np.cumsum(steps,axis=0)))
  contours.append(np.column_stack((xy[:-1],xy[1:])))
 geometry('integer'+str(seed),contours,'common_integer')
for name,scale in [('subpixel',.13),('fractional',.71)]:
 contours=[]
 for i in range(5):
  t=np.linspace(0,2.5+i*.2,12+i);xy=np.column_stack((30+i*20+scale*12*np.cos(t),40+scale*12*np.sin(t)))
  contours.append(np.column_stack((xy[:-1],xy[1:])))
 geometry(name,contours,'python_truncation')
# Diverse masks include nearest-neighbor ties, closed loops, endpoints and singletons.
for name in ['rectangle','ring','cross','diagonal','sparse','isolated']:
 h,w=25,31;rr,cc=np.indices((h,w));mask=np.zeros((h,w),bool)
 if name=='rectangle':mask[5:20,7]=1;mask[5:20,24]=1;mask[5,7:25]=1;mask[19,7:25]=1
 if name=='ring':mask=((rr-12)**2+(cc-15)**2>=49)&((rr-12)**2+(cc-15)**2<=81)
 if name=='cross':mask[4:21,15]=1;mask[12,4:27]=1
 if name=='diagonal':mask[np.arange(4,20),np.arange(5,21)]=1;mask[np.arange(4,20),np.arange(5,21)+5]=1
 if name=='sparse':mask=rng.random((h,w))<.035
 if name=='isolated':mask[7,8]=1;mask[16,22]=1
 save(name+'.input.csv',mask.astype(int));cases.append(dict(family='mat',name=name,domain='binary'))
 cmd=f"ink=dlmread('{out/(name+'.input.csv')}',',');M=computeMAT(1-ink,28);[D,I]=computeGradientVectorField(double(~bwmorph(ink,'skel',Inf)));"
 for key,expr in [('distance','M.distance_map'),('index','I'),('aof','M.AOF'),('skeleton','M.skeleton')]:
  cmd+=ocsv(name+'.octave.'+key+'.csv',expr)
 commands.append(cmd)
# EDT/AOF comparisons feed the same distance and nearest-label fields to each AOF kernel.
for name in ['features_random','features_ties']:
 mask=(rng.random((13,17))<.12) if name.endswith('random') else np.zeros((13,17),bool)
 if name.endswith('ties'):mask[3,4]=1;mask[3,12]=1;mask[9,4]=1;mask[9,12]=1
 save(name+'.input.csv',mask.astype(int));cases.append(dict(family='edt_aof',name=name,domain='nonempty'))
 cmd=f"F=dlmread('{out/(name+'.input.csv')}',',');[D,I]=bwdist(F);"+ocsv(name+'.octave.distance.csv','D')+ocsv(name+'.octave.index.csv','I')
 for eps in [1,2]:
  cmd+=f"alpha=(2*pi)/60;s=(0:59)';circle=[cos(alpha*s),sin(alpha*s)];A=computeAOF(D,I,circle,{eps});"+ocsv(f'{name}.octave.aof{eps}.csv','A')
 commands.append(cmd)
 def run_python(name=name,mask=mask):
  D,I=bwdist(mask);save(name+'.python.distance.csv',D);save(name+'.python.index.csv',I)
  D=np.loadtxt(out/(name+'.octave.distance.csv'),delimiter=',')
  labels=np.loadtxt(out/(name+'.octave.index.csv'),delimiter=',').astype(int)-1
  I=(labels%mask.shape[0])*mask.shape[1]+labels//mask.shape[0]
  alpha=(2*np.pi)/60;t=alpha*np.arange(60);circle=np.column_stack((np.cos(t),np.sin(t)))
  for eps in [1,2]:save(f'{name}.python.aof{eps}.csv',computeAOF(D,I,circle,eps))
 python_jobs.append((name,run_python))
for name in ['constant','varying','curved','terminal_bend','short']:
 n=21 if name!='short' else 2;x=np.arange(1,n+1,dtype=float);y=np.zeros(n);radius=np.full(n,3.)
 if name=='varying':radius=2+.07*x+.2*np.sin(x/3)
 if name=='curved':y=4*np.sin(x/4);radius=2+.1*x
 if name=='terminal_bend':y[-1]=100
 save(name+'.branch.csv',np.column_stack((x,y,radius)))
 cases.append(dict(family='branch',name=name,domain='fixed_movmean'))
 cmd=f"B=dlmread('{out/(name+'.branch.csv')}',',');b=struct('X',B(:,1),'Y',B(:,2),'Radius',B(:,3));"
 cmd+="G=[computeMATpropertyPerBranch(b,'parallelism'),computeMATpropertyPerBranch(b,'separation'),computeMATpropertyPerBranch(b,'mirror')];"
 commands.append(cmd+ocsv(name+'.octave.scores.csv','G'))
with tempfile.TemporaryDirectory(prefix='mlv-parity-') as td:
 td=Path(td)
 original=(source/'computeGradientVectorField.m').read_text()
 assert original.count('jj++m_Neighbors8')==2
 (td/'computeGradientVectorField.m').write_text(original.replace('jj++m_Neighbors8','jj + m_Neighbors8'))
 adapter="""function y=smoothdata(x,method,k)
if nargin~=3 || ~strcmp(method,'movmean'),error('Only explicit fixed-window movmean is adapted');end
y=zeros(size(x));n=numel(x);for i=1:n,lo=max(1,i-floor(k/2));hi=min(n,i+ceil(k/2)-1);y(i)=mean(x(lo:hi));end
end
"""
 (td/'smoothdata.m').write_text(adapter)
 program=f"pkg load image;addpath('{source}');addpath('{td}','-begin');\n"+'\n'.join(commands)
 (td/'run_fixtures.m').write_text(program)
 p=subprocess.run(['octave-cli','--quiet','--no-history',str(td/'run_fixtures.m')],capture_output=True,text=True,timeout=120)
 print(p.stdout);print(p.stderr,file=sys.stderr)
 if p.returncode:raise SystemExit(p.returncode)
 for name,job in python_jobs:
  try:job()
  except Exception as error:failures.append({'case':name,'exception':type(error).__name__,'message':str(error)})
with (out/'cases.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=['family','name','domain'],lineterminator='\n');writer.writeheader();writer.writerows(cases)
sources=json.loads((root/'vendor/SOURCES.json').read_text())
metadata={'seed':9201,'numpy':np.__version__,'python':platform.python_version(),
 'source_commits':{s['name']:s['commit'] for s in sources['repositories']},
 'adapters':{'gradient':'two jj++m_Neighbors8 -> jj + m_Neighbors8 substitutions',
 'smoothdata':adapter,'smoothdata_sha256':hashlib.sha256(adapter.encode()).hexdigest(),
 'python_length_histogram':'After preserving raw failure, reshape contourLengths from (n,1) to (n,) for the adapted length-histogram comparison only'},
 'python_failures':failures,
 'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*.csv')}}
(out/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps({'cases':len(cases),'files':len(metadata['files']),'python_failures':failures},indent=2))
