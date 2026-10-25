"""Verified consolidation of the two user-supplied reference clones."""
from pathlib import Path
import hashlib, json, os, shutil, subprocess, tempfile, datetime
ROOT=Path.cwd().resolve()
SPECS=[
 {'name':'MLV_toolbox','old':'MLV_toolbox','commit':'f12418a5e03c2f38e222fed1dffccf42b227cd29','tree':'74b4bc8f3be1b34823be9c3dbfde2491f87b6fc2','origin':'https://github.com/bwlabToronto/MLV_toolbox'},
 {'name':'pyMLV','old':'vendor/pyMLV','commit':'14b3183b91146b1021eb7ef97db75ad703b06832','tree':'f8f84ad328172e6ef3d2730339de4011fe57b9b0','origin':'https://github.com/bwlabToronto/pyMLV'}]
def git(base,*args):
 p=subprocess.run(['git','--no-optional-locks','-C',str(base),*args],capture_output=True,text=True,check=True)
 return p.stdout.strip()
def digest(p):
 with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def size(p):
 return sum(f.stat().st_size for f in p.rglob('*') if f.is_file())
data=ROOT/'reference-data/mlv'
if data.exists(): raise SystemExit('Existing canonical data collection; aborting.')
for spec in SPECS:
 base=ROOT/spec['old']
 assert git(base,'rev-parse','HEAD')==spec['commit']
 assert git(base,'rev-parse','HEAD^{tree}')==spec['tree']
 assert git(base,'config','--get','remote.origin.url')==spec['origin']
 assert git(base,'status','--porcelain=v1','--untracked-files=normal')==''
 spec['refs']=git(base,'for-each-ref','--format=%(refname) %(objectname)').splitlines()
 for ref in spec['refs']:
  if ref.startswith('refs/heads/'):
   assert ref.split()[1] in {r.split()[1] for r in spec['refs'] if r.startswith('refs/remotes/origin/')}
 spec['before_bytes']=size(base)
# Revalidate the reviewed duplicate files before consolidation.
review=json.loads((ROOT/'VENDOR_DUPLICATES.json').read_text())
for group in review['exact_duplicates']['groups']:
 for rel in group['paths']:
  p=ROOT/rel
  assert p.is_file() and p.stat().st_size==group['size_bytes'] and digest(p)==group['sha256'], rel
stage=Path(tempfile.mkdtemp(prefix='.mlv-cleanup-',dir=ROOT/'vendor'))
retained=[]
for spec in SPECS:
 old=ROOT/spec['old']; dest=stage/spec['name']; dest.mkdir()
 for p in old.rglob('*'):
  if not p.is_file() or '.git' in p.relative_to(old).parts: continue
  rel=p.relative_to(old); s=rel.as_posix()
  keep=False
  if spec['name']=='MLV_toolbox':
   keep=(s in ['README.md','setup.m'] or
         (s.startswith('MLVcode/') and p.suffix=='.m' and not s.startswith('MLVcode/edges-master/')) or
         s.startswith('3rdParty/contour_tracing/') or
         s in ['demos/Frontiers2023_Figures.m','demos/testPatterns.svg'])
  else:
   keep=(s in ['README.md','LICENSE','requirements.txt','mainpage.py'] or
         (s.startswith('MLVcode/') and p.suffix=='.py') or
         (s.startswith('demos/') and p.suffix=='.ipynb'))
  if not keep: continue
  q=dest/rel; q.parent.mkdir(parents=True,exist_ok=True)
  original=digest(p); transform=None
  if p.suffix=='.ipynb':
   nb=json.loads(p.read_text())
   for cell in nb.get('cells',[]):
    if cell.get('cell_type')=='code':
     cell['outputs']=[]; cell['execution_count']=None
   q.write_text(json.dumps(nb,indent=1,ensure_ascii=False)+'\n')
   transform='Cleared saved notebook outputs and execution counts; source cells retained.'
  else: shutil.copy2(p,q)
  retained.append({'repository':spec['name'],'path':s,'original_sha256':original,'snapshot_sha256':digest(q),'bytes':q.stat().st_size,'transformation':transform})
# Assemble a single retained copy of data; use renames, not a second bulk copy.
data.mkdir(parents=True)
orig_mat=ROOT/'MLV_toolbox'; orig_py=ROOT/'vendor/pyMLV'
data_records=[]
for p in (orig_mat/'dataSets').rglob('*'):
 if p.is_file():
  rel=p.relative_to(orig_mat)
  data_records.append({'source_repository':'MLV_toolbox','source_commit':SPECS[0]['commit'],'source_path':rel.as_posix(),'destination':str(data/'dataSets'/p.relative_to(orig_mat/'dataSets')),'bytes':p.stat().st_size,'sha256':digest(p)})
shutil.move(str(orig_mat/'dataSets'),str(data/'dataSets'))
examples=data/'examples'; examples.mkdir()
for spec,rel,name in [(SPECS[0],'demos/Forest_DBW.mat','Forest_DBW.mat'),(SPECS[1],'MLVcode/MAT.mat','MAT.mat'),(SPECS[1],'MLVcode/imgLD.mat','imgLD.mat')]:
 p=ROOT/spec['old']/rel; q=examples/name
 record={'source_repository':spec['name'],'source_commit':spec['commit'],'source_path':rel,'destination':str(q),'bytes':p.stat().st_size,'sha256':digest(p)}
 shutil.move(str(p),str(q)); data_records.append(record)
for r in data_records:
 p=Path(r['destination']); assert p.stat().st_size==r['bytes'] and digest(p)==r['sha256']
manifest={'schema_version':1,'created':datetime.datetime.now(datetime.timezone.utc).isoformat(),'remote_custody_verified':'GitHub commit API verified both source commits and root Git trees in this session before cleanup.','repositories':SPECS,'retained_source_files':retained,'shared_data':{'path':str(data),'filesystem_device':os.stat(data).st_dev,'volume':'APFS Data volume, disk3s5 at assessment','files':data_records},'notes':['Source snapshots exclude Git object stores, datasets, generated documentation/figures, model weights and saved notebook outputs.','Full originals are recoverable from pinned upstream commits using tools/restore_upstream.py.','Retained source files preserve their original copyright and license notices; root license claims do not override per-file notices.']}
(ROOT/'vendor/SOURCES.json').write_text(json.dumps(manifest,indent=2)+'\n')
(data/'MANIFEST.json').write_text(json.dumps(manifest['shared_data'],indent=2)+'\n')
# All original changes since preflight above are the explicitly recorded data moves.
quarantine=ROOT/'.cleanup-staging'; quarantine.mkdir(exist_ok=False)
for spec in SPECS:
 old=ROOT/spec['old']; parked=quarantine/spec['name']
 old.rename(parked)
 (stage/spec['name']).rename(ROOT/'vendor'/spec['name'])
stage.rmdir()
for r in retained:
 assert digest(ROOT/'vendor'/r['repository']/r['path'])==r['snapshot_sha256']
# Delete only the two exact, owned, verified replacement trees.
for spec in SPECS:
 parked=quarantine/spec['name']
 assert parked.parent==ROOT/'.cleanup-staging' and (parked/'.git').is_dir()
 assert git(parked,'rev-parse','HEAD')==spec['commit']
 shutil.rmtree(parked)
quarantine.rmdir()
after_source=size(ROOT/'vendor'); after_data=size(data)
receipt={'date':datetime.datetime.now(datetime.timezone.utc).isoformat(),'before_logical_bytes':sum(s['before_bytes'] for s in SPECS),'after_vendor_bytes':after_source,'after_shared_data_bytes':after_data,'logical_bytes_removed':sum(s['before_bytes'] for s in SPECS)-after_source-after_data,'source_files_verified':len(retained),'retained_data_files_verified':len(data_records),'source_manifest':'vendor/SOURCES.json','shared_data_manifest':str(data/'MANIFEST.json'),'restoration':'tools/restore_upstream.py','external_custody':[(s['origin'],s['commit'],s['tree']) for s in SPECS],'full_references_removed':True}
(ROOT/'CLEANUP_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
