#!/usr/bin/env python3
"""Verify retained reference bytes without modifying sources or data."""
from pathlib import Path
import hashlib,json,sys,datetime
root=Path(__file__).resolve().parents[1]
sources=json.loads((root/'vendor/SOURCES.json').read_text())
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
checks=[]
for entry in sources['retained_source_files']:
    checks.append((root/'vendor'/entry['repository']/entry['path'],
                   entry['snapshot_sha256'],entry['bytes'],'source'))
for entry in sources['shared_data']['files']:
    # Paths are recorded at cleanup. An optional argument relocates the data root.
    base=Path(sys.argv[1]) if len(sys.argv)>1 else Path(sources['shared_data']['path'])
    rel=Path(entry['destination']).relative_to(sources['shared_data']['path'])
    checks.append((base/rel,entry['sha256'],entry['bytes'],'data'))
failures=[]
for path,sha,size,kind in checks:
    if not path.is_file():failures.append({'path':str(path),'reason':'missing'})
    elif path.stat().st_size!=size or digest(path)!=sha:
        failures.append({'path':str(path),'reason':'size or SHA256 mismatch'})
result={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'sources':sum(z[3]=='source' for z in checks),
        'data_files':sum(z[3]=='data' for z in checks),
        'data_bytes':sum(z[2] for z in checks if z[3]=='data'),
        'failures':failures}
print(json.dumps(result,indent=2))
raise SystemExit(bool(failures))
