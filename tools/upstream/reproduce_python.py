#!/usr/bin/env python3
"""Synthetic Python upstream checks; leaves the reference snapshot untouched."""
from pathlib import Path
import sys, json, ast, importlib.util
import numpy as np
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'vendor/pyMLV'))
sys.path.insert(0,str(root/'vendor/pyMLV/MLVcode'))
from MLVcode.computeOrientation import computeOrientation
from MLVcode.getHorizontalVerticalStats import getHorizontalVerticalStats
contours=np.empty((1,1),dtype=object)
contours[0,0]=np.array([[1.1,1.1,1.9,1.9]])
v={'numContours':np.array([[1]]),'contours':contours}
result=computeOrientation(v)['orientations'][0][0]
print('MLV-P001 subpixel orientation:',result,'expected:',315.)
contours[0,0]=np.array([[1.,1.,11.,1.]])
try:
    getHorizontalVerticalStats({'numContours':np.array([[1]]),'contours':contours})
except Exception as e:
    print('MLV-P002 horizontal/vertical histogram:',type(e).__name__,str(e))
stems={p.stem for p in (root/'vendor/pyMLV/MLVcode').glob('*.py')}
missing=[]
for name in ['getBranchDerivative','computeMATpropertyPerBranch','traceSkeleton']:
    tree=ast.parse((root/'vendor/pyMLV/MLVcode'/f'{name}.py').read_text())
    for node in ast.walk(tree):
        if isinstance(node,ast.ImportFrom) and node.module in {'smoothData','diff','fitLineSegments','GetConSeg'} and node.module not in stems:
            missing.append((name,node.module))
print('MLV-P003 missing local dependencies:',json.dumps(missing))
# Array semantics witness for the alias operation present in the upstream source.
original=np.array([[0.,30.]])
flux=original.copy();skeleton=flux
skeleton[skeleton<18]=0;skeleton[skeleton>=18]=1
print('MLV-P004 alias semantics witness: retained_flux=',flux.tolist(),'original=',original.tolist())
