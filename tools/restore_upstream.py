#!/usr/bin/env python3
"""Restore an untouched full upstream checkout without overwriting existing data."""
import argparse,json,subprocess
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('repository',choices=['MLV_toolbox','pyMLV'])
parser.add_argument('--destination',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'vendor/SOURCES.json').read_text())
spec=next(r for r in manifest['repositories'] if r['name']==args.repository)
if args.destination.exists(): parser.error('destination already exists')
subprocess.run(['git','clone','--no-checkout',spec['origin'],str(args.destination)],check=True)
subprocess.run(['git','-C',str(args.destination),'checkout','--detach',spec['commit']],check=True)
tree=subprocess.check_output(['git','-C',str(args.destination),'rev-parse','HEAD^{tree}'],text=True).strip()
if tree!=spec['tree']: raise SystemExit('Restored Git tree does not match the recorded tree')
print('Restored verified upstream tree:',tree)
