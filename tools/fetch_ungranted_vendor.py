#!/usr/bin/env python3
"""Fetch retained MLV_toolbox files that lack an explicit license grant.

The public repository omits upstream files whose license is unresolved
(terms NOASSERTION in docs/evidence/attribution/vendor-license-inventory.csv).
Parity, benchmark and license-audit tools need the complete retained snapshot.
This script clones the pinned upstream commit into a temporary directory,
copies only the missing NOASSERTION files into vendor/MLV_toolbox, and verifies
every copied file against its SHA-256 in vendor/SOURCES.json. Existing files
are never overwritten. The fetched files are ignored by Git.
"""
import csv, hashlib, json, shutil, subprocess, sys, tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / "vendor/SOURCES.json").read_text())
spec = next(r for r in manifest["repositories"] if r["name"] == "MLV_toolbox")
hashes = {
    e["path"]: e["snapshot_sha256"]
    for e in manifest["retained_source_files"]
    if e["repository"] == "MLV_toolbox"
}
inventory = root / "docs/evidence/attribution/vendor-license-inventory.csv"
wanted = [
    r["path"]
    for r in csv.DictReader(inventory.open())
    if r["repository"] == "MLV_toolbox" and r["terms"] == "NOASSERTION"
]
dest = root / "vendor/MLV_toolbox"
missing = [p for p in wanted if not (dest / p).exists()]
if not missing:
    print("All", len(wanted), "ungranted files already present.")
    sys.exit(0)

with tempfile.TemporaryDirectory() as tmp:
    clone = Path(tmp) / "MLV_toolbox"
    subprocess.run(
        ["git", "clone", "--quiet", "--no-checkout", spec["origin"], str(clone)],
        check=True,
    )
    subprocess.run(
        ["git", "-C", str(clone), "checkout", "--quiet", "--detach", spec["commit"]],
        check=True,
    )
    tree = subprocess.check_output(
        ["git", "-C", str(clone), "rev-parse", "HEAD^{tree}"], text=True
    ).strip()
    if tree != spec["tree"]:
        raise SystemExit("Upstream tree does not match the recorded tree")
    for p in missing:
        data = (clone / p).read_bytes()
        if hashlib.sha256(data).hexdigest() != hashes[p]:
            raise SystemExit("Hash mismatch for " + p)
        target = dest / p
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
print("Fetched and verified", len(missing), "of", len(wanted), "ungranted files.")
