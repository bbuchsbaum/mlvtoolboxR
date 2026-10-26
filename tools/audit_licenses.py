#!/usr/bin/env python3
"""Inventory notices in the pinned reference snapshots; this is not legal clearance."""
import collections, csv, hashlib, json, re
from pathlib import Path
root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / "vendor/SOURCES.json").read_text())
out = root / "docs/evidence/attribution"
out.mkdir(parents=True, exist_ok=True)
rows = []
for entry in manifest["retained_source_files"]:
    path = root / "vendor" / entry["repository"] / entry["path"]
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != entry["snapshot_sha256"]:
        raise SystemExit("Vendored snapshot changed: " + str(path.relative_to(root)))
    text = data.decode("utf-8", errors="replace")
    if "GNU General Public License" in text and re.search(r"(?:version|Version) 3", text):
        terms, basis = "GPL-3.0-or-later", "explicit per-file GPL v3-or-later notice"
    elif "Permission is hereby granted" in text and "Peter Kovesi" in text:
        terms, basis = "LicenseRef-Kovesi-permissive", "explicit abbreviated permission notice; preserve verbatim"
    elif entry["repository"] == "pyMLV":
        terms, basis = "MIT-root", "repository LICENSE; retain any individual copyright notices"
    else:
        terms, basis = "NOASSERTION", "no explicit license grant found in retained file; no root license at pin"
    credits = [line.strip() for line in text.splitlines()
               if re.search(r"copyright|Peter Kovesi|Aravind Narayanan", line, re.I)]
    rows.append(dict(repository=entry["repository"], path=entry["path"],
                     sha256=digest, terms=terms, basis=basis,
                     notices=" | ".join(credits), package_distribution="excluded"))
with (out / "vendor-license-inventory.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
summary = dict(
    files=len(rows), all_snapshot_hashes_match=True,
    repositories={r["name"]:{"url":r["origin"],"commit":r["commit"]}
                  for r in manifest["repositories"]},
    classifications={repo:dict(collections.Counter(r["terms"] for r in rows if r["repository"]==repo))
                     for repo in ["MLV_toolbox","pyMLV"]},
    limitation="Notice inventory only. NOASSERTION requires permission/lineage review before redistribution or further code reuse; acknowledgments do not supply permission.")
(out / "audit.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
