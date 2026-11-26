"""Recover the curve indices of each photo-only M3 drawing (development only).

Applies the same rule as screen.py: complete curves by decreasing key until
the next would exceed the photo's FDoG ink. Checks the number of curves taken
against source-pool-photo-scores.csv.

Usage: python m3_selection.py SOURCE_POOL_DIR PHOTO_SCORES_CSV OUT_CSV
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

pool_dir, scores_csv, out_csv = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
scores = {r["id"]: r for r in csv.DictReader(open(scores_csv))}
rows = []
for pid, r in sorted(scores.items()):
    pool = json.loads((pool_dir / "pools" / f"{pid}.json").read_text())
    fd = [np.asarray(p, dtype=np.int64) for p in pool["fdog"]]
    cap = int(np.unique(np.concatenate(fd)).size) if fd else 0
    assert cap == int(r["fdog_ink_budget"]), (pid, cap, r["fdog_ink_budget"])
    for m, key in (("fdog512", "fdog512_length"), ("sed", "sed_strength")):
        paths = [np.asarray(p, dtype=np.int64) for p in pool[m]]
        n = max([cap] + [int(p.max()) + 1 for p in paths if len(p)])
        seen = np.zeros(n, bool)
        used = 0
        taken = []
        for j in np.argsort(-np.asarray(pool[key], float), kind="stable"):
            new = paths[j][~seen[paths[j]]]
            if used + len(new) > cap:
                break
            seen[new] = True
            used += len(new)
            taken.append(int(j) + 1)
        assert len(taken) == int(r[f"m3_{m}_curves"]), (
            pid,
            m,
            len(taken),
            r[f"m3_{m}_curves"],
        )
        assert used == int(r[f"m3_{m}_ink"]), (pid, m, used, r[f"m3_{m}_ink"])
        rows.append(dict(id=pid, source=m, curves=" ".join(map(str, taken))))
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id", "source", "curves"])
    w.writeheader()
    w.writerows(rows)
print(len(rows), "selections written; all curve and ink counts match the screen")
