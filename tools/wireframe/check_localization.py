"""Compare observed extra matches to shifts of the same additional detections."""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial.distance import cdist

from score import positions, read_csv, write_csv


def cardinality(refs, predictions):
    if not len(refs) or not len(predictions):
        return 0
    graph = csr_matrix(cdist(refs, predictions) <= 4)
    return int((maximum_bipartite_matching(graph, perm_type="column") >= 0).sum())


def check(folder, panel, repeats=199):
    sites = read_csv(folder / f"support-{panel}-sites.csv")
    events = read_csv(folder / f"support-{panel}-predictions.csv")
    scored = read_csv(folder / f"support-{panel}-scores.csv")
    expected = {(r["id"], r["type"], r["method"]): int(r["matched"]) for r in scored}
    methods = ("mlsd_mask", "supported")
    from_sites = {k: 0 for k in ("T", "X")}
    shifted = {k: np.zeros(repeats, dtype=int) for k in from_sites}
    rng = np.random.default_rng(731099)
    for photo in read_csv(folder / "panel.csv"):
        i = photo["id"]
        if (i, "X", "supported") not in expected:
            continue
        size = np.array([int(photo["width"]), int(photo["height"])])
        # A common translation preserves clustering and relative T/X positions
        # within a photo. Wrap at each image edge to retain all additions.
        offsets = rng.uniform(0.2, 0.8, size=(repeats, 2)) * size
        for k in ("T", "X"):
            pts = positions([r for r in sites if r["id"] == i and r["kind"] == k])
            pred = [
                positions(
                    [
                        r
                        for r in events
                        if (r["id"], r["type"], r["method"]) == (i, k, m)
                    ]
                )
                for m in methods
            ]
            base, proposed = pred
            remaining = Counter(map(tuple, proposed))
            remaining.subtract(map(tuple, base))
            if any(n < 0 for n in remaining.values()):
                raise ValueError("A baseline detection was lost")
            extras = np.array(list(remaining.elements())).reshape(-1, 2)
            counts = [cardinality(pts, z) for z in pred]
            assert counts == [expected[i, k, m] for m in methods]
            from_sites[k] += counts[1] - counts[0]
            if len(extras) and len(pts):
                for n, offset in enumerate(offsets):
                    moved = ((extras - 1 + offset) % size) + 1
                    shifted[k][n] += (
                        cardinality(pts, np.concatenate((base, moved))) - counts[0]
                    )
    return [
        dict(
            panel=panel,
            kind=k,
            extra_matches=from_sites[k],
            shifts=repeats,
            null_median=float(np.median(v)),
            null_p95=int(np.quantile(v, 0.95, method="higher")),
        )
        for k, v in shifted.items()
    ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "panel", choices=("development", "issue", "replication", "confirmation")
    )
    args = parser.parse_args()
    rows = check(args.directory, args.panel)
    write_csv(
        args.directory / f"support-{args.panel}-shift-check.csv", rows, list(rows[0])
    )
