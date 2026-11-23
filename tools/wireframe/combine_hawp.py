"""Evaluate corroborated additions to the M-LSD junction table on fresh photographs."""
import argparse
from pathlib import Path

from scipy.spatial.distance import cdist

from evaluate_hawp import compare, selection_key
from score import (
    ROOT,
    positions,
    original,
    read_csv,
    score,
    supplement,
    totals,
    write_csv,
)

METHODS = ("mlsd_mask", "corroborated", "union")


def combine(primary, enhanced, vectors, traced):
    # A secondary detector may add a new site, but cannot assign a second
    # label to a primary site within the benchmark's 4-pixel resolution.
    primary = [r for r in primary if r["type"] in ("T", "X")]
    additions = [r for r in traced if r["type"] in ("T", "X")]
    a = positions(additions)
    close_to_primary = (cdist(a, positions(primary)) <= 4).any(axis=1)
    witnesses = enhanced + vectors
    nearby = cdist(a, positions(witnesses)) <= 4
    supplemented = list(primary)
    union = list(primary)
    for i, candidate in enumerate(additions):
        if close_to_primary[i]:
            continue
        union.append(candidate)
        if any(
            near and w["type"] == candidate["type"]
            for near, w in zip(nearby[i], witnesses)
        ):
            supplemented.append(candidate)
    return [
        dict(r, method=name)
        for name, rows in zip(METHODS, (primary, supplemented, union))
        for r in rows
    ]


def benchmark(folder, split, photos, reviews, cutoff, evaluate_methods):
    baseline = read_csv(folder / "junctions.csv")
    hawp = read_csv(folder / f"hawp-{split}-junctions.csv")
    events = []
    for i, _, _ in photos:

        def method(rows, name):
            return [r for r in rows if r["id"] == i and r["method"] == name]

        events.extend(
            combine(
                method(baseline, "mlsd_mask"),
                method(baseline, "mlsd_clahe"),
                method(hawp, f"hawp_vector_{cutoff:.2f}"),
                method(hawp, f"hawp_mask_{cutoff:.2f}"),
            )
        )
    scores, negatives = score(
        folder,
        photos,
        reviews,
        prefix=f"combined-{split}-",
        methods=evaluate_methods,
        detections=events,
    )
    return [totals(split, name, scores, negatives) for name in evaluate_methods]


def main(folder, manifest, replication):
    chosen = [
        r
        for r in read_csv(folder / "hawp-summary.csv")
        if r["panel"] == "development" and r["selected"] == "1"
    ]
    (result,) = chosen
    cutoff = float(result["cutoff"])
    issue, negatives = original(folder, manifest)
    issue_ids = {i for i, _, _ in issue}
    dev = [r for r in supplement(folder, groups=("sample",)) if r[0] not in issue_ids]
    dev_ids = {r[0] for r in dev}
    reviewed = [
        r
        for r in read_csv(ROOT / "docs/evidence/wireframe/near-miss-review.csv")
        if r["id"] in dev_ids
    ]
    development = benchmark(folder, "development", dev, reviewed, cutoff, METHODS)
    baseline = development[0]

    def count(row):
        return row["T_detected"] + row["X_detected"]

    eligible = [
        r
        for r in development
        if r["flagged"] <= baseline["flagged"] and count(r) <= 1.20 * count(baseline)
    ]
    selected = max(eligible, key=selection_key)["method"]
    print("Development selection:", selected, development, flush=True)
    test_methods = tuple(dict.fromkeys(("mlsd_mask", selected)))
    frozen = benchmark(folder, "issue", issue, negatives, cutoff, test_methods)

    fresh = supplement(replication, groups=("replication",))
    fresh_ids = {r[0] for r in fresh}
    if len(fresh_ids) != 24 or fresh_ids & (dev_ids | issue_ids):
        raise ValueError("Replication roster must be 24 distinct unseen photos")
    fresh_negatives = read_csv(
        ROOT / "docs/evidence/wireframe/replication-gap-review.csv"
    )
    compare(replication, "replication", fresh, fresh_negatives, [cutoff])
    heldout = benchmark(
        replication, "replication", fresh, fresh_negatives, cutoff, test_methods
    )
    all_rows = development + frozen + heldout
    for row in all_rows:
        row["selected"] = int(row["method"] == selected)
    write_csv(folder / "combined-summary.csv", all_rows, list(all_rows[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("replication", type=Path)
    args = parser.parse_args()
    main(args.directory.resolve(), args.manifest, args.replication.resolve())
