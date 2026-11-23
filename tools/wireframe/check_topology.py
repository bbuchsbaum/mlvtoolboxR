"""Independently check same-type matches and bootstrap the paired photo recall change."""
import argparse
from pathlib import Path

import numpy as np

from check_localization import cardinality
from score import positions, read_csv, write_csv


def paired(folder, panel, combine=()):
    prefix = f"topology-{panel}-"
    sources = (
        [(folder / name, name) for name in combine] if combine else [(folder, panel)]
    )
    scores, events, sites = [], [], []
    ids = []
    for source, label in sources:
        name = f"topology-{label}-"
        table = read_csv(source / (name + "scores.csv"))
        photo_ids = list(dict.fromkeys(r["id"] for r in table))
        if set(ids) & set(photo_ids):
            raise ValueError("Duplicated photo in paired comparison")
        ids.extend(photo_ids)
        scores.extend(table)
        events.extend(read_csv(source / (name + "predictions.csv")))
        sites.extend(read_csv(source / (name + "sites.csv")))
    rng = np.random.default_rng(812711)
    sample = rng.integers(0, len(ids), size=(10000, len(ids)))
    rows = []
    for kind in ("T", "X"):
        counts = []
        for i in ids:
            refs = positions([r for r in sites if r["id"] == i and r["kind"] == kind])
            pair = []
            for method in ("supported", "topology"):
                pred = positions(
                    [
                        r
                        for r in events
                        if r["id"] == i and r["type"] == kind and r["method"] == method
                    ]
                )
                actual = cardinality(refs, pred)
                (logged,) = (
                    r
                    for r in scores
                    if (r["id"], r["type"], r["method"]) == (i, kind, method)
                )
                if [len(refs), len(pred), actual] != [
                    int(logged[k]) for k in ("reference", "detected", "matched")
                ]:
                    raise ValueError(f"Disagreement on {i} {method} {kind}")
                pair.append(actual)
            counts.append([len(refs), *pair])
        ref, old, new = np.array(counts).T
        diff = new - old
        denominators = ref[sample].sum(axis=1)
        valid = denominators > 0
        delta = (diff[sample].sum(axis=1)[valid] / denominators[valid]) * 100
        lo, hi = np.quantile(delta, [0.025, 0.975])
        recall = new[sample].sum(axis=1)[valid] / denominators[valid] * 100
        recall_lo, recall_hi = np.quantile(recall, [0.025, 0.975])
        rows.append(
            dict(
                panel=panel,
                kind=kind,
                photos=len(ids),
                reference=ref.sum(),
                previous=old.sum(),
                candidate=new.sum(),
                delta_pp=100 * diff.sum() / ref.sum(),
                lo_pp=lo,
                hi_pp=hi,
                recall_lo_pct=recall_lo,
                recall_hi_pct=recall_hi,
                photos_gained=(diff > 0).sum(),
                photos_lost=(diff < 0).sum(),
            )
        )
    write_csv(folder / (prefix + "paired.csv"), rows, list(rows[0]))
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("panel")
    parser.add_argument(
        "--combine",
        nargs="+",
        default=[],
        help="Combine disjoint named panels under the given folder",
    )
    args = parser.parse_args()
    print(paired(args.folder, args.panel, args.combine))
