"""Replay the York gap challenge, including a ranking comparison at equal T recall."""

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label
from scipy.spatial.distance import cdist

from contour_metrics import junction_counts
from score import ROOT, cross, positions, read_csv, write_csv

EVIDENCE = ROOT / "docs/evidence/wireframe/native-arms"


def drawings_bridge(mask, stem, end, target):
    """Require one connected ink component from the designated stem into the bar.

    Constrain connectivity to narrow corridors around the annotated stem and bar.
    End and target are one-based points at the nominal gap ends.
    """
    u = end - stem
    u = u / np.linalg.norm(u)
    start = end - 3 * u
    a, b = target
    v = b - a
    v /= np.linalg.norm(v)
    center = end + u * (cross(a - end, v) / cross(u, v))
    yy, xx = np.mgrid[1 : mask.shape[0] + 1, 1 : mask.shape[1] + 1]

    def axes(origin, direction):
        dx, dy = xx - origin[0], yy - origin[1]
        return dx * direction[0] + dy * direction[1], abs(
            dx * direction[1] - dy * direction[0]
        )

    along, off = axes(start, u)
    along_bar, off_bar = axes(center, v)
    corridor = (
        (along >= -1) & (along <= np.linalg.norm(center - start) + 1) & (off <= 1.5)
    )
    corridor |= (abs(along_bar) <= 6) & (off_bar <= 1.5)
    comp = label(mask & corridor, structure=np.ones((3, 3), bool))[0]
    starts = np.unique(comp[(abs(along) <= 1) & (off <= 1.5)])
    ends = np.unique(
        comp[(abs(along_bar) >= 3) & (abs(along_bar) <= 6) & (off_bar <= 1.5)]
    )
    return bool(set(starts[starts > 0]) & set(ends))


def at_hits(ranked, truth, goal):
    def prefixes(n):
        kept = defaultdict(list)
        for _, i, event in ranked[:n]:
            kept[i].append(event)
        return kept

    lo, hi = 0, len(ranked)
    while lo < hi:
        n = (lo + hi) // 2
        partial = prefixes(n)
        hits = sum(
            junction_counts(ref, positions(partial[i])) for i, ref in truth.items()
        )
        if hits >= goal:
            hi = n
        else:
            lo = n + 1
    return prefixes(lo)


def york(folder, outputs, report):
    from PIL import Image

    roster = read_csv(folder / "panel.csv")
    controls = defaultdict(list)
    for r in read_csv(EVIDENCE / "gap-review.csv"):
        if r["include"] == "1":
            controls[r["id"]].append(r)
    frozen = {
        r["id"]: r["pixel_sha256"]
        for r in read_csv(EVIDENCE / "photos.csv")
        if r["split"] == "fresh"
    }
    if {p["id"]: p["pixel_sha256"] for p in roster} != frozen or len(roster) != len(
        frozen
    ):
        raise ValueError("York evaluation requires the frozen 51-photo roster")
    archived = defaultdict(list)
    refs = defaultdict(list)
    for r in read_csv(EVIDENCE / "inputs.csv.gz"):
        if r["role"] == "reference":
            refs[r["id"], r["type"]].append(r)
        elif r["role"] == "context":
            archived[r["id"]].append(r)
    truth = {p["id"]: positions(refs[p["id"], "T"]) for p in roster}
    rankings = defaultdict(list)
    all_events, rows, gaps = {}, [], []
    for p in roster:
        i = p["id"]
        results = json.loads((outputs / (i + ".json")).read_text())
        original = archived[i]
        exact = [
            (float(r["x"]), float(r["y"]), r["type"], float(r["score"]))
            for r in original
        ]
        replay = [
            (r["x"], r["y"], r["type"], r["confidence"]) for r in results["context"]
        ]
        if exact != replay:
            raise ValueError(f"Context replay differs from frozen predictions for {i}")
        for method, events in results.items():
            all_events[i, method] = events
            row = dict(id=i, method=method)
            for kind in ("T", "X"):
                chosen = [r for r in events if r["type"] == kind]
                if kind == "T":
                    rankings[method].extend((r["confidence"], i, r) for r in chosen)
                ref = positions(refs[i, kind])
                if kind == "X":
                    budget = sum(r["type"] == "X" for r in results["context"])
                    ranked = sorted(chosen, key=lambda r: -r["confidence"])
                    row["X_at_context_count"] = junction_counts(
                        ref, positions(ranked[:budget])
                    )
                    row["X_capped_count"] = min(budget, len(ranked))
                row.update(
                    {
                        kind + "_ref": len(ref),
                        kind + "_pred": len(chosen),
                        kind + "_hit": junction_counts(ref, positions(chosen)),
                    }
                )
            rows.append(row)
    goal = min(sum(r["T_hit"] for r in rows if r["method"] == m) for m in rankings)
    for method, ranking in rankings.items():
        ranking.sort(key=lambda x: (-x[0], x[1], x[2]["x"], x[2]["y"]))
        equal = at_hits(ranking, truth, goal)
        for row in (r for r in rows if r["method"] == method):
            i = row["id"]
            events = all_events[i, method]
            clipped = equal[i] + [r for r in events if r["type"] == "X"]
            row["equal_T_pred"] = len(equal[i])
            row["equal_T_hit"] = junction_counts(truth[i], positions(equal[i]))
            if not controls[i]:
                continue
            lines = np.loadtxt(
                folder / (i + ".annot.csv"), delimiter=",", skiprows=1
            ).reshape(-1, 2, 2)
            mask = (
                np.asarray(Image.open(outputs / (i + "." + method + ".drawing.png")))
                > 0
            )
            for g in controls[i]:
                stem = lines[int(g["stem"])][:: 1 if g["reverse"] == "0" else -1]
                item = dict(id=i, token=g["token"], method=method)
                for label_, predictions in [("flag", events), ("equal_flag", clipped)]:
                    item[label_] = int(
                        (cdist(positions([g]), positions(predictions)) <= 4).any()
                    )
                item["connected"] = int(
                    drawings_bridge(mask, *stem, lines[int(g["target"])])
                )
                gaps.append(item)
    for m in rankings:
        assert sum(r["equal_T_hit"] for r in rows if r["method"] == m) == goal
    for name, records in (("york-photos.csv", rows), ("york-gaps.csv", gaps)):
        write_csv(report / name, records, list(records[0]))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("images", "outputs", "report"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    york(a.images, a.outputs, a.report)
