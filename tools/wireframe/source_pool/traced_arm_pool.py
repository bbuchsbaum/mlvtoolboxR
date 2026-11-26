"""Traced-branch four-arm opportunity with additional source curves in the pool.

A variant of tools/wireframe/traced_arm_opportunity.py. It reuses that
script's branch tracing, templates and checkpoint rule unchanged, and adds
the curves of each screened source to the combined FDoG and M-LSD pool.
With no added source it reproduces the original per-reference table.

    py tools/wireframe/source_pool/traced_arm_pool.py \
        .agent-work/toronto/staged .agent-work/contour-study/draw-toronto \
        .agent-work/crossing-groups/baseline .agent-work/source-pool/draw \
        .agent-work/source-pool/arms
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from arm_opportunity import ROSTER, maximum, read, sample  # noqa: E402
from score import read_csv, write_csv  # noqa: E402
from score_crossing_groups import _four_distinct_directions, four_arm_pairs  # noqa: E402
from traced_arm_opportunity import Curves, edges, templates  # noqa: E402

SOURCES = ("fdog512", "sed")
CONFIGS = {
    "combined": (),
    "with_fdog512": ("fdog512",),
    "with_sed": ("sed",),
    "with_both": ("fdog512", "sed"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("staged", "base", "baseline", "extra", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    rot = np.array([[1, 1], [-1, 1]]) / np.sqrt(2)
    roster = read_csv(ROSTER)
    assert len(roster) == 60 and {r["panel"] for r in roster} == {"development"}
    for row in roster:
        name = row["id"]
        a = read(args.staged / (name + ".json"))
        fd = read(args.base / (name + ".fdog.contours.json"))
        ml = read(args.base / (name + ".mlsd.contours.json"))
        pr = read(args.baseline / (name + ".curve_rank.contours.json"))
        extra = {k: read(args.extra / f"{name}.{k}.contours.json") for k in SOURCES}
        artist = {int(c["id"]): sample([c]) for c in a["polylines"]}
        drawings = [("fdog", fd), ("mlsd", ml), ("selected", pr)]
        drawings += list(extra.items())
        pools = {k: Curves(sample(draw["contours"])) for k, draw in drawings}
        strict = {i for draw in (fd, ml) for i, j in four_arm_pairs(a, draw)}
        declared = np.array(
            [
                [s["x"], s["y"]]
                for draw in (fd, ml)
                for s in draw["sites"]
                if s["type"] == "X"
            ]
        ).reshape(-1, 2)
        for i, s in enumerate(s for s in a["junctions"] if s["type"] == "X"):
            if (
                _four_distinct_directions(s["incident_directions"], observed=False)
                is None
            ):
                continue
            center = np.asarray(s["xy"])
            dirs = np.asarray(s["incident_directions"])
            local = [
                p for idx in set(s["incident_polyline_ids"]) for p in artist[int(idx)]
            ]
            br = Curves(local).branches(center)
            q = templates(center, dirs, br)
            rec = dict(
                id=name,
                category=row["category"],
                reference=i,
                usable=int(q is not None),
                source_strict=int(i in strict),
                declared=int(np.any(np.linalg.norm(declared - center, axis=1) <= 4)),
            )
            if q is not None:
                rec["artist"] = maximum(edges(q, br))
                assert rec["artist"] == 4
                branches = {k: p.branches(center) for k, p in pools.items()}
                es = {k: edges(q, b) for k, b in branches.items()}
                rec.update({k: maximum(x) for k, x in es.items()})
                for label, added in CONFIGS.items():
                    hits = es["fdog"] + es["mlsd"]
                    for k in added:
                        hits = hits + es[k]
                    rec[label] = maximum(hits)
                rec["rotated"] = maximum(
                    edges(
                        center + (q - center) @ rot, branches["fdog"] + branches["mlsd"]
                    )
                )
            rows.append(rec)
    fields = list(next(r for r in rows if r["usable"]))
    write_csv(args.output / "traced-pool-development.csv", rows, fields)
    t = [r for r in rows if r["usable"]]
    print("Artist gives four measurable branches", len(t), "/", len(rows))
    keys = ["artist", "fdog", "mlsd", "fdog512", "sed", "selected", "rotated"]
    print({k: sum(r[k] == 4 for r in t) for k in keys + list(CONFIGS)})
    for label in CONFIGS:
        hit = [r for r in t if r[label] == 4]
        gained = [r for r in hit if r["combined"] != 4]
        print(
            label,
            "four-arm",
            len(hit),
            "photos",
            len({r["id"] for r in hit}),
            "gained",
            len(gained),
            "on photos",
            len({r["id"] for r in gained}),
            "arms",
            dict(sorted(Counter(r[label] for r in t).items())),
        )


if __name__ == "__main__":
    main()
