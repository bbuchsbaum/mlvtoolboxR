"""Reference-guided development screen for locally available crossing branches."""

import argparse
import json
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

from score_crossing_groups import _four_distinct_directions, four_arm_pairs
from score import ROOT, read_csv, write_csv

ROSTER = ROOT / "docs/evidence/wireframe/curve-rank/development.csv"


def input_directories(description):
    parser = argparse.ArgumentParser(description=description)
    for name in ("staged", "base", "baseline", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    return args


radii = np.array([3, 5, 7, 9])
cos = np.cos(np.deg2rad(25))
rot = np.array([[1, 1], [-1, 1]]) / np.sqrt(2)


def sample(lines):
    paths = []
    for c in lines:
        if "points" in c:
            xy = np.array(c["points"])
            segs = np.stack((xy[:-1], xy[1:]), axis=1)
        else:
            segs = np.asarray(c["segments"]).reshape(-1, 2, 2)
        path = []
        last = None
        for a, b in segs:
            if last is None or np.linalg.norm(a - last) > 1e-7:
                if path:
                    paths.append(np.vstack(path))
                path = [a[None, :]]
            last = b
            n = int(np.ceil(np.linalg.norm(b - a) * 2))
            if n:
                path.append(a + (np.arange(1, n + 1) / n)[:, None] * (b - a))
        if path:
            paths.append(np.vstack(path))
    return paths


class Pool:
    def __init__(self, paths):
        self.paths = paths
        xy = np.vstack(paths) if paths else np.empty((0, 2))
        self.index = cKDTree(xy)
        self.ids = np.repeat(np.arange(len(paths)), [len(p) for p in paths])

    def edges(self, center, directions):
        result = []
        wanted = center + directions[:, None, :] * radii[None, :, None]
        for index in np.unique(self.ids[self.index.query_ball_point(center, 4)]):
            path = self.paths[index]
            d2 = np.sum((path - center) ** 2, axis=1)
            j = int(np.argmin(d2))
            for br in (path[j:], path[: j + 1][::-1]):
                far = np.flatnonzero(np.sum((br - center) ** 2, axis=1) > 11.01**2)
                if len(far):
                    br = br[: far[0]]
                if len(br) < 2:
                    continue
                # Traverse one continuous outward branch, matching 4 checkpoints in order.
                dist = ((br[None, :, None, :] - wanted[:, None, :, :]) ** 2).sum(
                    axis=-1
                )
                steps = br - br[0]
                length = np.linalg.norm(steps, axis=1)
                align = directions @ steps.T >= cos * length[None, :]
                hits = []
                for d in range(4):
                    start = 0
                    for k in range(4):
                        ok = dist[d, :, k] <= 4
                        if k == 3:
                            ok &= align[d] & (length > 0)
                        where = np.flatnonzero(ok[start:])
                        if not len(where):
                            break
                        start += where[0] + 1
                    else:
                        hits.append(d)
                if hits:
                    result.append(hits)
        return result


def maximum(columns):
    edges = np.zeros((4, len(columns)), bool)
    for j, hits in enumerate(columns):
        edges[hits, j] = True
    i, j = linear_sum_assignment(~edges)
    return int(edges[i, j].sum())


def check():
    center = np.array([14.0, 14.0])
    dirs = np.array([[1, 0], [-1, 0], [0, 1], [0, -1.0]])
    h = {"points": [center + [-10, 0], center + [10, 0]]}
    v = {"points": [center + [0, -10], center + [0, 10]]}
    t = {"points": [center, center + [0, 10]]}
    gap = {"points": [center + [0, -10], center + [0, -5]]}

    def count(paths, ds=dirs):
        return maximum(Pool(sample(paths)).edges(center, ds))

    assert (
        count([h, v]) == 4
        and count([h, t]) == 3
        and count([h, h]) == 2
        and count([h, t, gap]) == 3
        and count([h, v], dirs @ rot) == 0
    )


def read(p):
    return json.loads(p.read_text())


def main(args):
    check()
    records = []
    for row in read_csv(ROSTER):
        name = row["id"]
        artist = read(args.staged / (name + ".json"))
        fd = read(args.base / (name + ".fdog.contours.json"))
        ml = read(args.base / (name + ".mlsd.contours.json"))
        pred = read(args.baseline / (name + ".curve_rank.contours.json"))
        pools = {
            k: Pool(sample(raw))
            for k, raw in [
                ("fdog", fd["contours"]),
                ("mlsd", ml["contours"]),
                ("retained", pred["contours"]),
                ("artist", artist["polylines"]),
            ]
        }
        src = {i for draw in (fd, ml) for i, j in four_arm_pairs(artist, draw)}
        sites = {
            k: np.array(
                [[s["x"], s["y"]] for s in draw["sites"] if s["type"] == "X"]
            ).reshape(-1, 2)
            for k, draw in [("fdog", fd), ("mlsd", ml)]
        }
        for i, s in enumerate(s for s in artist["junctions"] if s["type"] == "X"):
            d = s["incident_directions"]
            rec = {
                "id": name,
                "category": row["category"],
                "artist_x_index": i,
                "eligible": int(
                    _four_distinct_directions(d, observed=False) is not None
                ),
            }
            if rec["eligible"]:
                dirs = np.asarray(d)
                dirs = dirs / np.linalg.norm(dirs, axis=1)[:, None]
                center = np.asarray(s["xy"])
                es = {k: p.edges(center, dirs) for k, p in pools.items()}
                rec.update({k: maximum(e) for k, e in es.items()})
                rec["combined"] = maximum(es["fdog"] + es["mlsd"])
                rec["rotated_combined"] = maximum(
                    pools["fdog"].edges(center, dirs @ rot)
                    + pools["mlsd"].edges(center, dirs @ rot)
                )
                rec["source_x_near"] = int(
                    any(
                        np.any(np.linalg.norm(xy - center, axis=1) <= 4)
                        for xy in sites.values()
                    )
                )
                rec["matched_source_x"] = int(i in src)
            records.append(rec)
    fields = list(next(r for r in records if r["eligible"]))
    write_csv(args.output / "development-probe.csv", records, fields)
    el = [r for r in records if r["eligible"]]
    print("Eligible / total", len(el), len(records))
    print(
        "Four supported arms",
        {
            k: sum(r[k] == 4 for r in el)
            for k in [
                "artist",
                "fdog",
                "mlsd",
                "combined",
                "retained",
                "rotated_combined",
            ]
        },
    )
    print(
        "Paired opportunities",
        Counter((r["combined"] == 4, r["retained"] == 4) for r in el),
    )
    print(
        "Declared site within 4px / source strict among supported",
        Counter(
            (r["source_x_near"], r["matched_source_x"])
            for r in el
            if r["combined"] == 4
        ),
    )
    print(
        "Photo counts",
        len({r["id"] for r in el if r["combined"] == 4}),
        len({r["id"] for r in el if r["retained"] == 4}),
    )
    print(
        "Real vs rotated",
        Counter((r["combined"] == 4, r["rotated_combined"] == 4) for r in el),
    )
    print(
        "Strict sources covered?",
        Counter(r["combined"] for r in el if r["matched_source_x"]),
    )
    print(
        "Corrected for ideal geometry subset",
        Counter((r["artist"] == 4, r["combined"] == 4) for r in el),
    )


if __name__ == "__main__":
    main(input_directories(__doc__))
