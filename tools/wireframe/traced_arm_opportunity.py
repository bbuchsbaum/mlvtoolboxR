"""Trace the artist's actual incident branches and pair them to source branches."""

from arm_opportunity import ROSTER, cos, input_directories, maximum, radii, read, sample
from score import read_csv, write_csv
from score_crossing_groups import _four_distinct_directions, four_arm_pairs
from scipy.spatial import cKDTree
from scipy.optimize import linear_sum_assignment
import numpy as np
from collections import Counter


class Curves:
    def __init__(self, path_list):
        self.paths = path_list
        self.index = cKDTree(np.vstack(path_list) if path_list else np.empty((0, 2)))
        self.ids = np.repeat(np.arange(len(path_list)), [len(p) for p in path_list])

    def branches(self, center):
        result = []
        for idx in np.unique(self.ids[self.index.query_ball_point(center, 4)]):
            path = self.paths[idx]
            m = np.sum((path - center) ** 2, axis=1).argmin()
            for br in (path[m:], path[: m + 1][::-1]):
                # Only inspect near-centre portions of curves.
                outside = np.flatnonzero(np.linalg.norm(br - center, axis=1) > 13.01)
                if len(outside):
                    br = br[: outside[0]]
                if len(br) > 1:
                    result.append(br)
        return result


def templates(center, dirs, branches):
    pts = []
    rays = []
    for br in branches:
        along = np.r_[0, np.linalg.norm(np.diff(br, axis=0), axis=1).cumsum()]
        if along[-1] < 9:
            continue
        q = np.column_stack([np.interp(radii, along, axis) for axis in br.T])
        initial = q[0] - br[0]
        n = np.linalg.norm(initial)
        if n < 1.0e-8:
            continue
        pts.append(q)
        rays.append(initial / n)
    if len(pts) < 4:
        return None
    sim = dirs @ np.asarray(rays).T
    i, j = linear_sum_assignment(np.where(sim >= cos, -sim, 100))
    return np.array([pts[k] for k in j]) if np.all(sim[i, j] >= cos) else None


def edges(target, bs):
    # Compare directions along the traced arm between the first and last checkpoints.
    net = target[:, -1] - target[:, 0]
    length = np.linalg.norm(net, axis=1)
    result = []
    for b in bs:
        eligible = []
        delta = b[None, :, None, :] - target[:, None, :, :]
        close = np.sum(delta**2, axis=-1) <= 4 + 1e-12
        for d in range(4):
            for i in np.flatnonzero(close[d, :, 0]):
                starts = np.flatnonzero(close[d, i + 1 :, 1])
                if not len(starts):
                    continue
                j = i + 1 + starts[0]
                ks = np.flatnonzero(close[d, j + 1 :, 2])
                if not len(ks):
                    continue
                k = j + 1 + ks[0]
                finals = b[k + 1 :][close[d, k + 1 :, 3]] - b[i]
                norm = np.linalg.norm(finals, axis=1)
                if np.any(
                    (norm > 1e-7)
                    & (length[d] > 1e-7)
                    & (finals @ net[d] >= cos * norm * length[d])
                ):
                    eligible.append(d)
                    break
        if eligible:
            result.append(eligible)
    return result


def main(args):
    rows = []
    rot = np.array([[1, 1], [-1, 1]]) / np.sqrt(2)
    for row in read_csv(ROSTER):
        name = row["id"]
        a = read(args.staged / (name + ".json"))
        fd = read(args.base / (name + ".fdog.contours.json"))
        ml = read(args.base / (name + ".mlsd.contours.json"))
        pr = read(args.baseline / (name + ".curve_rank.contours.json"))
        artist = {int(c["id"]): sample([c]) for c in a["polylines"]}
        pools = {
            k: Curves(sample(draw["contours"]))
            for k, draw in [("fdog", fd), ("mlsd", ml), ("selected", pr)]
        }
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
                rec["combined"] = maximum(es["fdog"] + es["mlsd"])
                rec["rotated"] = maximum(
                    edges(
                        center + (q - center) @ rot, branches["fdog"] + branches["mlsd"]
                    )
                )
            rows.append(rec)
    fields = list(next(r for r in rows if r["usable"]))
    write_csv(args.output / "traced-development.csv", rows, fields)
    t = [r for r in rows if r["usable"]]
    print("Artist gives four measurable branches", len(t), "/", len(rows))
    print(
        {
            k: sum(r[k] == 4 for r in t)
            for k in ["artist", "fdog", "mlsd", "combined", "rotated", "selected"]
        }
    )
    print("Raw arms", Counter(r["combined"] for r in t))
    print("pool only across photos", len({r["id"] for r in t if r["combined"] == 4}))
    print(
        "Declared vs strict for complete source opportunities",
        Counter((r["declared"], r["source_strict"]) for r in t if r["combined"] == 4),
    )
    print(
        "Source strict usable each / all",
        [(r["usable"], r.get("combined")) for r in rows if r["source_strict"]],
    )
    print(
        "Source 4 vs selected",
        Counter((r["combined"] == 4, r["selected"] == 4) for r in t),
    )


if __name__ == "__main__":
    main(input_directories(__doc__))
