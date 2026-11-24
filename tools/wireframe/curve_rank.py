"""Fit and apply a photo-only marginal-utility model for complete vector paths."""

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from sklearn.ensemble import HistGradientBoostingRegressor
from threadpoolctl import threadpool_limits

from curve_observations import COLUMN_NAMES, observations
from fit_topology import pack
from topology import tree_response
from score import ROOT, read_csv, write_csv

EVIDENCE = ROOT / "docs/evidence/wireframe/curve-rank"
SOURCES = ("fdog", "mlsd")
JUNCTION_CREDIT = 8.0
FEATURE_NAMES = (
    *COLUMN_NAMES,
    "fraction_fresh_ink",
    "fraction_fresh_neighborhood",
    "fresh_neighborhood_per_ink",
    "endpoints_near_selected",
    "source_sites",
    "source_site_progress",
    "source_sites_completable",
    "budget_remaining_fraction",
)


class FrozenModel:
    def __init__(self, artifact):
        if tuple(artifact["features"]) != FEATURE_NAMES:
            raise ValueError("The fitted feature order differs from the observations")
        self.ensemble = artifact["ensemble"]
        self.ensemble["trees"] = [
            {key: np.asarray(values) for key, values in tree.items()}
            for tree in self.ensemble["trees"]
        ]

    def predict(self, features):
        return tree_response(features, self.ensemble)


def mask(path):
    with Image.open(path) as image:
        return np.asarray(image) > 0


def tube(pixels, shape):
    """Indices within two Euclidean pixels of a path, in R matrix order."""
    h, w = shape
    y, x = pixels % h, pixels // h
    pieces = []
    for dy in range(-2, 3):
        for dx in range(-2, 3):
            if dx * dx + dy * dy <= 4:
                yy, xx = y + dy, x + dx
                ok = (yy >= 0) & (yy < h) & (xx >= 0) & (xx < w)
                pieces.append(yy[ok] + h * xx[ok])
    return np.unique(np.concatenate(pieces))


def membership(groups, n):
    counts = np.array([len(g) for g in groups], dtype=int)
    indices = np.concatenate(groups) if groups else np.empty(0, int)
    return csr_matrix(
        (np.ones(len(indices), dtype=np.int32), indices, np.r_[0, np.cumsum(counts)]),
        shape=(len(groups), n),
    ).tocsc()


@dataclass
class Pool:
    shape: tuple
    paths: list
    pixels: list
    tubes: list
    ends: np.ndarray
    obs: np.ndarray
    junctions: list
    members: object
    ink_index: object
    area_index: object
    budget: int


def read_pool(staged, drawings, raster, row):
    name = row["id"]
    with Image.open(staged / (name + ".small.png")) as img:
        rgb = np.asarray(img.convert("RGB"))
    if hashlib.sha256(rgb.tobytes()).hexdigest() != row["pixel_sha256"]:
        raise ValueError(f"Photo checksum differs for {name}")
    shape = rgb.shape[:2]
    h, w = shape
    inputs = json.loads((raster / (name + ".json")).read_text())
    paths, pixels, sources, groups, sites = [], [], [], [], []
    for source in SOURCES:
        prefix = drawings / f"{name}.{source}"
        drawing = json.loads(Path(str(prefix) + ".contours.json").read_text())
        ids = {c["id"]: len(paths) + j for j, c in enumerate(drawing["contours"])}
        full = np.zeros(h * w, bool)
        for c, p in zip(drawing["contours"], inputs[source], strict=True):
            path = np.asarray(c["segments"], float).reshape(-1, 4)
            ix = np.unique(np.asarray(p, int))
            if not len(path) or not len(ix) or ix.min() < 0 or ix.max() >= h * w:
                raise ValueError(f"Invalid contour in {prefix}")
            full[ix] = True
            paths.append(path)
            pixels.append(ix)
            sources.append(source)
        if not np.array_equal(
            full, mask(Path(str(prefix) + ".drawing.png")).ravel(order="F")
        ):
            raise ValueError(f"Path raster differs from package drawing for {prefix}")
        if source == "fdog":
            budget = int(full.sum())
        for site in drawing["sites"]:
            members = sorted({ids[m["contour_id"]] for m in site["members"]})
            if members and site["type"] in ("T", "X"):
                groups.append(np.array(members, int))
                sites.append((site["type"], site["x"], site["y"]))
    ends = np.array([[p[0, :2], p[-1, 2:]] for p in paths]).reshape(-1, 2, 2)
    xy = np.clip(np.rint(ends - 1).astype(int), [0, 0], [w - 1, h - 1])
    tubes = [tube(p, shape) for p in pixels]
    obs = observations(rgb, pixels, paths, sources)
    if not np.all(np.isfinite(obs)):
        raise ValueError("Photo-only observations must be finite")
    return Pool(
        shape,
        paths,
        pixels,
        tubes,
        xy[:, :, 1] + h * xy[:, :, 0],
        obs,
        sites,
        membership(groups, len(paths)).T.tocsr(),
        membership(pixels, w * h),
        membership(tubes, w * h),
        budget,
    )


def label_junctions(pool, annotations):
    """Training credit for source sites near a same-type artist reference."""
    kinds = np.array([s[0] for s in pool.junctions])
    found = np.zeros(len(kinds), bool)
    for kind in ("T", "X"):
        truth = [s["xy"] for s in annotations["junctions"] if s["type"] == kind]
        indices = np.flatnonzero(kinds == kind)
        if len(truth) and len(indices):
            xy = [pool.junctions[j][1:] for j in indices]
            found[indices] = cKDTree(truth).query(xy)[0] <= 4
    return found


class Selection:
    def __init__(self, pool, artist=None, matched=None):
        self.pool = pool
        self.order = []
        self.remaining = pool.budget
        self.used = np.zeros(np.prod(pool.shape), bool)
        self.covered = np.zeros_like(self.used)
        self.ink = np.array([len(p) for p in pool.pixels])
        self.area = np.array([len(p) for p in pool.tubes])
        self.initial_ink = self.ink.copy()
        self.initial_area = self.area.copy()
        self.degree = np.asarray(pool.members.sum(axis=0)).ravel()
        self.missing = self.degree.copy()
        self.artist = artist
        self.matched = matched
        if artist is not None:
            self.gain = np.array([artist[p].sum() for p in pool.tubes], float)
            self.member_credit = np.asarray(pool.members @ (matched / self.degree))

    def eligible(self):
        return (self.ink > 0) & (self.ink <= self.remaining)

    def scale(self):
        return np.maximum(self.area, self.ink)

    def features(self):
        touches = self.covered[self.pool.ends].sum(axis=1)
        neighbors = self.pool.members
        progress = 1 - self.missing / self.degree
        present = np.asarray(neighbors @ progress)
        completion = np.asarray(neighbors @ (self.missing == 1).astype(float))
        return np.column_stack(
            [
                self.pool.obs,
                self.ink / self.initial_ink,
                self.area / self.initial_area,
                self.area / np.maximum(self.ink, 1),
                touches,
                np.asarray(neighbors.sum(axis=1)).ravel(),
                present,
                completion,
                np.full(len(self.ink), self.remaining / max(self.pool.budget, 1)),
            ]
        )

    def reward(self):
        completion = self.pool.members @ (self.matched & (self.missing == 1)).astype(
            float
        )
        return self.gain + JUNCTION_CREDIT / 2 * (self.member_credit + completion)

    def take(self, j, *, allow_covered=False):
        if not self.eligible()[j] and not (
            allow_covered and self.ink[j] == 0 and j not in self.order
        ):
            raise ValueError("Chosen path is outside the remaining ink budget")
        p = self.pool
        ink = p.pixels[j][~self.used[p.pixels[j]]]
        area = p.tubes[j][~self.covered[p.tubes[j]]]
        self.remaining -= len(ink)
        self.used[ink] = True
        self.covered[area] = True
        self.ink -= np.bincount(p.ink_index[:, ink].indices, minlength=len(self.ink))
        self.area -= np.bincount(
            p.area_index[:, area].indices, minlength=len(self.area)
        )
        if self.artist is not None:
            hit = area[self.artist[area]]
            self.gain -= np.bincount(
                p.area_index[:, hit].indices, minlength=len(self.gain)
            )
        self.missing[p.members[j].indices] -= 1
        self.order.append(int(j))


def best(state, rewards, valid):
    candidates = np.flatnonzero(valid)
    values = rewards[candidates] / state.ink[candidates]
    tied = candidates[values == values.max()]
    return tied[np.argmax(rewards[tied])]


def training_rows(pool, artist, matched, photo_id):
    chunks = []
    labels = []
    weights = []
    seed = int.from_bytes(hashlib.sha256(photo_id.encode()).digest()[:8], "little")
    perm = np.random.default_rng(seed).permutation(len(pool.paths))
    for policy in ("oracle", "shuffle"):
        state = Selection(pool, artist, matched)
        marker = 0 if policy == "oracle" else 1
        while state.eligible().any() and marker < 5:
            valid = state.eligible()
            reward = state.reward()
            fraction = 1 - state.remaining / pool.budget
            if fraction + 1e-12 >= marker / 5:
                chunks.append(state.features()[valid])
                labels.append(reward[valid] / state.scale()[valid])
                weights.append(np.sqrt(state.ink[valid]))
                marker = int(fraction * 5 + 1e-9) + 1
                if marker >= 5:
                    break
            chosen = (
                best(state, reward, valid)
                if policy == "oracle"
                else perm[valid[perm]][0]
            )
            state.take(chosen)
    if not chunks:
        return None
    w = np.concatenate(weights)
    return np.vstack(chunks), np.concatenate(labels), w / w.sum()


def predict_paths(pool, model):
    """Select with source geometry and RGB observations; no artist inputs."""
    return fill_singles(Selection(pool), model)


def fill_singles(state, model):
    """Spend the remaining budget on the best feasible paths."""
    while state.eligible().any():
        valid = state.eligible()
        reward = np.zeros(len(valid))
        reward[valid] = (
            np.maximum(0, model.predict(state.features()[valid])) * state.scale()[valid]
        )
        state.take(best(state, reward, valid))
    return state


def fit(records):
    data = [r for r in records if r is not None]
    x, y, w = (np.concatenate(v) for v in zip(*data, strict=True))
    model = HistGradientBoostingRegressor(
        max_iter=80,
        max_leaf_nodes=12,
        min_samples_leaf=50,
        learning_rate=0.065,
        l2_regularization=2,
        early_stopping=False,
        random_state=20261124,
    )
    model.fit(x, y, sample_weight=w / w.mean())
    return model


def emit(pool, state, row, output):
    name = row["id"]
    with (output / (name + ".curve_rank.segments.csv")).open("w") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["group", "x1", "y1", "x2", "y2"])
        for group, idx in enumerate(state.order, 1):
            writer.writerows([group, *s] for s in pool.paths[idx])
    return dict(
        id=name,
        fold=row["fold"],
        paths=len(pool.paths),
        selected=len(state.order),
        budget=pool.budget,
        ink=int(state.used.sum()),
        segments_sha256=hashlib.sha256(
            (output / (name + ".curve_rank.segments.csv")).read_bytes()
        ).hexdigest(),
    )


def roster(name):
    rows = read_csv(EVIDENCE / (name + ".csv"))
    if len({r["id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate photograph")
    return rows


def experiment(args):
    advancing = args.command in ("refit", "predict")
    if advancing:
        gate = json.loads(args.gate.read_text())
        if not gate["gate_pass"] or gate["values"]["curve_rank"]["photos"] != 60:
            raise ValueError("The frozen development gate did not pass")
    rows = roster("reserve-test" if args.command == "predict" else "development")
    pools = {
        r["id"]: read_pool(args.staged, args.drawings, args.raster, r) for r in rows
    }
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = []
    if args.command == "predict":
        model = FrozenModel(json.loads(args.model.read_text()))
        for r in rows:
            p = pools[r["id"]]
            summaries.append(emit(p, predict_paths(p, model), r, args.output))
    else:
        records = {}
        for r in rows:
            name = r["id"]
            p = pools[name]
            artist = mask(args.drawings / (name + ".artist.drawing.png")).ravel(
                order="F"
            )
            annotations = json.loads((args.staged / (name + ".json")).read_text())
            records[name] = training_rows(
                p, artist, label_junctions(p, annotations), name
            )
        if args.command == "refit":
            model = fit(list(records.values()))
            artifact = dict(
                features=FEATURE_NAMES,
                ensemble=pack(model),
                development_sha256=hashlib.sha256(
                    (EVIDENCE / "development.csv").read_bytes()
                ).hexdigest(),
            )
            encoded = (
                json.dumps(artifact, separators=(",", ":"), allow_nan=False) + "\n"
            )
            frozen = FrozenModel(json.loads(encoded))
            training_features = np.concatenate(
                [record[0] for record in records.values() if record is not None]
            )
            np.testing.assert_allclose(
                frozen.predict(training_features),
                model.predict(training_features),
                rtol=0,
                atol=1e-14,
            )
            args.model.write_text(encoded)
            return
        for f in range(1, 6):
            training = [r for r in rows if int(r["fold"]) != f]
            held = [r for r in rows if int(r["fold"]) == f]
            if (
                len(training) != 48
                or len(held) != 12
                or len({r["category"] for r in held}) != 6
            ):
                raise ValueError("The frozen five-fold roster has changed")
            model = fit([records[r["id"]] for r in training])
            for r in held:
                p = pools[r["id"]]
                summaries.append(emit(p, predict_paths(p, model), r, args.output))
            print(f"Completed fold {f} with {len(held)} held-out photos", flush=True)
    write_csv(args.output / "selection.csv", summaries, list(summaries[0]))
    write_csv(
        args.output / "r-plan.csv",
        [dict(id=r["id"], method="curve_rank") for r in rows],
        ["id", "method"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["crossfit", "refit", "predict"])
    for name in ("staged", "drawings", "raster", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument(
        "--model",
        type=Path,
        help="Numeric JSON fit artifact, required for refit or prediction",
    )
    parser.add_argument(
        "--gate",
        type=Path,
        help="Development gate JSON, required for refit or prediction",
    )
    args = parser.parse_args()
    if args.command != "crossfit" and (args.gate is None or args.model is None):
        parser.error("Refit and prediction require --gate and --model")
    with threadpool_limits(limits=2):
        experiment(args)
