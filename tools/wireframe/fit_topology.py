"""Fit photographic junction scores on the official training split only."""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.spatial import cKDTree

from score import read_csv, references
from topology import cross, model_features, probabilities


def negative_sites(lines, refs):
    s = lines[:, 0]
    e = lines[:, 1]
    v = e - s
    leng = np.linalg.norm(v, axis=1)
    # Near-miss projections strictly outside one segment and inside the second.
    dirs = v / np.maximum(1e-9, leng[:, None])
    a, b = np.nonzero(np.ones((len(s), len(s)), bool))
    den = cross(dirs[a], dirs[b])
    keep = (abs(den) >= 0.65) & (leng[a] >= 10) & (leng[b] >= 20)
    a, b, den = a[keep], b[keep], den[keep]
    diff = s[b] - s[a]
    ua = cross(diff, dirs[b]) / den
    ub = cross(diff, dirs[a]) / den
    beyond = np.maximum(-ua, ua - leng[a])
    keep = (beyond >= 4) & (beyond <= 10) & (ub >= 8) & (ub <= leng[b] - 8)
    points = s[a[keep]] + ua[keep, None] * dirs[a[keep]]
    if len(refs) and len(points):
        points = points[cKDTree([p for _, p in refs]).query(points)[0] >= 8]
    return points


def feature_file(folder, i, variant):
    path = folder / (i + "." + variant + ".npz")
    if not path.exists():
        with Image.open(folder / (i + ".small.png")) as photo:
            rgb = np.asarray(photo.convert("RGB"))
        wire = np.loadtxt(
            folder / (i + ".hawp.csv"), delimiter=",", skiprows=1
        ).reshape(-1, 5)
        detail = (
            np.loadtxt(
                folder / (i + ".hawp-tiles.csv"), delimiter=",", skiprows=1
            ).reshape(-1, 5)
            if variant == "multiscale"
            else None
        )
        xy, f = model_features(rgb, wire, variant, detail)
        np.savez_compressed(path, xy=xy, f=f)
    return path


def cache(folder, variant, training=False):
    for row in read_csv(folder / "panel.csv"):
        i = row["id"]
        path = feature_file(folder, i, variant)
        if training:
            with np.load(path) as record:
                xy, f = record["xy"], record["f"]
            lines = np.loadtxt(
                folder / (i + ".annot.csv"), delimiter=",", skiprows=1
            ).reshape(-1, 2, 2)
            refs = references(lines, [int(row["width"]), int(row["height"])])
            labels = []
            for kind in ("T", "X"):
                pts = [p for k, p in refs if k == kind]
                dist = cKDTree(pts).query(xy)[0] if pts else np.full(len(xy), np.inf)
                labels.append(dist <= 3.5)
            neg = negative_sites(lines, refs) if variant != "arms" else []
            hard = (
                cKDTree(neg).query(xy)[0] <= 4 if len(neg) else np.zeros(len(xy), bool)
            )
            yield f, np.column_stack(labels), hard


def pack(model):
    # Store only numerical splits and leaf values. Inference needs neither a
    # pickle nor scikit-learn private classes.
    trees = []
    for (predictor,) in model._predictors:
        nodes = predictor.nodes
        trees.append(
            {
                k: nodes[col].tolist()
                for k, col in [
                    ("feature", "feature_idx"),
                    ("threshold", "num_threshold"),
                    ("left", "left"),
                    ("right", "right"),
                    ("leaf", "is_leaf"),
                    ("value", "value"),
                ]
            }
        )
    return dict(bias=model._baseline_prediction.item(), trees=trees)


def fit(folder, output, variant):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from threadpoolctl import threadpool_limits

    if {r["panel"] for r in read_csv(folder / "panel.csv")} != {"train"}:
        raise ValueError("Training requires the official training roster")
    records = list(cache(folder, variant, training=True))
    f, y, hard = [np.concatenate(values) for values in zip(*records)]
    models = []
    with threadpool_limits(limits=4):
        for labels in y.T:
            model = HistGradientBoostingClassifier(
                max_iter=100 if variant == "arms" else 130,
                max_leaf_nodes=15,
                learning_rate=0.09,
                l2_regularization=2,
                min_samples_leaf=45,
                early_stopping=False,
                random_state=7311,
            )
            model.fit(
                f,
                labels,
                sample_weight=None
                if variant == "arms"
                else np.where(hard & ~labels, 8.0, 1.0),
            )
            packed = pack(model)
            actual = probabilities(f[::23], [packed])[:, 0]
            np.testing.assert_allclose(
                actual, model.predict_proba(f[::23])[:, 1], atol=1e-12
            )
            models.append(packed)
    output.write_text(json.dumps(models, separators=(",", ":")) + "\n")
    print(
        "training photos, rows, type-positive rows",
        len(records),
        len(f),
        y.sum(0),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--fit", type=Path, help="Fit and write a compact JSON model")
    parser.add_argument(
        "--variant", choices=("arms", "photoarms", "multiscale"), default="multiscale"
    )
    args = parser.parse_args()
    if args.fit:
        fit(args.folder, args.fit, args.variant)
    else:
        list(cache(args.folder, args.variant))
