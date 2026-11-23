"""Fit and compare fixed-location arm verifiers on prepared photo directories."""
import argparse
import gzip
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image, ImageOps
from sklearn.ensemble import HistGradientBoostingClassifier
from threadpoolctl import threadpool_limits

from arm_verification import baseline_features, exchange, measurements
from check_localization import cardinality
from fit_topology import pack
from score import ROOT, match, positions, read_csv, references, write_csv
from topology import probabilities


def stage(split, images, output):
    """Recreate a split from the frozen sites, scores and user-supplied photographs."""
    evidence = ROOT / "docs/evidence/wireframe/native-arms"
    photos = [r for r in read_csv(evidence / "photos.csv") if r["split"] == split]
    if not photos:
        raise ValueError(f"Unknown split {split}")
    inputs = defaultdict(list)
    for r in read_csv(evidence / "inputs.csv.gz"):
        inputs[r["id"], r["role"]].append(r)
    output.mkdir(parents=True, exist_ok=False)
    for r in photos:
        i = r["id"]
        for suffix, key in ((".png", "native_sha256"), (".small.png", "pixel_sha256")):
            src = images / (i + suffix)
            with Image.open(src) as image:
                pixels = ImageOps.exif_transpose(image).convert("RGB").tobytes()
            if hashlib.sha256(pixels).hexdigest() != r[key]:
                raise ValueError(f"Pixel checksum mismatch at {src}")
            shutil.copyfile(src, output / src.name)
        lines = images / (i + ".annot.csv")
        if lines.exists():
            shutil.copyfile(lines, output / lines.name)
        for role, suffix, fields in (
            ("context", ".events.csv", ["x", "y", "type", "score"]),
            ("reference", ".sites.csv", ["x", "y", "type"]),
        ):
            rows = [{k: v[k] for k in fields} for v in inputs[i, role]]
            write_csv(output / (i + suffix), rows, fields)
    write_csv(output / "panel.csv", photos, list(photos[0]))


def sites(folder, photo):
    name = photo["id"]
    exact = folder / (name + ".sites.csv")
    if exact.exists():
        rows = read_csv(exact)
        return [positions([r for r in rows if r["type"] == k]) for k in "TX"]
    lines = np.loadtxt(
        folder / (name + ".annot.csv"), skiprows=1, delimiter=","
    ).reshape(-1, 2, 2)
    truth = references(lines, [int(photo["width"]), int(photo["height"])])
    return [
        np.asarray([p for k, p in truth if k == kind]).reshape(-1, 2) for kind in "TX"
    ]


def photo_data(folder, photo):
    name = photo["id"]
    events = read_csv(folder / (name + ".events.csv"))
    xy = positions(events)
    size = np.array([int(photo["width"]), int(photo["height"])])
    inputs = {"scores": baseline_features(events)}
    for mode, suffix in (("working", ".small.png"), ("native", ".png")):
        with Image.open(folder / (name + suffix)) as image:
            pixels = np.asarray(ImageOps.exif_transpose(image).convert("RGB"))
        if (
            mode == "working"
            and hashlib.sha256(pixels.tobytes()).hexdigest() != photo["pixel_sha256"]
        ):
            raise ValueError(f"Working pixels differ for {name}")
        inputs[mode] = np.column_stack(
            [inputs["scores"], measurements(pixels, xy, size)]
        )
    return xy, inputs, sites(folder, photo)


def training(folders, output):
    # Ambiguous (both) and unmatched (neither) locations provide no T-versus-X label.
    xs = {m: [] for m in ("scores", "working", "native")}
    ys = []
    for folder in folders:
        for row in read_csv(folder / "panel.csv"):
            xy, inputs, refs = photo_data(folder, row)
            labels = np.zeros((len(xy), 2), bool)
            for j, ref in enumerate(refs):
                for _, k in match(ref, xy):
                    labels[k, j] = True
            unique = labels.sum(1) == 1
            ys.extend(labels[unique, 1])
            for m in xs:
                xs[m].append(inputs[m][unique])
    fitted = {}
    for m, pieces in xs.items():
        x = np.vstack(pieces)
        clf = HistGradientBoostingClassifier(
            max_iter=65,
            max_leaf_nodes=7,
            learning_rate=0.06,
            l2_regularization=4,
            min_samples_leaf=30,
            early_stopping=False,
            random_state=913,
        )
        with threadpool_limits(limits=2):
            clf.fit(x, ys)
            fitted[m] = pack(clf)
            np.testing.assert_allclose(
                probabilities(x, [fitted[m]])[:, 0],
                clf.predict_proba(x)[:, 1],
                atol=1e-12,
            )
    output.write_bytes(
        gzip.compress(json.dumps(fitted, separators=(",", ":")).encode(), mtime=0)
    )
    print(
        "Unambiguous matched training locations, X, T",
        len(ys),
        sum(ys),
        len(ys) - sum(ys),
        flush=True,
    )


def evaluate(folder, model, output):
    fitted = json.loads(gzip.decompress(model.read_bytes()))
    scores, decisions = [], []
    for photo in read_csv(folder / "panel.csv"):
        i = photo["id"]
        xy, inputs, refs = photo_data(folder, photo)
        is_x = inputs["scores"][:, 0].astype(bool)
        kinds = {"context": is_x}
        for m, tree in fitted.items():
            kinds[m] = exchange(is_x, probabilities(inputs[m], [tree])[:, 0])
        for m, assigned in kinds.items():
            if assigned.sum() != is_x.sum():
                raise ValueError(f"Type counts changed on {i}")
            row = dict(
                id=i,
                method=m,
                n_T=int((~assigned).sum()),
                n_X=int(assigned.sum()),
                swaps=int(np.sum(assigned & ~is_x)),
            )
            for j, (k, idx) in enumerate((("T", ~assigned), ("X", assigned))):
                c = cardinality(refs[j], xy[idx])
                if c != len(match(refs[j], xy[idx])):
                    raise ValueError(f"Matching disagreement on {i}")
                row[k] = c
                row[k + "_reference"] = len(refs[j])
            scores.append(row)
        for j, (x, y) in enumerate(xy):
            decisions.append(
                dict(
                    id=i,
                    x=x,
                    y=y,
                    **{m: "X" if labels[j] else "T" for m, labels in kinds.items()},
                )
            )
    write_csv(output.with_suffix(".csv"), scores, list(scores[0]))
    write_csv(
        output.with_suffix(".sites.csv"),
        decisions,
        ["id", "x", "y", "context", "scores", "working", "native"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("stage")
    p.add_argument("split")
    p.add_argument("images", type=Path)
    p.add_argument("output", type=Path)
    p = sub.add_parser("fit")
    p.add_argument("output", type=Path)
    p.add_argument("folders", nargs="+", type=Path)
    p = sub.add_parser("score")
    for field in ("folder", "model", "output"):
        p.add_argument(field, type=Path)
    a = parser.parse_args()
    if a.command == "stage":
        stage(a.split, a.images, a.output)
    elif a.command == "fit":
        training(a.folders, a.output)
    else:
        evaluate(a.folder, a.model, a.output)
