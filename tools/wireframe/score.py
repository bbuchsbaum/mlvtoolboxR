"""Match junction types derived from the human line annotations to detected junctions."""
import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

METHODS = (
    "fdog0",
    "fdog6",
    "mlsd_mask",
    "mlsd_vector",
    "mlsd_clahe",
    "annotation_mask",
    "annotation_vector",
)


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def classify(rays):
    angles = []
    for x, y in rays:
        angle = np.degrees(np.arctan2(y, x)) % 360
        if all(abs((angle - other + 180) % 360 - 180) > 10 for other in angles):
            angles.append(angle)
    angles.sort()
    if len(angles) not in (3, 4) or min(np.diff(angles + [angles[0] + 360])) < 30:
        return None
    if len(angles) == 3:
        if any(abs(b - a - 180) <= 20 for a, b in itertools.combinations(angles, 2)):
            return "T"
    elif all(abs(angles[k + 2] - angles[k] - 180) <= 20 for k in (0, 1)):
        return "X"
    return None


def references(lines, size):
    points = list(lines.reshape(-1, 2))
    origins = lines[:, 0]
    vectors = lines[:, 1] - origins
    lengths = np.linalg.norm(vectors, axis=1)
    unit = vectors / np.maximum(lengths, 1e-9)[:, None]
    for i, j in itertools.combinations(range(len(lines)), 2):
        den = cross(vectors[i], vectors[j])
        if abs(den) < 1e-9:
            continue
        delta = origins[j] - origins[i]
        t, u = cross(delta, vectors[j]) / den, cross(delta, vectors[i]) / den
        if 0 < t < 1 and 0 < u < 1:
            points.append(origins[i] + t * vectors[i])
    accepted = []
    for point in points:
        if np.any(point < 5) or np.any(point > np.array(size) - 4):
            continue
        along = ((point - origins) * unit).sum(axis=1)
        projected = origins + along[:, None] * unit
        touches = (
            (np.linalg.norm(projected - point, axis=1) <= 0.5)
            & (along >= -0.5)
            & (along <= lengths + 0.5)
            & (lengths > 0)
        )
        if sum(touches) < 2:
            continue
        rays = np.concatenate(
            [
                -vectors[touches & (along >= 4)],
                vectors[touches & (lengths - along >= 4)],
            ]
        )
        kind = classify(rays)
        if kind and not any(
            k == kind and np.linalg.norm(point - q) <= 1 for k, q in accepted
        ):
            accepted.append((kind, point))
    return accepted


def match(reference, predicted):
    if not len(reference) or not len(predicted):
        return []
    distances = cdist(reference, predicted)
    # A forbidden pair costs more than the sum of all feasible distances.
    costs = np.where(distances <= 4, distances, 5 * (min(distances.shape) + 1))
    refs, preds = linear_sum_assignment(costs)
    return [(i, j) for i, j in zip(refs, preds) if distances[i, j] <= 4]


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def supplement(folder):
    panel = [
        p for p in read_csv(folder / "panel.csv") if p["panel"] in ("sample", "witness")
    ]
    photos = []
    for photo in panel:
        lines = np.loadtxt(
            folder / (photo["id"] + ".annot.csv"), delimiter=",", skiprows=1
        ).reshape(-1, 2, 2)
        refs = references(lines, [int(photo["width"]), int(photo["height"])])
        photos.append((photo["id"], photo["panel"], refs))
    return photos


def original(folder, manifest):
    prepared = {p["id"]: p["pixel_sha256"] for p in read_csv(folder / "panel.csv")}
    if manifest.is_file():
        photos, reviews = [], []
        for photo in json.loads(manifest.read_text())["photos"]:
            i = Path(photo["image"]).stem
            if prepared.get(i) != photo["rgb_sha256"]:
                raise ValueError(f"Missing or different frozen image {i}")
            refs = [
                (kind, xy)
                for kind, sites in photo["reference_sites"].items()
                for xy in sites
            ]
            photos.append((i, "issue", refs))
            reviews.extend(
                dict(id=i, x=r["x"], y=r["y"], include="1")
                for r in photo["reviewed_near_misses"]
            )
        return photos, reviews
    labels = read_csv(manifest / "reference_junctions.csv")
    photos = []
    for photo in read_csv(manifest / "images.csv"):
        i = Path(photo["image"]).stem
        if prepared.get(i) != photo["rgb_sha256"]:
            raise ValueError(f"Missing or different frozen image {i}")
        refs = [
            (r["type"], (float(r["x"]), float(r["y"])))
            for r in labels
            if r["image"] == photo["image"]
        ]
        photos.append((i, "issue", refs))
    decisions = {
        r["id"]: r["status"] for r in read_csv(manifest / "near_miss_review.csv")
    }
    reviews = [
        dict(id=Path(r["image"]).stem, x=r["x"], y=r["y"], include="1")
        for r in read_csv(manifest / "near_miss_candidates.csv")
        if decisions[r["id"]] == "negative"
    ]
    return photos, reviews


def score(folder, photos, reviews, prefix=""):
    detections = read_csv(folder / "junctions.csv")
    scores, sites, negatives = [], [], []
    for image_id, group, refs in photos:
        observed = [r for r in detections if r["id"] == image_id]
        for kind in ("T", "X"):
            points = [q for k, q in refs if k == kind]
            rows = [dict(id=image_id, kind=kind, x=q[0], y=q[1]) for q in points]
            for method in METHODS:
                pred = [
                    (float(r["x"]), float(r["y"]))
                    for r in observed
                    if r["method"] == method and r["type"] == kind
                ]
                matched = {i for i, _ in match(points, pred)}
                for i, row in enumerate(rows):
                    row[method] = int(i in matched)
                scores.append(
                    dict(
                        id=image_id,
                        panel=group,
                        type=kind,
                        method=method,
                        reference=len(points),
                        detected=len(pred),
                        matched=len(matched),
                    )
                )
            sites.extend(rows)
        for review in reviews:
            if review["id"] != image_id or review["include"] != "1":
                continue
            row = dict(id=image_id, x=review["x"], y=review["y"])
            center = np.array([float(row["x"]), float(row["y"])])
            for method in METHODS:
                flags = [
                    r["type"]
                    for r in observed
                    if r["method"] == method
                    and r["type"] in ("T", "X")
                    and np.linalg.norm(center - [float(r["x"]), float(r["y"])]) <= 4
                ]
                row[method] = "+".join(sorted(set(flags))) or "-"
            negatives.append(row)
    if len(negatives) != sum(r["include"] == "1" for r in reviews):
        raise ValueError("A reviewed negative is missing from the photo panel")
    write_csv(
        folder / (prefix + "scores.csv"),
        scores,
        ["id", "panel", "type", "method", "reference", "detected", "matched"],
    )
    write_csv(
        folder / (prefix + "sites.csv"), sites, ["id", "kind", "x", "y", *METHODS]
    )
    write_csv(
        folder / (prefix + "near-miss-scores.csv"),
        negatives,
        ["id", "x", "y", *METHODS],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("review_csv", type=Path)
    parser.add_argument(
        "--manifest",
        "--frozen",
        type=Path,
        help="Published issue #3 benchmark JSON or original archive directory",
    )
    args = parser.parse_args()
    score(args.directory, supplement(args.directory), read_csv(args.review_csv))
    if args.manifest:
        score(args.directory, *original(args.directory, args.manifest), prefix="issue-")
