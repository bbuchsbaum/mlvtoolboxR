"""Evaluate frozen photo detections and traced drawings against Toronto artists' vectors."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import json
from pathlib import Path

import numpy as np
from PIL import Image

from contour_incidence import incident_correspondence
from contour_metrics import (
    covered,
    junction_counts,
    point_scores,
    reference_sites,
)
from score import positions, read_csv, write_csv


def measure(photo, staged, folder, pilot=False):
    i = photo["id"]
    geom = json.loads((staged / (i + ".json")).read_text())
    detected = json.loads((folder / (i + ".json")).read_text())
    artist = np.asarray(Image.open(folder / (i + ".artist.drawing.png"))) > 0
    incidence_cache = {}
    methods = (
        ["union", "union_vector", "fdog_retrace", "mlsd", "fdog"]
        if pilot
        else [*detected, "mlsd", "fdog"]
    )
    for method in methods:
        row = dict(
            id=i, panel=photo["panel"], category=photo["category"], method=method
        )
        # Localizations can be queried separately from the junctions formed by ink.
        point = point_scores(detected.get(method, []), geom, detected["context"])
        row.update(
            {
                k: v if method in detected or k.endswith("_ref") else np.nan
                for k, v in point.items()
            }
        )
        path = folder / f"{i}.{method}"
        drawing = np.asarray(Image.open(str(path) + ".drawing.png")) > 0
        sites = read_csv(Path(str(path) + ".sites.csv"))
        raw = Path(str(path) + ".contours.json").read_bytes()
        if raw not in incidence_cache:
            observed = json.loads(raw)
            match = incident_correspondence(
                geom, observed["sites"], observed["contours"]
            )
            incidence_cache[raw] = {k: match[k]["matched"] for k in ("T", "X")}
        for kind, hits in incidence_cache[raw].items():
            row[kind + "_incidence_hit"] = hits
        row.update(
            ink=int(drawing.sum()),
            artist_ink=int(artist.sum()),
            covered_pred=covered(drawing, artist),
            covered_artist=covered(artist, drawing),
        )
        for kind, narm in [("T", 3), ("X", 4)]:
            xy, angles = reference_sites(geom, kind)
            pts = [s for s in sites if s["type"] == kind]
            pred = positions(pts)
            observed = [
                np.deg2rad(np.fromstring(s["directions"], sep=" ")) for s in pts
            ]
            row[kind + "_draw_pred"] = len(pred)
            row[kind + "_draw_hit"] = junction_counts(xy, pred)
            eligible = np.array([len(a) == narm for a in angles], dtype=bool)
            row[kind + "_direction_ref"] = int(eligible.sum())
            row[kind + "_direction_hit"] = junction_counts(
                xy[eligible], pred, [a for a, k in zip(angles, eligible) if k], observed
            )
        yield row


def metrics(table, columns):
    """Return pooled micro-averages, using photographs as the resampling units."""
    total = {k: table[:, v].sum(axis=0) for k, v in columns.items()}

    def rate(a, b):
        return np.divide(a, b, where=b != 0, out=np.full(np.shape(a), np.nan))

    out = {}
    for k in ("T", "X"):
        for variant, num, den in (
            ("recall", "_hit", "_ref"),
            ("precision", "_hit", "_pred"),
            ("capped_recall", "_at_context_count", "_ref"),
            ("drawing_recall", "_draw_hit", "_ref"),
            ("incidence_recall", "_incidence_hit", "_ref"),
            ("drawing_precision", "_draw_hit", "_draw_pred"),
            ("direction_recall", "_direction_hit", "_direction_ref"),
        ):
            out[f"{k}_{variant}"] = rate(total[k + num], total[k + den])
    prec = rate(total["covered_pred"], total["ink"])
    rec = rate(total["covered_artist"], total["artist_ink"])
    out.update(
        ink_precision=prec,
        ink_recall=rec,
        ink_F1=rate(2 * prec * rec, prec + rec),
        mean_ink=total["ink"] / len(table),
    )
    return out


def summarize(rows, report, pilot=False):
    baseline = "fdog_retrace" if pilot else "context"
    fields = list(rows[0])[4:]
    columns = {c: j for j, c in enumerate(fields)}
    summaries, paired = [], []
    rng = np.random.default_rng(20261124)
    for split in sorted({r["panel"] for r in rows}):
        panel = [r for r in rows if r["panel"] == split]
        methods = sorted({r["method"] for r in panel})
        tables = {
            m: sorted(
                (r for r in panel if r["method"] == m),
                key=lambda r: (r["category"], r["id"]),
            )
            for m in methods
        }
        groups = np.array([r["category"] for r in tables[baseline]])
        ids = [r["id"] for r in tables[baseline]]
        samples = np.concatenate(
            [
                rng.choice(ix, (2000, len(ix)), replace=True)
                for category in sorted(set(groups))
                for ix in [np.flatnonzero(groups == category)]
            ],
            axis=1,
        )
        stats = {}
        for method, t in tables.items():
            if ids != [r["id"] for r in t]:
                raise ValueError("Unpaired photo comparison")
            arr = np.array([[r[k] for k in fields] for r in t], dtype=float)
            for category in ["all", *sorted(set(groups))]:
                sub = arr if category == "all" else arr[groups == category]
                summaries.append(
                    dict(
                        panel=split,
                        category=category,
                        method=method,
                        photos=len(sub),
                        **metrics(sub, columns),
                    )
                )
            # Chunk bootstrap samples to keep the largest temporary array small.
            chunks = [
                metrics(arr[block].transpose(1, 2, 0), columns)
                for block in np.array_split(samples, 10)
            ]
            stats[method] = {
                k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]
            }
        for method in (["union","union_vector"] if pilot else ["stem", "tracks", "both"]):
            for k in stats[method]:
                comparator = "fdog" if method == "union_vector" else baseline
                deltas = stats[method][k] - stats[comparator][k]
                if not np.isfinite(deltas).all():
                    continue
                lo, hi = np.quantile(deltas, [0.025, 0.975])
                paired.append(dict(panel=split, method=method, comparator=comparator, metric=k, lo=lo, hi=hi))
    for filename, items in (
        ("summary.csv", summaries),
        ("paired-intervals.csv", paired),
    ):
        write_csv(
            report / (("pilot-" if pilot else "") + filename), items, list(items[0])
        )


def score_photo(photo, staged, folder, pilot):
    return list(measure(photo, staged, folder, pilot))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("staged", "outputs", "report"):
        p.add_argument(name, type=Path)
    p.add_argument(
        "--panel", choices=["development", "confirmation", "all"], default="all"
    )
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--development-pilot", action="store_true")
    a = p.parse_args()
    if a.development_pilot:
        if a.panel == "confirmation":
            p.error("The coverage pilot is restricted to development photos")
        a.panel = "development"
    photos = read_csv(a.staged / "panel.csv")
    selected = [photo for photo in photos if a.panel in ("all", photo["panel"])]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = [
            r
            for group in pool.map(
                score_photo,
                selected,
                repeat(a.staged),
                repeat(a.outputs),
                repeat(a.development_pilot),
            )
            for r in group
        ]
    name = "pilot-scores.csv" if a.development_pilot else "toronto-scores.csv"
    write_csv(a.report / name, rows, list(rows[0]))
    summarize(rows, a.report, a.development_pilot)
    if a.development_pilot:
        raise SystemExit(0)
    audits = [
        dict(
            id=p["id"],
            panel=p["panel"],
            category=p["category"],
            **json.loads((a.outputs / (p["id"] + ".audit.json")).read_text()),
        )
        for p in photos
        if a.panel in ("all", p["panel"])
    ]
    write_csv(a.report / "candidate-evidence.csv", audits, list(audits[0]))
