"""Score a frozen whole-curve ranking trial against staged artist drawings.

The inputs deliberately separate the public R drawing stage from this scorer:
``EXPERIMENT_OUTPUTS`` contains only ``<id>.curve_rank`` drawings, while the
artist reference and the three fixed controls live elsewhere.  This command
never selects, fits, or ranks curves; it is an outcome-only, reproducible
measurement step.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from PIL import Image
from contour_metrics import covered, junction_counts

from contour_incidence import incident_correspondence
from score import read_csv, write_csv

METHODS = ("fdog", "union_vector", "mlsd", "curve_rank")
BOOTSTRAPS = 2000
SEED = 20261124


def blank(shape):
    return np.zeros(shape, dtype=bool)


def read_drawing(prefix, shape):
    """Read a R export, making a missing or empty export an explicit zero row."""
    image_path = Path(str(prefix) + ".drawing.png")
    sites_path = Path(str(prefix) + ".sites.csv")
    contours_path = Path(str(prefix) + ".contours.json")
    present = [p.is_file() for p in (image_path, sites_path, contours_path)]
    if not any(present):
        return blank(shape), {"sites": [], "contours": []}, "missing"
    if not all(present):
        return blank(shape), {"sites": [], "contours": []}, "incomplete"
    try:
        drawing = np.asarray(Image.open(image_path).convert("L")) > 0
        if drawing.shape != shape:
            return blank(shape), {"sites": [], "contours": []}, "wrong_shape"
        # Parse both exports: sites.csv is the declared public tabular drawing
        # export, contours.json supplies the member geometry needed for the
        # incidence check.  An empty sites CSV is valid.
        with sites_path.open(newline="") as f:
            list(csv.DictReader(f))
        observed = json.loads(contours_path.read_text())
        if not isinstance(observed, dict):
            raise ValueError("contours JSON is not an object")
        if not isinstance(observed.get("sites"), list) or not isinstance(
            observed.get("contours"), list
        ):
            raise ValueError("contours JSON has invalid sites or contours")
        return drawing, observed, "ok"
    except (OSError, ValueError, json.JSONDecodeError):
        return blank(shape), {"sites": [], "contours": []}, "invalid"


def read_artist(prefix, shape):
    """Artist references intentionally have only their raster drawing export."""
    path = Path(str(prefix) + ".drawing.png")
    if not path.is_file():
        return blank(shape), "missing"
    try:
        drawing = np.asarray(Image.open(path).convert("L")) > 0
        return (
            (drawing, "ok") if drawing.shape == shape else (blank(shape), "wrong_shape")
        )
    except OSError:
        return blank(shape), "invalid"


def raster_metrics(drawing, artist):
    ink, artist_ink = int(drawing.sum()), int(artist.sum())
    pred_near_artist = covered(drawing, artist)
    artist_near_pred = covered(artist, drawing)
    precision = pred_near_artist / ink if ink else 0.0
    coverage = artist_near_pred / artist_ink if artist_ink else 0.0
    f1 = (
        2 * precision * coverage / (precision + coverage)
        if precision + coverage
        else 0.0
    )
    return dict(
        ink=ink,
        artist_ink=artist_ink,
        pixel_precision=precision,
        pixel_coverage=coverage,
        pixel_F1=f1,
    )


def incidence_metrics(geometry, observed):
    result = incident_correspondence(geometry, observed["sites"], observed["contours"])
    out = {}
    for kind in ("T", "X"):
        item = result[kind]
        out[f"{kind}_ref"] = item["artist"]
        out[f"{kind}_incidence_match"] = item["matched"]
    return out


def positional_metrics(geometry, observed):
    """Same-type, 4-pixel, one-to-one junction position matches for the gate."""
    out = {}
    for kind in ("T", "X"):
        truth = np.asarray(
            [j["xy"] for j in geometry["junctions"] if j["type"] == kind], float
        ).reshape(-1, 2)
        predicted = np.asarray(
            [[s["x"], s["y"]] for s in observed["sites"] if s.get("type") == kind],
            float,
        ).reshape(-1, 2)
        out[f"{kind}_position_pred"] = len(predicted)
        out[f"{kind}_position_match"] = junction_counts(truth, predicted)
    return out


def validate_roster(path):
    rows = read_csv(path)
    required = {"id", "category", "panel", "fold"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Roster must have id, category, panel, and fold columns")
    ids = [r["id"] for r in rows]
    if any(not r[k] for r in rows for k in required - {"fold"}):
        raise ValueError("Roster has an empty required value")
    if any(not r["fold"] and r["panel"] != "reserve-test" for r in rows):
        raise ValueError("Only reserve-test roster rows may have an empty fold")
    if len(ids) != len(set(ids)):
        raise ValueError("Roster IDs must be distinct")
    return sorted(rows, key=lambda r: (r["category"], r["panel"], r["fold"], r["id"]))


def score(roster, staged, base, experiment):
    staged_panel = {r["id"]: r for r in read_csv(staged / "panel.csv")}
    rows = []
    for photo in roster:
        ident = photo["id"]
        source = staged / f"{ident}.json"
        if not source.is_file():
            raise ValueError(f"Roster/source mismatch: absent {source}")
        staged_row = staged_panel.get(ident)
        source_panel = "reserve" if photo["panel"] == "reserve-test" else photo["panel"]
        if (
            not staged_row
            or staged_row.get("category") != photo["category"]
            or staged_row.get("panel") != source_panel
        ):
            raise ValueError(f"Roster/staged panel mismatch for {ident}")
        if (
            photo.get("pixel_sha256")
            and staged_row.get("pixel_sha256") != photo["pixel_sha256"]
        ):
            raise ValueError(f"Roster/staged pixel hash mismatch for {ident}")
        try:
            geometry = json.loads(source.read_text())
            if geometry.get("id") != ident:
                raise ValueError("source id differs")
            expected = {"category": photo["category"], "panel": source_panel}
            for field, value in expected.items():
                if field in geometry and geometry[field] != value:
                    raise ValueError(f"source {field} differs from roster")
            working = geometry["working"]
            shape = (int(working["height"]), int(working["width"]))
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError(f"Invalid staged geometry for {ident}") from exc
        artist, status = read_artist(base / f"{ident}.artist", shape)
        if status != "ok":
            raise ValueError(
                f"Artist drawing must be present and valid for {ident}: {status}"
            )
        for method in METHODS:
            folder = experiment if method == "curve_rank" else base
            drawing, observed, status = read_drawing(
                folder / f"{ident}.{method}", shape
            )
            row = dict(
                id=ident,
                category=photo["category"],
                panel=photo["panel"],
                fold=photo["fold"],
                method=method,
                output_status=status,
            )
            row.update(raster_metrics(drawing, artist))
            row.update(positional_metrics(geometry, observed))
            row.update(incidence_metrics(geometry, observed))
            rows.append(row)
    return rows


def aggregate(rows):
    total = {
        k: sum(float(r[k]) for r in rows)
        for k in (
            "ink",
            "artist_ink",
            "T_ref",
            "X_ref",
            "T_position_pred",
            "X_position_pred",
            "T_position_match",
            "X_position_match",
            "T_incidence_match",
            "X_incidence_match",
        )
    }
    precision = (
        sum(r["pixel_precision"] * r["ink"] for r in rows) / total["ink"]
        if total["ink"]
        else 0.0
    )
    coverage = (
        sum(r["pixel_coverage"] * r["artist_ink"] for r in rows) / total["artist_ink"]
        if total["artist_ink"]
        else 0.0
    )
    f1 = (
        2 * precision * coverage / (precision + coverage)
        if precision + coverage
        else 0.0
    )
    out = dict(
        photos=len(rows),
        ink=int(total["ink"]),
        artist_ink=int(total["artist_ink"]),
        pixel_precision=precision,
        pixel_coverage=coverage,
        pixel_F1=f1,
    )
    for kind in ("T", "X"):
        out[f"{kind}_ref"] = int(total[f"{kind}_ref"])
        out[f"{kind}_position_pred"] = int(total[f"{kind}_position_pred"])
        for variant in ("position", "incidence"):
            out[f"{kind}_{variant}_match"] = int(total[f"{kind}_{variant}_match"])
            out[f"{kind}_{variant}_recall"] = (
                total[f"{kind}_{variant}_match"] / total[f"{kind}_ref"]
                if total[f"{kind}_ref"]
                else 0.0
            )
    out["non_ok_outputs"] = sum(r["output_status"] != "ok" for r in rows)
    return out


def summaries(rows):
    result = []
    for category in ["all", *sorted({r["category"] for r in rows})]:
        group = (
            rows
            if category == "all"
            else [r for r in rows if r["category"] == category]
        )
        for method in METHODS:
            item = aggregate([r for r in group if r["method"] == method])
            result.append(dict(category=category, method=method, **item))
    return result


def bootstrap(rows):
    by_method = {m: {r["id"]: r for r in rows if r["method"] == m} for m in METHODS}
    ids = sorted(by_method["fdog"])
    if any(set(table) != set(ids) for table in by_method.values()):
        raise ValueError("Methods are not aligned to exactly the roster photos")
    category = {r["id"]: r["category"] for r in by_method["fdog"].values()}
    rng = np.random.default_rng(SEED)
    metrics = (
        "pixel_coverage",
        "pixel_F1",
        "T_position_match",
        "X_position_match",
        "T_incidence_match",
        "X_incidence_match",
    )
    sampled = []
    for value in sorted(set(category.values())):
        group = [ident for ident in ids if category[ident] == value]
        sampled.append(rng.choice(group, size=(BOOTSTRAPS, len(group)), replace=True))
    samples = np.concatenate(sampled, axis=1)
    comparisons = {"fdog": metrics[:2], "union_vector": metrics[2:]}
    result = []
    for comparator, measures in comparisons.items():
        differences = []
        for draw in samples:
            candidate, control = [
                aggregate([by_method[m][i] for i in draw])
                for m in ("curve_rank", comparator)
            ]
            differences.append([candidate[k] - control[k] for k in measures])
        candidate, control = [
            aggregate(list(by_method[m].values())) for m in ("curve_rank", comparator)
        ]
        bounds = np.quantile(differences, [0.025, 0.975], axis=0)
        for metric, (lo, hi) in zip(measures, bounds.T, strict=True):
            result.append(
                dict(
                    method="curve_rank",
                    comparator=comparator,
                    metric=metric,
                    estimate=candidate[metric] - control[metric],
                    lo=lo,
                    hi=hi,
                    resamples=BOOTSTRAPS,
                    seed=SEED,
                )
            )
    return result


def gates(rows):
    tables = {m: aggregate([r for r in rows if r["method"] == m]) for m in METHODS}
    candidate, fdog, union = (
        tables["curve_rank"],
        tables["fdog"],
        tables["union_vector"],
    )
    per_image = {r["id"]: r for r in rows if r["method"] == "curve_rank"}
    fdog_rows = {r["id"]: r for r in rows if r["method"] == "fdog"}
    checks = {
        "pooled_pixel_coverage_gain_at_least_0_02": candidate["pixel_coverage"]
        >= fdog["pixel_coverage"] + 0.02,
        "pooled_pixel_F1_at_least_fdog": candidate["pixel_F1"] >= fdog["pixel_F1"],
        "T_positional_matches_at_least_union_vector": candidate["T_position_match"]
        >= union["T_position_match"],
        "X_positional_matches_at_least_union_vector": candidate["X_position_match"]
        >= union["X_position_match"],
        "every_photo_ink_at_most_fdog": all(
            per_image[i]["ink"] <= fdog_rows[i]["ink"] for i in per_image
        ),
        "all_required_outputs_present": all(r["output_status"] == "ok" for r in rows),
    }
    return {
        "gate_pass": all(checks.values()),
        "checks": checks,
        "values": {
            "curve_rank": tables["curve_rank"],
            "fdog": fdog,
            "union_vector": union,
        },
        "thresholds": {
            "pixel_coverage_gain": 0.02,
            "pixel_tolerance": 2,
            "junction_tolerance": 4,
            "bootstrap_seed": SEED,
        },
    }


def main(args):
    roster = validate_roster(args.roster)
    args.report.mkdir(parents=True, exist_ok=True)
    rows = score(roster, args.staged, args.base_outputs, args.experiment_outputs)
    write_csv(args.report / "curve-rank-photo-scores.csv", rows, list(rows[0]))
    summary = summaries(rows)
    write_csv(args.report / "curve-rank-summary.csv", summary, list(summary[0]))
    intervals = bootstrap(rows)
    write_csv(
        args.report / "curve-rank-paired-bootstrap.csv", intervals, list(intervals[0])
    )
    (args.report / "curve-rank-gate.json").write_text(
        json.dumps(gates(rows), indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("staged", type=Path)
    parser.add_argument("base_outputs", type=Path)
    parser.add_argument("experiment_outputs", type=Path)
    parser.add_argument("roster", type=Path)
    parser.add_argument("report", type=Path)
    main(parser.parse_args())
