"""Score the frozen 60-photo development crossing-group experiment.

The four prediction arms are deliberately kept separate: FDoG and artist
drawings live in ``BASE``; the historical curve ranker, photo-only grouped
candidate, and reference-assisted diagnostic each expose the ordinary
``<id>.curve_rank.*`` R export in their own directory.  The diagnostic is
reported, but is never eligible to pass the gate.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from contour_incidence import incident_correspondence
from score import read_csv, write_csv
from score_curve_rank import (
    BOOTSTRAPS,
    SEED,
    aggregate as _aggregate,
    incidence_metrics,
    positional_metrics,
    raster_metrics,
    read_artist,
    read_drawing,
)


METHODS = ("fdog", "baseline", "grouped", "diagnostic")
FROZEN_ROSTER = (
    Path(__file__).parents[2] / "docs/evidence/wireframe/curve-rank/development.csv"
)
ARCHIVED_CONTROL = FROZEN_ROSTER.parent / "development/curve-rank-photo-scores.csv"


def aggregate(rows):
    """The established aggregate plus the new explicit strict-X total."""
    output = _aggregate(rows)
    output["X_four_arm_match"] = sum(int(row["X_four_arm_match"]) for row in rows)
    output["X_four_arm_recall"] = (
        output["X_four_arm_match"] / output["X_ref"] if output["X_ref"] else 0.0
    )
    return output


def _four_distinct_directions(value, *, observed):
    """Return four distinct ray angles, or ``None`` for an invalid arm set."""
    directions = np.asarray(value, dtype=float)
    if observed and directions.shape == (4,):
        if not np.all(np.isfinite(directions)):
            return None
        # R exports observed arms in mathematical degrees (y up).
        angles = np.deg2rad(directions % 360)
    elif not observed and directions.shape == (4, 2):
        if not np.all(np.isfinite(directions)):
            return None
        lengths = np.linalg.norm(directions, axis=1)
        if np.any(lengths <= 1e-12):
            return None
        angles = np.arctan2(directions[:, 1], directions[:, 0])
    else:
        return None
    ordered = np.sort(angles % (2 * np.pi))
    if np.min(np.diff(np.r_[ordered, ordered[0] + 2 * np.pi])) <= 1e-12:
        return None
    return angles


def four_arm_pairs(artist, observed, result=None):
    """Return maximum one-to-one strict X pairs as ``(artist, observed)``.

    Unlike legacy incidence, a strict crossing requires explicitly exported
    and mutually matched directions: exactly four artist rays, exactly four
    observed rays, and ``pair_evidence``'s direction check to be True.
    ``result`` permits callers that already computed incidence to reuse it.
    """
    if result is None:
        result = incident_correspondence(
            artist, observed["sites"], observed["contours"]
        )
    artist_x = [j for j in artist["junctions"] if j["type"] == "X"]
    observed_x = [s for s in observed["sites"] if s.get("type") == "X"]
    edges = [[] for _ in artist_x]
    for evidence in result["X"]["evidence"]:
        i, j = evidence["artist_index"], evidence["observed_index"]
        if not (evidence.get("accepted") and evidence.get("direction_agrees") is True):
            continue
        if (
            _four_distinct_directions(
                artist_x[i].get("incident_directions", []), observed=False
            )
            is None
            or _four_distinct_directions(observed_x[j].get("arms", []), observed=True)
            is None
        ):
            continue
        edges[i].append(j)
    assigned = {}

    def augment(i, seen):
        for j in edges[i]:
            if j in seen:
                continue
            seen.add(j)
            if j not in assigned or augment(assigned[j], seen):
                assigned[j] = i
                return True
        return False

    for i in range(len(edges)):
        augment(i, set())
    return sorted((i, j) for j, i in assigned.items())


def validate_roster(path):
    rows = read_csv(path)
    required = {"id", "category", "panel", "fold", "pixel_sha256"}
    if len(rows) != 60 or not rows or not required.issubset(rows[0]):
        raise ValueError("Crossing-group scorer requires the frozen 60-photo roster")
    if any(not r[k] for r in rows for k in required):
        raise ValueError("Frozen roster has an empty required value")
    if any(r["panel"] != "development" for r in rows):
        raise ValueError("Crossing-group scorer may score development photos only")
    if len({r["id"] for r in rows}) != 60:
        raise ValueError("Frozen roster IDs must be distinct")
    frozen = read_csv(FROZEN_ROSTER)
    if {(r["id"], r["pixel_sha256"]) for r in rows} != {
        (r["id"], r["pixel_sha256"]) for r in frozen
    }:
        raise ValueError("Roster differs from frozen development.csv")
    return sorted(rows, key=lambda r: (r["category"], int(r["fold"]), r["id"]))


def _geometry(staged, photo):
    ident = photo["id"]
    panel = {r["id"]: r for r in read_csv(staged / "panel.csv")}.get(ident)
    source = staged / f"{ident}.json"
    if (
        not source.is_file()
        or not panel
        or panel.get("category") != photo["category"]
        or panel.get("panel") != "development"
    ):
        raise ValueError(f"Roster/staged panel mismatch for {ident}")
    if panel.get("pixel_sha256") != photo["pixel_sha256"]:
        raise ValueError(f"Roster/staged pixel hash mismatch for {ident}")
    try:
        data = json.loads(source.read_text())
        if (
            data.get("id") != ident
            or data.get("category") not in (None, photo["category"])
            or data.get("panel") not in (None, "development")
        ):
            raise ValueError("source metadata differs")
        working = data["working"]
        return data, (int(working["height"]), int(working["width"]))
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid staged geometry for {ident}") from exc


def score(roster, staged, base, baseline, grouped, diagnostic):
    rows = []
    arms = {
        "fdog": (base, "fdog"),
        "baseline": (baseline, "curve_rank"),
        "grouped": (grouped, "curve_rank"),
        "diagnostic": (diagnostic, "curve_rank"),
    }
    for photo in roster:
        ident = photo["id"]
        geometry, shape = _geometry(staged, photo)
        artist, status = read_artist(base / f"{ident}.artist", shape)
        if status != "ok":
            raise ValueError(
                f"Artist drawing must be present and valid for {ident}: {status}"
            )
        for method, (folder, export) in arms.items():
            drawing, observed, status = read_drawing(
                folder / f"{ident}.{export}", shape
            )
            incidence = incident_correspondence(
                geometry, observed["sites"], observed["contours"]
            )
            row = dict(
                id=ident,
                category=photo["category"],
                panel="development",
                fold=photo["fold"],
                method=method,
                output_status=status,
            )
            row.update(raster_metrics(drawing, artist))
            row.update(positional_metrics(geometry, observed))
            row.update(incidence_metrics(geometry, observed))
            row["X_four_arm_match"] = len(four_arm_pairs(geometry, observed, incidence))
            rows.append(row)
    return rows


def summaries(rows):
    output = []
    for category in ["all", *sorted({r["category"] for r in rows})]:
        selected = (
            rows
            if category == "all"
            else [r for r in rows if r["category"] == category]
        )
        for method in METHODS:
            output.append(
                dict(
                    category=category,
                    method=method,
                    **aggregate([r for r in selected if r["method"] == method]),
                )
            )
    return output


def bootstrap(rows):
    by_method = {m: {r["id"]: r for r in rows if r["method"] == m} for m in METHODS}
    ids = sorted(by_method["grouped"])
    if any(set(table) != set(ids) for table in by_method.values()):
        raise ValueError("Methods are not aligned to exactly the roster photos")
    categories = {r["id"]: r["category"] for r in by_method["grouped"].values()}
    rng = np.random.default_rng(SEED)
    samples = np.concatenate(
        [
            rng.choice(
                [i for i in ids if categories[i] == c],
                size=(BOOTSTRAPS, sum(categories[i] == c for i in ids)),
                replace=True,
            )
            for c in sorted(set(categories.values()))
        ],
        axis=1,
    )
    metrics = (
        "pixel_coverage",
        "pixel_F1",
        "T_position_match",
        "X_position_match",
        "T_incidence_match",
        "X_incidence_match",
        "X_four_arm_match",
    )
    diffs = []
    for draw in samples:
        selected = {
            method: aggregate([by_method[method][ident] for ident in draw])
            for method in ("grouped", "baseline")
        }
        diffs.append(
            [selected["grouped"][key] - selected["baseline"][key] for key in metrics]
        )
    diffs = np.asarray(diffs)
    candidate, control = (
        aggregate(list(by_method[m].values())) for m in ("grouped", "baseline")
    )
    bounds = np.quantile(diffs, [0.025, 0.975], axis=0)
    return [
        dict(
            method="grouped",
            comparator="baseline",
            metric=key,
            estimate=candidate[key] - control[key],
            lo=lo,
            hi=hi,
            resamples=BOOTSTRAPS,
            seed=SEED,
        )
        for key, (lo, hi) in zip(metrics, bounds.T, strict=True)
    ]


def gates(rows):
    tables = {m: aggregate([r for r in rows if r["method"] == m]) for m in METHODS}
    grouped, fdog, baseline = tables["grouped"], tables["fdog"], tables["baseline"]
    grouped_rows = {r["id"]: r for r in rows if r["method"] == "grouped"}
    fdog_rows = {r["id"]: r for r in rows if r["method"] == "fdog"}
    baseline_rows = {r["id"]: r for r in rows if r["method"] == "baseline"}
    strict_gain_photos = sum(
        grouped_rows[i]["X_four_arm_match"] > baseline_rows[i]["X_four_arm_match"]
        for i in grouped_rows
    )
    checks = {
        "all_required_outputs_present": all(r["output_status"] == "ok" for r in rows),
        "every_grouped_photo_ink_at_most_fdog": all(
            grouped_rows[i]["ink"] <= fdog_rows[i]["ink"] for i in grouped_rows
        ),
        "pooled_pixel_coverage_gain_at_least_0_02": grouped["pixel_coverage"]
        >= fdog["pixel_coverage"] + 0.02,
        "pooled_pixel_F1_at_least_baseline_minus_0_005": grouped["pixel_F1"]
        >= baseline["pixel_F1"] - 0.005,
        "T_positional_matches_at_least_baseline": grouped["T_position_match"]
        >= baseline["T_position_match"],
        "X_positional_matches_at_least_baseline": grouped["X_position_match"]
        >= baseline["X_position_match"],
        "strict_X_four_arm_gain_at_least_2": grouped["X_four_arm_match"]
        >= baseline["X_four_arm_match"] + 2,
        "strict_X_four_arm_gain_on_at_least_2_photos": strict_gain_photos >= 2,
    }
    return {
        "gate_pass": all(checks.values()),
        "checks": checks,
        "values": tables,
        "strict_X_four_arm_gain_photos": strict_gain_photos,
        "thresholds": {
            "pixel_coverage_gain": 0.02,
            "pixel_F1_loss": 0.005,
            "strict_X_four_arm_gain": 2,
            "strict_X_four_arm_photos": 2,
            "bootstrap_seed": SEED,
        },
        "diagnostic_qualifies": False,
    }


def verify_control(rows):
    archived = {(r["id"], r["method"]): r for r in read_csv(ARCHIVED_CONTROL)}
    fields = (
        "ink",
        "artist_ink",
        "pixel_precision",
        "pixel_coverage",
        "pixel_F1",
        "T_position_pred",
        "T_position_match",
        "X_position_pred",
        "X_position_match",
        "T_ref",
        "T_incidence_match",
        "X_ref",
        "X_incidence_match",
    )
    for row in rows:
        if row["method"] != "baseline":
            continue
        expected = archived.get((row["id"], "curve_rank"))
        if not expected or any(
            not np.isclose(float(row[k]), float(expected[k]), rtol=0, atol=1e-12)
            for k in fields
        ):
            raise ValueError(
                f"Baseline control does not reproduce archived development score for {row['id']}"
            )


def main(args):
    roster = validate_roster(args.roster)
    rows = score(
        roster, args.staged, args.base, args.baseline, args.grouped, args.diagnostic
    )
    verify_control(rows)
    args.report.mkdir(parents=True, exist_ok=True)
    write_csv(args.report / "crossing-group-photo-scores.csv", rows, list(rows[0]))
    summary = summaries(rows)
    write_csv(args.report / "crossing-group-summary.csv", summary, list(summary[0]))
    intervals = bootstrap(rows)
    write_csv(
        args.report / "crossing-group-paired-bootstrap.csv",
        intervals,
        list(intervals[0]),
    )
    (args.report / "crossing-group-gate.json").write_text(
        json.dumps(gates(rows), indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("staged", type=Path)
    parser.add_argument("base", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("grouped", type=Path)
    parser.add_argument("diagnostic", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("--roster", type=Path, default=FROZEN_ROSTER)
    main(parser.parse_args())
