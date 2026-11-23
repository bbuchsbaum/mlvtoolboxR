"""Join the two masked reviews to their key and summarize by detector and type."""
import argparse
from collections import Counter
from pathlib import Path

from scipy.spatial.distance import cdist

from review_topology import X_REVIEW_CUTOFF
from score import match, positions, read_csv, write_csv


def audit(key, ratings, detections, output, bands=False):
    sites = read_csv(key)
    decisions = read_csv(ratings)
    labels = {r["token"]: r for r in decisions}
    if (
        len(decisions) != len(labels)
        or len(labels) != len(sites)
        or set(labels) != {r["token"] for r in sites}
    ):
        raise ValueError("Each masked site needs one row of decisions")
    rows = []
    for site in sites:
        row = dict(
            site, **{k: labels[site["token"]][k] for k in ("reviewer_a", "reviewer_b")}
        )
        if any(
            row[k] not in ("T", "X", "both", "other", "unclear")
            for k in ("reviewer_a", "reviewer_b")
        ):
            raise ValueError(f"Unknown visual label for {site['token']}")
        rows.append(row)
    events = read_csv(detections)
    if bands:
        counts = Counter(
            ("upper" if float(r["confidence"]) >= X_REVIEW_CUTOFF else "lower", "X")
            for r in events
            if r["type"] == "X" and r["method"] == "topology"
        )
    else:
        counts = Counter(
            ({"supported": "previous", "topology": "candidate"}[r["method"]], r["type"])
            for r in events
            if r["type"] in ("T", "X")
        )
    summary = []
    sources = (
        ("upper", "lower", "reference")
        if bands
        else ("previous", "candidate", "reference")
    )
    for source in sources:
        for kind in ("T", "X"):
            subset = [r for r in rows if r["source"] == source and r["type"] == kind]
            if not subset:
                continue

            def confirmed(r, reviewer):
                return r[reviewer] in (kind, "both")

            summary.append(
                dict(
                    source=source,
                    type=kind,
                    population=counts.get((source, kind), ""),
                    reviewed=len(subset),
                    reviewer_a=sum(confirmed(r, "reviewer_a") for r in subset),
                    reviewer_b=sum(confirmed(r, "reviewer_b") for r in subset),
                    both_confirm=sum(
                        confirmed(r, "reviewer_a") and confirmed(r, "reviewer_b")
                        for r in subset
                    ),
                    same_classification=sum(
                        r["reviewer_a"] == r["reviewer_b"] for r in subset
                    ),
                    unclear_a=sum(r["reviewer_a"] == "unclear" for r in subset),
                    unclear_b=sum(r["reviewer_b"] == "unclear" for r in subset),
                )
            )
    write_csv(output, summary, list(summary[0]))
    return summary


def score_retention(folder, gap_ledger, output):
    """Count reference matches and reviewed-gap flags retained at the secondary X cutoff."""
    events = [
        r
        for r in read_csv(folder / "topology-remainder-predictions.csv")
        if r["method"] == "topology" and r["type"] == "X"
    ]
    sites = [
        r for r in read_csv(folder / "topology-remainder-sites.csv") if r["kind"] == "X"
    ]
    gaps = [r for r in read_csv(gap_ledger) if r["include"] == "1"]
    photos = read_csv(folder / "panel.csv")
    rows = []
    for band, cutoff in [("all", 0), ("upper", X_REVIEW_CUTOFF)]:
        predicted = [r for r in events if float(r["confidence"]) >= cutoff]
        matches = flags = 0
        for photo in photos:
            image_id = photo["id"]
            xy = positions([r for r in predicted if r["id"] == image_id])
            ref = positions([r for r in sites if r["id"] == image_id])
            neg = positions([r for r in gaps if r["id"] == image_id])
            matches += len(match(ref, xy))
            flags += int((cdist(neg, xy) <= 4).any(axis=1).sum())
        rows.append(
            dict(
                band=band,
                cutoff=cutoff,
                detected=len(predicted),
                matched=matches,
                reference=len(sites),
                X_flags=flags,
            )
        )
    write_csv(output, rows, list(rows[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("key", "ratings", "detections", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument(
        "--bands", action="store_true", help="Summarize high and low X score bands"
    )
    parser.add_argument(
        "--retention",
        type=Path,
        help="Write a secondary X-cutoff comparison for the remainder panel",
    )
    parser.add_argument(
        "--gap-ledger",
        type=Path,
        help="Reviewed gap decisions for the retention comparison",
    )
    args = parser.parse_args()
    if args.retention and not (args.bands and args.gap_ledger):
        parser.error("--retention requires --bands and --gap-ledger")
    print(audit(args.key, args.ratings, args.detections, args.output, args.bands))
    if args.retention:
        score_retention(args.detections.parent, args.gap_ledger, args.retention)
