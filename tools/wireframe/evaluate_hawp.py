"""Choose a HAWPv2 configuration on disjoint development photos, then evaluate once."""
import argparse
import subprocess
from pathlib import Path

from score import ROOT, original, read_csv, score, supplement, totals, write_csv

CUTOFFS = (0.05, 0.50, 0.90, 0.95)
VIEWS = ("vector", "mask")


def compare(folder, split, photos, reviews, cutoffs):
    plan = folder / f"hawp-{split}-plan.csv"
    events = folder / f"hawp-{split}-junctions.csv"
    jobs = [dict(id=i, cutoff=c) for i, _, _ in photos for c in cutoffs]
    write_csv(plan, jobs, ["id", "cutoff"])
    subprocess.run(
        [
            "Rscript",
            "--vanilla",
            "-e",
            'pkgload::load_all(".", quiet=TRUE); source("tools/wireframe/trace_hawp.R")',
            str(folder),
            str(plan),
            str(events),
        ],
        cwd=ROOT,
        check=True,
    )
    methods = [f"hawp_{view}_{c:.2f}" for c in cutoffs for view in VIEWS]
    scores, negatives = score(
        folder,
        photos,
        reviews,
        prefix=f"hawp-{split}-",
        detections=read_csv(events),
        methods=methods,
    )
    summary = []
    for cutoff in cutoffs:
        for view in VIEWS:
            method = f"hawp_{view}_{cutoff:.2f}"
            row = totals(split, method, scores, negatives)
            row["cutoff"] = cutoff
            summary.append(row)
    return summary


def selection_key(row):
    recall = [row[f"{k}_matched"] / row[f"{k}_reference"] for k in ("T", "X")]
    return (
        -row["flagged"],
        min(recall),
        sum(recall),
        -row["T_detected"] - row["X_detected"],
    )


def evaluate(folder, manifest, reviews):
    issue, issue_negatives = original(folder, manifest)
    issue_ids = {i for i, _, _ in issue}
    dev = [p for p in supplement(folder) if p[1] == "sample" and p[0] not in issue_ids]
    dev_ids = {p[0] for p in dev}
    dev_negatives = [r for r in read_csv(reviews) if r["id"] in dev_ids]
    development = compare(folder, "development", dev, dev_negatives, CUTOFFS)
    chosen = max(development, key=selection_key)
    print("Selected on disjoint development photos:", chosen, flush=True)
    frozen = compare(folder, "issue", issue, issue_negatives, [chosen["cutoff"]])
    for row in development + frozen:
        row["selected"] = int(row["method"] == chosen["method"])
    write_csv(folder / "hawp-summary.csv", development + frozen, list(chosen))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    evaluate(
        args.directory.resolve(),
        args.manifest,
        ROOT / "docs/evidence/wireframe/near-miss-review.csv",
    )
