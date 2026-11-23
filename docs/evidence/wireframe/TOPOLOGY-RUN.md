# Running the experimental photo junction scorer

## One image

The scorer takes an ordinary photo and returns one row per detected T or X. The CSV gives one-based pixel centers both on the 256-pixel-long-side working image and on the original after EXIF orientation. The method can remove or retype a previous prediction. It produces junction locations; the additional Xs do not constitute a unified contour drawing.

Use the same R package checkout and Python environment as the [previous experiment](FOLLOWUP.md#run-on-one-photo), including the supervised HAWPv2 dependencies. Download the pinned M-LSD model, HAWPv2 source and HAWPv2 checkpoint as described in the [M-LSD setup](REPORT.md#method-and-reproduction) and [HAWP setup](FOLLOWUP.md#reproduce-the-benchmark). Use the current CPU patch on pristine HAWP source. In addition to CPU compatibility, it returns an empty line set when fewer than two graph vertices exist. The small fitted T/X scorer is included in this checkout. From the repository root, run

```sh
py="$work/venv/bin/python"
"$py" tools/wireframe/photo_junctions.py photo.jpg junctions.csv --method topology \
  --weights "$work/mlsd.tflite" --hawp-source "$work/hawp-src" \
  --hawp-weights "$work/hawpv2-edb9b23f.pth"
```

## Reproduce the training and comparison

Start with the preceding experiment's [reproduction commands](FOLLOWUP.md#reproduce-the-benchmark). They create the 23, 24 and 48 calibration photo results and the separate original-panel results. Use a fresh working directory for a different code revision. Intermediate files include inputs, line predictions, feature matrices and scores. Only training requires scikit-learn.

```sh
uv pip install --python "$py" scikit-learn==1.7.2
"$py" tools/wireframe/prepare_topology_train.py "$work/topology-train"
"$py" tools/wireframe/prepare.py "$work/prospective" --manifest "$work/issue-3.json" --extra-panel prospective

"$py" tools/wireframe/predict_tiles.py "$work/prospective"
Rscript --vanilla -e 'pkgload::load_all(".", quiet=TRUE); source("tools/wireframe/trace_tiles.R")' "$work/prospective"
"$py" tools/wireframe/fragment_crossings.py "$work/prospective"
Rscript --vanilla -e 'pkgload::load_all(".", quiet=TRUE); source("tools/wireframe/compare.R")' "$work/prospective"

for folder in "$work/topology-train" "$work" "$work/replication" "$work/confirmation" "$work/prospective"; do
  "$py" tools/wireframe/predict_hawp.py "$folder" "$work/hawp-src" "$work/hawpv2-edb9b23f.pth"
  "$py" tools/wireframe/predict_hawp.py "$folder" "$work/hawp-src" "$work/hawpv2-edb9b23f.pth" --tiles
done

for variant in photoarms multiscale; do
  "$py" tools/wireframe/fit_topology.py "$work/topology-train" --variant "$variant" --fit "$work/$variant-model.json"
  for folder in "$work" "$work/replication" "$work/confirmation" "$work/prospective"; do
    "$py" tools/wireframe/fit_topology.py "$folder" --variant "$variant"
  done
  "$py" tools/wireframe/evaluate_topology.py "$work" "$work/$variant-model.json" "$work/$variant-cutoffs.json" \
    --variant "$variant" --select --panel issue
done

# Read the variant actually frozen in this checkout.
variant=$("$py" -c 'import gzip,json; print(json.load(gzip.open("tools/wireframe/junction-scorer.json.gz"))["variant"])')
for panel in development issue replication confirmation prospective; do
  "$py" tools/wireframe/evaluate_topology.py "$work" "$work/$variant-model.json" "$work/$variant-cutoffs.json" \
    --variant "$variant" --panel "$panel"
  case "$panel" in development|issue) folder="$work" ;; *) folder="$work/$panel" ;; esac
  "$py" tools/wireframe/check_topology.py "$folder" "$panel"
done
```

The training script checks the training metadata digest before selecting the 320 images. The prospective evaluation checks the original public manifest hash and checks the cardinality, identities and resized RGB hashes of all 503 photos across the five input directories. The scorer export compares numerical tree predictions with scikit-learn on every 23rd training candidate. The independent matching check uses sparse bipartite matching. Uncertainty estimates resample whole photographs in pairs 10,000 times using seed 812711. [Training code](../../../tools/wireframe/fit_topology.py) · [Paired checks](../../../tools/wireframe/check_topology.py)

To reproduce the reported predictions with the shipped fit, prepare the photos and multiscale features as above and use these model files in the evaluation commands. This skips fitting and threshold selection. The scorer's SHA-256 is `69e647dfdb0133620b15072bc85c6b061e820fda09719c7201407aa44812f436`.

```sh
"$py" - "$work" <<'PY'
import gzip, json, sys
from pathlib import Path

work = Path(sys.argv[1])
with gzip.open("tools/wireframe/junction-scorer.json.gz", "rt") as stream:
    scorer = json.load(stream)
for suffix, key in [("model", "models"), ("cutoffs", "cutoffs")]:
    (work / f"multiscale-{suffix}.json").write_text(json.dumps(scorer[key]))
PY
```

## Reproduce the masked review

Gap candidates are sampled by hash from annotated pairs with a 4–10 pixel gap and a clear target interior. Screening took place on the RGB crops. The four-pixel competing-junction radius is measured around the point on the target. Raters of the prediction cards received the cropped images, opaque identifiers and the same classification instructions. The sampling key contains the coordinates and method labels.

```sh
folder="$work/prospective"
"$py" tools/wireframe/review_topology.py "$folder" "$work/blind-gaps" "$folder/gap-key.csv"
"$py" tools/wireframe/review_topology.py "$folder" "$work/blind-predictions" "$folder/prediction-key.csv" \
  --predictions "$folder/topology-prospective-predictions.csv"
"$py" tools/wireframe/audit_topology.py "$folder/prediction-key.csv" \
  docs/evidence/wireframe/prospective-ratings.csv "$folder/topology-prospective-predictions.csv" "$folder/visual-summary.csv"
```

For the remaining test split, use the same prospective preparation, M-LSD, HAWP, feature and scoring commands with `remainder` in place of `prospective`. The roster contains 279 IDs. Evaluation excludes the one byte-identical copy of a calibration photo before matching, gap scoring and masked sampling. The 278 fresh photos contribute 19 eligible gaps. [Gap decisions](remainder-gap-review.csv) · [Duplicate-image check](TOPOLOGY-PLAN.md)

Verify the pooled photo comparison with

```sh
"$py" tools/wireframe/check_topology.py "$work" pooled --combine prospective remainder
```

To render the matched and missed reference Xs on the original floor photo, run

```sh
"$py" tools/wireframe/witnesses.py "$work" --floor-method topology --floor-recovery "$work/topology-floor.png"
```

The topology runner also writes a `score` for each predicted type. At earlier tracer positions this combines a nearby model score with a small same-type contribution from the tracer. This ranking score has no calibration to human visual certainty. The secondary X cutoff, fixed before the 278-photo extension, can be applied in R with `subset(read.csv("junctions.csv"), type == "T" | score >= 0.186)`.

To reconstruct the additional blinded samples, run

```sh
"$py" tools/wireframe/review_topology.py "$work/prospective" "$work/blind-references" \
  "$work/prospective/reference-key.csv" --references
"$py" tools/wireframe/review_topology.py "$work/remainder" "$work/blind-score-bands" \
  "$work/remainder/score-band-key.csv" --score-bands
"$py" tools/wireframe/audit_topology.py "$work/remainder/score-band-key.csv" \
  docs/evidence/wireframe/remainder-score-ratings.csv "$work/remainder/topology-remainder-predictions.csv" \
  "$work/remainder/visual-score-bands.csv" --bands --retention "$work/remainder/score-retention.csv" \
  --gap-ledger docs/evidence/wireframe/remainder-gap-review.csv
```

The extended comparison also checks that the 462 distinct test IDs cover the complete pinned test metadata and that their resized pixel hashes are distinct from those of all 320 training photos.
