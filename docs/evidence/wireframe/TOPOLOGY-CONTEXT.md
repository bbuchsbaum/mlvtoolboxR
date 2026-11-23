# Using image context to recover photographic junctions

The experimental `--method context` recovered more annotated T junctions than `--method topology` on both Wireframe and York Urban. On Wireframe, the T advantage persisted when both methods could report the same number of detections per photo. On York, the T advantage depended on reporting more detections. The evidence for recognizing visible four-arm crossings remains weak. [Single-photo runner](../../../tools/wireframe/photo_junctions.py) · [Per-photo counts](context-by-photo.csv)

## Frozen comparison

The new scorer uses the same candidate locations, contour anchors and four-pixel suppression rule as the topology scorer. It adds 146 measurements from the pretrained HAWPv2 image features before the detector converts them to lines. It also summarizes line evidence within 12 and 48 pixels and across the image, adding 89 descriptors. Neither inference nor training uses photo identifiers as features. [Feature extraction](../../../tools/wireframe/context_features.py) · [HAWP sampling](../../../tools/wireframe/predict_hawp.py)

We fitted separate T and X scorers on the same 320 Wireframe training photos. Four variants used image features, neighborhood summaries, or both in boosted trees, and both in a neural network. Each was tested alone and as a 50/50 blend with the existing scores. We calibrated four thresholds for each on the same 95 earlier photos. Total detections and gap flags could not exceed the baseline, and each type could produce up to 105% of its baseline count. Selection maximized the smaller of the two recall gains, breaking ties on their mean and then on fewer detections. The unblended trees with both feature sets won. They recovered 1,059/1,688 Ts and 232/363 Xs, up from 1,002 and 215. Detections fell from 5,079 to 5,009 and gap flags from four to three. We froze this fit and its thresholds before scoring the panels below. [Eight calibration results and cutoffs](context-calibration.csv) · [Training roster](topology-train-panel.csv)

| Evaluation panel | T matches, topology → context | X matches, topology → context | All detections | Reviewed gaps flagged |
| --- | ---: | ---: | ---: | ---: |
| Original issue, 24 photos | 33 → 37 / 62 | 53 → 53 / 101 | 1,495 → 1,502 | 0 → 0 / 15 |
| Wireframe, 342 previously evaluated photos | 3,540 → 3,781 / 6,292 | 649 → 658 / 1,243 | 18,274 → 18,508 | 5 → 5 / 49 |
| York Urban, 51 external test photos | 422 → 477 / 1,180 | 271 → 287 / 665 | 4,676 → 5,313 | 4 → 5 / 26 |

The 95% paired photo-bootstrap intervals on the 342 photos were **+2.73 to +4.93 percentage points** for the T recall gain of 3.83 points and **−1.00 to +2.50** for the X gain of 0.72 points. On York they were **+2.45 to +6.71** for T (+4.66) and **−0.77 to +5.16** for X (+2.41). Both matchers enforce one-to-one same-type matches within four pixels. The bootstrap resampled whole photos 10,000 times with seed 812711. The original floor photo went from six to eight of its 33 Xs. The pooled Wireframe gap count masks one fewer flag in the first 64 photos and one extra in the other 278. [Per-photo scores](context-by-photo.csv) · [Original issue manifest](https://github.com/bbuchsbaum/mlvtoolboxR/issues/3#issuecomment-7101349061)

## What survives the same detection budget?

York produced 582 more T detections and 55 more X detections. To measure ranking at equal output size, we then kept the highest exported scores from each method up to the **smaller count for that photo and type**. Both methods were truncated where appropriate. This is a post-evaluation comparison of their rankings. The model, thresholds and public runner use the frozen setting above.

| Panel | Type | Detections per method | Topology matches | Context matches | Recall difference in points (95% paired interval) |
| --- | --- | ---: | ---: | ---: | ---: |
| Wireframe, 342 | T | 11,102 | 3,484 | 3,596 | +1.78 (+0.98 to +2.60) |
| Wireframe, 342 | X | 5,657 | 633 | 646 | +1.05 (−0.63 to +2.73) |
| York, 51 | T | 2,444 | 418 | 419 | +0.08 (−2.04 to +1.94) |
| York, 51 | X | 1,968 | 266 | 284 | +2.71 (0.00 to +5.19) |

At equal counts both methods flagged four of 49 Wireframe gaps and four of 26 York gaps. The extra York flag at the frozen operating point missed our criterion for replacing the previous scorer. The context scorer is offered as an option for recovering more annotated junctions. [Equal-count scoring code](../../../tools/wireframe/score_context.py)

## What the York photographs show

We used all 51 photos in the official test split listed by Denis, Elder and Estrada (2008), *Efficient Edge-Based Methods for Estimating Manhattan Frames in Urban Imagery*. These are predominantly architectural scenes. The authors annotated line segments aligned with the dominant scene directions. We derived the T/X references with the same finite-segment rule as on Wireframe. York's annotations are selective and an intersection of annotated lines can be partly occluded. Resized RGB hashes were distinct from all prior Wireframe panels and the training images. A masked Astra reviewer accepted 26 of 45 proposed gap controls before the comparison. [York Urban Line Segment Database](https://www.elderlab.yorku.ca/resources/york-urban-line-segment-database-information/) · [Roster and hashes](york-panel.csv) · [Gap decisions](york-gap-review.csv)

Two independent Codex Astra sessions classified the same 108 shuffled photo crops with method, type, score and reference status hidden. Among 24 predicted Ts from each method, both reviewers saw a T on **17/24** for topology and **19/24** for context. Among 24 predicted Xs they both saw an X on **5/24** and **4/24**, respectively. They jointly confirmed only **3/6** of the annotated X control sites as visibly four-armed. This small same-model review gives no evidence that the context scorer improves visual X precision. Human-rated precision and performance outside these architectural datasets remain open. [Masked sample and individual decisions](context-york-ratings.csv)

## Reproduce

Use the pinned Python/R environment and model files in the earlier run instructions. To run the new scorer on one image or refit its trees from the prepared training data, use the commands below. The four operating thresholds are in the selected calibration row. The included compressed fit has SHA-256 `69a4ef9fb5715524b84385b1a3aae4b74eaff1ac461ed2501da397c0b4a92d6e`. The scores use this fit's own scale. The topology scorer's secondary X cutoff was calibrated for its own outputs. [Dependency setup](TOPOLOGY-RUN.md) · [Refit code](../../../tools/wireframe/fit_context.py)

```sh
py="$work/venv/bin/python"
"$py" tools/wireframe/photo_junctions.py photo.jpg junctions.csv --method context \
  --weights "$work/mlsd.tflite" --hawp-source "$work/hawp-src" \
  --hawp-weights "$work/hawpv2-edb9b23f.pth"

"$py" tools/wireframe/fit_context.py "$work/topology-train" "$work/hawp-src" \
  "$work/hawpv2-edb9b23f.pth" "$work/context-model.json"
```

To reproduce the York comparison from the authors' archive and the frozen predictions, run the following commands. These regenerate the official photo roster and its references, then check the full and equal-count results with two matching algorithms. For direct inference, run either scorer on the corresponding prepared PNG. Its working coordinates match the archived detection table. [Archive](https://www.elderlab.yorku.ca/?sdm_process_download=1&download_id=8288) · [Frozen York predictions](context-york-events.csv.gz)

```sh
"$py" tools/wireframe/prepare_york.py "$work/YorkUrbanDB.zip" "$work/york51"
"$py" tools/wireframe/score_context.py "$work/york51" \
  docs/evidence/wireframe/context-york-events.csv.gz "$work/york-scores.csv"
```

Refitting from the training images reproduced all split indices, thresholds and leaf values of the frozen context model exactly. The single-photo CLI reproduced all 74 context predictions and all 81 topology predictions on York image P1040798, including scores. The new runner returned an empty table with the expected columns for a blank input. Scoring the archived York detections against a fresh extraction of the authors’ annotations reproduced every full-count and equal-count row. [Scoring ledger](context-by-photo.csv)
