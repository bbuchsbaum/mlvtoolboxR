# Recovering photographic T and X junctions

Across **342 previously unused Wireframe photographs**, the frozen scorer recovered **56.3% of T junctions** and **52.2% of X junctions**. The preceding crop-based detector recovered 40.8% and 36.7% respectively. The matched increases were +15.5 points (95% interval 13.5 to 17.4) for T and +15.5 points (95% interval 11.2 to 19.9) for X. The intervals resample whole photographs in pairs. The 95% interval for absolute X recall is 45.3% to 59.5%. [Paired results](topology-paired.csv)

| Photo set | T matches / reference | X matches / reference | Total detections | Reviewed gaps flagged |
| --- | ---: | ---: | ---: | ---: |
| Original issue, 24 | 20 → 33 / 62 (32.3% → 53.2%) | 42 → 53 / 101 (41.6% → 52.5%) | 1,540 → 1,495 | 0 → 0 / 15 |
| First prospective test, 64 | 457 → 664 / 1091 (41.9% → 60.9%) | 44 → 62 / 150 (29.3% → 41.3%) | 3,452 → 3,232 | 0 → 1 / 30 |
| Extension, 278 | 2109 → 2876 / 5201 (40.5% → 55.3%) | 412 → 587 / 1093 (37.7% → 53.7%) | 16,423 → 15,042 | 4 → 4 / 19 |
| **All fresh photos, 342** | 2566 → 3540 / 6292 (40.8% → 56.3%) | 456 → 649 / 1243 (36.7% → 52.2%) | 19,875 → 18,274 | 4 → 5 / 49 |

Each cell compares the preceding detector with the frozen scorer. A match requires the correct type within four pixels and each reference and detection can match only once. Each gap counts once if either predicted type falls within four pixels of its projected intersection. On the 95 calibration photographs, T matches rose from 762 to 1,002 of 1,688 and X matches from 157 to 215 of 363. Detections fell from 5,635 to 5,079, with four of 42 gaps flagged by each method. [Panel totals](topology-summary.csv) · [Per-photo results](topology-by-photo.csv) · [Per-gap results](topology-gaps.csv)

The original issue photographs exceed the 50% target for each type with zero of 15 reviewed gaps flagged. On the floor photograph, X matches rose from three to six of 33 and T matches from zero to two of three. Both detectors recovered six of the seven annotated window crossings, including all four sites marked in the issue (0.42–1.61 pixels from the marked centers). The first prospective test fell short of 50% for X. [Original frozen manifest](https://github.com/bbuchsbaum/mlvtoolboxR/issues/3#issuecomment-7101349061) · [Photo scores](topology-by-photo.csv)

We extended the evaluation after the first 64 results and used the same weights and cutoffs throughout. The remaining 279 official test IDs contained one byte-identical copy of a calibration photo, with a different annotation. We detected the duplicate by hash before scoring and excluded its image and eligible gap from the extension. The remaining 278 photos and 19 gaps enter the table. The 64-photo and 278-photo rows are disjoint. The pooled row combines them. [Decision and timing](TOPOLOGY-PLAN.md) · [Extension roster](remainder-panel.csv) · [Screened gaps](remainder-gap-review.csv)

## What changed

The older method retained every full-image M-LSD junction, then added some crop detections and four-arm crossings. The new method can reject a spurious site, retype a junction, or add a site where the image and predicted lines support it. Its candidate positions include vertices from HAWPv2 and intersections between nonparallel finite segments from HAWPv2 and OpenCV's line segment detector. At each candidate, it measures evidence in 36 directions. A T needs two roughly opposite arms and a third arm crossing that axis. An X needs a second pair of opposite arms. Photographic features also measure the image gradient perpendicular to each candidate arm at distances from 1.5 to 17 pixels. [Feature extraction and site selection](../../../tools/wireframe/topology.py)

We fitted two gradient-boosted tree scorers to those measurements on 320 photos from the official Wireframe **training** split. A training candidate was positive for a type when it lay within 3.5 pixels of a reference of that type. To discourage false bridges, we weighted negatives eightfold when they were close to the projected intersection of a line ending 4–10 pixels short of another line. These automatically generated training negatives are distinct from the visually reviewed evaluation sites. One pair of thresholds applies at locations proposed by the previous tracer. A second pair applies to newly proposed positions. All methods use the same one-to-one same-type scoring radius of four pixels. [Training code](../../../tools/wireframe/fit_topology.py) · [Training roster](topology-train-panel.csv)

The first vector-only model lost too much recall when we tightened its threshold to control gap flags. We then used vector evidence to type locations found by the earlier tracer, added image-gradient measurements, and tested four contrast-enhanced overlapping views. We tuned on the 95 earlier supplementary photos, maximizing the lower of T and X recall within the previous detection-count and gap-flag limits. The original issue panel entered finalist selection. Its numbers measure engineering progress on this specific problem. The fresh 64 supply the prospective estimate. All photos are from Wireframe, which the upstream HAWP authors also used for validation. [Decision plan](TOPOLOGY-PLAN.md) · [Finalist comparisons](topology-finalists.csv) · [Upstream split configuration](https://github.com/cherubicXN/hawp/blob/92ae446875ef9296fc95de7313519aabd638ebb5/configs/hawpv2.yaml)

## Which predictions look like junctions?

Two independent Codex Astra sessions saw the same shuffled crops from the first 64 photos. They saw the location, raw context and a magnification. Method names, predicted types, confidence and reference labels were withheld. We sampled 24 predictions per type and method. Unclear ratings count as unconfirmed. The two sessions assigned the same visual class on 84/108 cards, including the 12 interspersed references. Both sessions used the same model family. [Masked review key](prospective-rating-key.csv) · [Individual ratings](prospective-ratings.csv)

| Sample | Reviewer A confirmed | Reviewer B confirmed | Both confirmed |
| --- | ---: | ---: | ---: |
| Previous T | 15/24 | 14/24 | 11/24 |
| New T | 20/24 | 20/24 | 19/24 |
| Previous X | 4/24 | 1/24 | 1/24 |
| New X | 5/24 | 3/24 | 3/24 |
| Reference T controls | 4/6 | 4/6 | 4/6 |
| Reference X controls | 2/6 | 1/6 | 1/6 |

The visual gain is clearest for Ts. The absolute X confirmation rate is low even on reference controls. Two photographs account for 48 of the 88 missed reference Xs in the first 64. One contains overlapping ceiling beams and the other window muntins. Multiple reference intersections can lie a few pixels apart on the same wide bar. Some visible boundaries stop against foreground trim while the corresponding annotated straight lines intersect. [Per-photo scores](topology-by-photo.csv) · [Ceiling photo](https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test/images/00130203.png) · [Window photo](https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test/images/00113162.png)

A third masked Astra session checked this explanation on random reference sites from those 64 photos. Of **36 missed reference Xs**, it judged **5 X, 23 T, 7 other and 1 unclear**. On separate controls it called **7/12 matched reference Xs** an X and **11/12 matched reference Ts** a T. Geometric recall and agreement about what is visible therefore answer different questions. The frozen labels stay intact in all reported recall numbers. [Reference-site audit](prospective-reference-audit.csv)

## Using the X score to prioritize visual review

The 0.186 cutoff was set from the 95 calibration photos before reading the extension results. On random predictions from the 278 fresh extension photos, reviewer A confirmed an X on **17/36** above the cutoff and **4/36** below it. Reviewer B confirmed **16/36** and **4/36**. They agreed on the exact visual class for 74/90 cards. On the nine X reference controls they saw 4 and 3 Xs, respectively. [Masked sample key](remainder-score-key.csv) · [Individual ratings](remainder-score-ratings.csv) · [Visual summary](remainder-visual-summary.csv)

Applying the cutoff reduces the X review list from 5,467 to 1,235 predictions (23%) while retaining **407 of the 587 matched X references** (69%). Two X proximity flags at the 19 reviewed gaps also drop out. The CLI exports this score for ranking. Its calibration to human judgments and its performance on other datasets remain unmeasured. The experimental runner returns typed point locations. Producing a contour drawing that connects the newly recovered junctions remains an open problem. [Score retention](remainder-score-retention.csv) · [Run on one photograph](TOPOLOGY-RUN.md)

## Verification

The training and test rosters account for all 462 official test IDs and 320 distinct training inputs. We checked the resized RGB hashes of all 782 entries. A second bipartite matcher reproduced the per-photo counts. Single-photo and batch line tracing agreed on all 69,083 typed events across the 279 prepared extension IDs. Single-photo inference reproduced the typed batch positions and their original-image coordinates for the floor witness and a prospective photo. Blank-image inference produced a header-only CSV. The numerical controls cover three-arm and four-arm inputs, a six-pixel gap and score export. The five Python controls and the full R package test suite passed. [Verifier](../../../tools/wireframe/check_topology.py) · [Numerical controls](../../../tools/wireframe/test_topology.py) · [Reproduction commands](TOPOLOGY-RUN.md)
