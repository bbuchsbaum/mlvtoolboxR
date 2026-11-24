# Do contour-aware decisions improve the drawing?

Keep the existing context detector. On the 90 Toronto confirmation photos, distinguishing nearby Xs by the line tracks involved recovered one additional annotated X. All 90 output drawings remained pixel-identical to the context drawings. Checking T stems gave up 17 annotated Ts on Toronto and left the same two connected gaps on York. The experiment also found a way to recover more junctions than FDoG for the same quantity of ink, with a cost in artist-contour coverage. That result comes from development photos. [Fixed comparison and follow-up design](PLAN.md) · [Reproduction](REPRODUCE.md)

## Artist drawings

The Toronto Scenes authors paired photographs with salient contour drawings traced by trained artists. Junction classes in the MATLAB files describe the geometry of those drawings. They do not provide an exhaustive human inventory of visible junctions in each photograph. Photo-level hash ranking within six categories assigned 60 images to development, 90 to confirmation, and 325 to reserve. All 475 source dimensions matched, and all 475 drawings followed stronger image gradients at the recorded positions than after a fixed shift. No exact working-RGB or 64-bit difference-hash matches within four bits were found across Toronto splits or against the locally available Wireframe and York panels. This check covers 320 Wireframe fitting images and the enumerated comparison panels. The complete pretraining corpora were unavailable for duplicate screening. [Authors’ archive](https://osf.io/9squn/) · [Roster](photos.csv) · [Input checks](input-checks.json)

Each method starts with the same context proposals. The stem rule follows a proposed stem toward a measured bar and can reject a local interruption with distal contrast support. The crossing rule can retain two nearby Xs when each has a different pair of identifiable nonparallel through-lines. For the common experimental drawing, extend the actual finite segments only when their projected intersection is close. Cut a measured gap only on the implicated stem. Run the R package's rasterizer, contour tracer and junction detector on each resulting drawing.

| 90 confirmation photos | Context | Stem check | Distinct X tracks | Both |
|---|---:|---:|---:|---:|
| Matched T proposals / 11,748 artist Ts | 627 | 610 | 627 | 610 |
| Matched X proposals / 1,748 artist Xs | 53 | 53 | 54 | 54 |
| X proposals emitted | 840 | 840 | 851 | 851 |
| Matched Xs at the context count per photo | 53 | 53 | 54 | 54 |
| Matched T / X junctions found in the final drawing | 687 / 40 | 686 / 40 | 687 / 40 | 686 / 40 |
| T / X matches with local incident-contour agreement | 28 / 2 | 27 / 2 | 28 / 2 | 27 / 2 |
| Contour-pixel F1 within two pixels | 0.5126 | 0.5123 | 0.5126 | 0.5123 |

The extra X occurred in a city scene. Baseline and candidate each recovered zero Xs in the 15 beaches, 15 forests and 15 mountains. Only 213 of the 1,077 above-cutoff X candidates had two identifiable crossing tracks. The 95% paired interval for the change in X recall at equal count was 0.00 to +0.18 percentage points. For the change in T recall from checking stems it was −0.21 to −0.09 points. Intervals use 2,000 photo-level bootstrap resamples within category. A separate comparison with raw M-LSD drawings shows that the shared short-join renderer increased ordinary drawing matches from 660 to 687 Ts and 38 to 40 Xs. Its strict incident-contour matches went from 30 to 28 Ts and stayed at two Xs. [Per-photo counts](toronto-scores.csv) · [By category](summary.csv) · [Candidate evidence](candidate-evidence.csv) · [Paired intervals](paired-intervals.csv)

Point matches require the same type and one-to-one assignment within four working pixels. The local incident-contour check compares the artist's named curves with the R junction's actual member curves. In the three-to-nine-pixel annulus around each site, each sampled curve point must lie within two pixels of its counterpart's curve set, in both directions. Available arm directions must pair within 25 degrees. This is a strict local test. It does not establish that entire contours correspond across the scene. Contour precision and recall measure the proportion of predicted and reference ink pixels within two pixels of the other drawing. The F1 score is their harmonic mean. [Point and ink metrics](../../../../tools/wireframe/contour_metrics.py) · [Incident geometry](../../../../tools/wireframe/contour_incidence.py)

## Reviewed York gaps

All 5,964 context predictions on the 51 previously evaluated York photos reproduced in order, with identical coordinates, classes and scores. The stem check reduced proximity flags from 10 to nine of 24 reviewed gaps, including at the common operating point of 667 matched Ts. Reaching those 667 matches required 3,231 context predictions or 3,360 stem-checked predictions of type T. At the full operating point the methods matched 683 and 667 Ts. For the two methods, a connected ink path ran from the designated stem to the target at the same two reviewed sites. The path test requires eight-connected ink inside narrow corridors around the annotated stem and target. Proximity and connectivity measure different errors. [York gap ledger](york-gaps.csv) · [Per-photo scores](york-photos.csv)

The crossing change matched 398 versus 373 York Xs at its own output count, and 377 versus 373 at equal per-photo X counts. York was a reused panel. The Toronto results quantify how little of that gain transferred.

## Combining straight and curved contours at a fixed ink budget

In the 60-photo development follow-up, each drawing received the ink budget of FDoG on the same photo. Retain complete M-LSD paths, then add FDoG paths ranked by how much new territory they cover. Preserve the selected paths as vectors or rasterize and retrace them. Neither option uses artist vectors to select a path. Separate FDoG controls measure the effect of retracing by itself.

| 60 development photos | Ink pixels | Artist ink covered | Contour precision | Contour F1 | Drawing T / X matches | Incident-contour T / X matches |
|---|---:|---:|---:|---:|---:|---:|
| FDoG | 134,708 | 42.65% | 75.75% | 0.5457 | 28 / 2 | 0 / 0 |
| Combined, kept as vectors | 134,680 | 41.34% | 77.63% | 0.5395 | 163 / 13 | 9 / 0 |
| Retraced FDoG | 129,439 | 39.48% | 73.13% | 0.5128 | 60 / 0 | 0 / 0 |
| Combined, retraced | 128,693 | 39.46% | 75.94% | 0.5194 | 207 / 9 | 6 / 0 |
| M-LSD | 100,588 | 33.74% | 83.16% | 0.4800 | 430 / 26 | 23 / 3 |

The vector combination recovered 135 more T and 11 more X junctions than FDoG using 28 fewer ink pixels overall. It retained most of FDoG's contour coverage but fell 1.31 percentage points short. F1 fell 0.0062 (paired 95% interval for the change −0.0174 to +0.0050). M-LSD still recovered the most junctions. These are development comparisons motivated by the first experiment. This pilot left 325 reserve photos unscored. The subsequent [whole-contour comparison](../curve-rank/REPORT.md) used 120 of them. [Pilot counts](pilot-scores.csv) · [Category results](pilot-summary.csv) · [Paired intervals](pilot-paired-intervals.csv)
