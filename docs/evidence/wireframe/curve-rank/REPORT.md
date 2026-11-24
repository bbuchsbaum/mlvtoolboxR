# Selecting whole contours improved drawings on 120 held-out Toronto photographs

At the same per-photo ink ceiling, the frozen selector covered **46.27% of artist contour pixels**, up from **43.77% with FDoG**. Contour F1 rose from **0.557 to 0.591**. The finished drawings matched **511 T junctions**, compared with **330** for the previous vector blend. The evidence for better crossings is weaker. The selector matched 37 X locations, but none had all the required incident geometry. [Held-out scores](reserve-test/curve-rank-summary.csv) · [Decision criteria](reserve-test/curve-rank-gate.json)

## Held-out comparison

| Method | Distinct ink pixels | Artist coverage | Contour F1 | T matches / detections | X matches / detections | Matches with incident geometry, T / X |
|---|---:|---:|---:|---:|---:|---:|
| FDoG | 269,683 | 43.77% | 0.557 | 76 / 279 | 0 / 23 | 1 / 0 |
| Previous vector blend | 269,666 | 43.05% | 0.559 | 330 / 867 | 28 / 390 | 28 / 3 |
| M-LSD alone | 216,797 | 38.15% | 0.530 | 1,008 / 2,595 | 86 / 1,054 | 80 / 4 |
| **Whole-contour selector** | **269,625** | **46.27%** | **0.591** | **511 / 1,030** | **37 / 517** | **32 / 0** |

The reference contains 704,761 artist ink pixels, 17,368 Ts and 2,604 Xs. The selector used 11,543 of 22,217 available curves. All 120 drawings stayed within their FDoG budgets. M-LSD alone found the most junctions and covered less of the artist drawing. [Per-photo counts](reserve-test/curve-rank-photo-scores.csv) · [Prediction hashes and budgets](reserve-test/selection.csv)

The coverage gain over FDoG was **2.49 percentage points** (paired 95% interval **2.00 to 2.98** points, with an unrounded lower bound of 1.997). The F1 gain was **0.0335** (0.0278 to 0.0391). Coverage increased in every category, from 0.28 points in mountains to 4.07 on highways. Compared with the previous blend, the difference was **+181 T matches** (103 to 261) and **+9 X matches** (−19 to 33). The X difference is compatible with a loss. Under the incidence test the differences were +4 Ts (−10 to 17) and −3 Xs (−9 to 0). [Paired photo intervals](reserve-test/curve-rank-paired-bootstrap.csv) · [Category results](reserve-test/curve-rank-summary.csv)

The 120 test IDs were fixed in `ef43e66`. The development fit and scoring rule were frozen in `61b5753`. All 120 predicted vector files were hashed and recorded in `0eb2007` before artist masks or geometry were staged for scoring. The other **205 reserve photographs have no method predictions or scores** from this trial. [Test roster](reserve-test.csv) · [Full source roster](../toronto/photos.csv)

## What the comparison measures

The selector can keep or discard a whole curve extracted from the photograph by FDoG or M-LSD. It sees local contrast, path shape, nearby curves, and how much new area each path would reach. Artist tracings supply its development targets. At inference the only inputs are the photograph and the two source drawings. Each photo has its own ceiling on ink, equal to the number of distinct pixels in its FDoG drawing. The R package renders the selected vectors and detects junctions on the completed drawing. [Frozen protocol](PLAN.md)

Coverage is the fraction of artist ink within two pixels of predicted ink. F1 also accounts for the fraction of predictions near artist ink. A positional match requires the same junction type within four pixels, with one-to-one matching. The stricter test also checks the paths entering the junction. These references describe junctions in artist tracings of salient contours. They cover a subset of visible photographic structure. [Artist reference and input checks](../toronto/REPORT.md) · [Incidence test](../../../../tools/wireframe/contour_incidence.py)

The 60 development photos were divided into five folds of 12, each with two photos per category. Every fold’s 12 drawings were predicted using only the other 48 photos for training. The development comparison passed all predeclared criteria. An independent rerun reproduced all 60 vector exports byte for byte. Erasing all artist pixels and junction annotations for one held-out fold produced the same 12 drawings. After refitting on all 60 photos, the numeric model reproduced the scikit-learn fit to within 1e-14 on every training example. Its SHA-256 is `d658e154644a0677509f0f460d099afe36c9b5ebc93e2344c36c47343492dade`. [Development comparison](development/curve-rank-gate.json) · [Fitted model](model.json)

Individual R rasterizations reproduced the stored pixel membership of all 9,843 development curves and all 22,217 test curves. Independent dense assignment and nearest-pixel calculations reproduced every reported positional match count and pixel coverage/precision numerator for all four methods on all 180 photographs. The development controls reproduced the earlier pilot. The model has been tested on held-out Toronto photos from these six categories. Accuracy on a different photograph-and-tracing corpus remains open. [Reproduction](REPRODUCE.md) · [Earlier pilot](../toronto/pilot-scores.csv)

## Where crossings are lost

After scoring, we checked the development drawings to see which paths the selector had discarded. M-LSD alone made 26 positional X matches. At only **2 of 26** did the selector keep all the M-LSD member curves. At 23 it kept some members and at one it kept none. All three development Xs where M-LSD met the strict incidence criterion lost at least one member. This check tracks exact source geometry. Other paths in the final drawing can create substitute junctions. [Development member counts](development/source-members.csv) · [Audit code](../../../../tools/wireframe/audit_curve_members.py)

The next experiment should let the selector take a pair or group of paths supporting a crossing together, paying only for the distinct ink they add. First test on development photos whether completing those groups can recover strict X matches while retaining the measured contour and T gains. Require incident-geometry improvement before spending another test panel. The present result supports the selector for contour coverage and T localization. Reliable recovery of four-arm crossings remains open.

The [crossing-group development experiment](../crossing-groups/REPORT.md) subsequently tested joint selection and the reference-assisted comparison.
