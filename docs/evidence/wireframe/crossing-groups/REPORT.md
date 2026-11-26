# Selecting crossing curves together recovered one four-arm X and lost 18 T matches

Keep the previous whole-contour selector. On the 60 Toronto development photos, joint selection increased correctly located Xs from 23 to 39. Only one of the 39 reproduced the artist’s four incident arms. T matches fell from 242 to 224. The frozen rule required two additional four-arm matches across at least two photos and preservation of the earlier T count. The 205 remaining reserve photos were left unscored. [Frozen protocol](PLAN.md) · [Decision and checks](crossing-group-gate.json)

## Finished drawings

| Method | Artist coverage | Contour F1 | T matches / detections | X matches / detections | Four-arm X matches |
|---|---:|---:|---:|---:|---:|
| FDoG | 42.65% | 0.5457 | 28 / 142 | 2 / 14 | 0 |
| Prior whole-contour selector | 44.91% | 0.5779 | 242 / 495 | 23 / 293 | 0 |
| Photo-only joint selection | 44.79% | 0.5770 | 224 / 472 | 39 / 415 | 1 |
| Reference-assisted diagnostic | 44.91% | 0.5779 | 242 / 496 | 26 / 294 | 3 |

All three selection policies used 134,701 distinct ink pixels in total and stayed within the FDoG budget on every photo. Against the previous selector, joint selection gained 16 positional X matches (paired 95% interval +8 to +24), lost 18 positional T matches (−31 to −5), and gained one four-arm X match (0 to +3). Contour coverage fell by 0.12 percentage points. It still exceeded FDoG by 2.14 points. The estimates resample whole photos within scene category. [Photo and category results](crossing-group-summary.csv) · [Paired intervals](crossing-group-paired-bootstrap.csv) · [Per-photo measurements](crossing-group-photo-scores.csv)

## What the diagnostic resolves

The rule gave a group of source curves extra weight when the photograph showed contrast on both sides of each of four distinct rays. Across the 60 photos, the two extractors supplied 396 X sites. Of these, 279 supplied a group with positive photographic support. The selector took a group action 159 times on 44 photos, and the completed drawings retained all members at 159 source X sites, up from 35. Many of those additional complete groups still disagreed with the artist. [Per-drawing group counts, budgets and vector hashes](selections.csv)

Three M-LSD sites on two photos already satisfied the four-arm criterion in their original source drawings. The earlier selector retained all the members of none of them. The photo-only rule retained one complete correct group. When the diagnostic used the artist labels to choose the three groups first, all three survived at their original artist references in the finished drawings, with the original 242 T matches. Thus the geometry and ink budget permitted all three recoveries. The tested contrast and angle rule recovered one. This diagnostic measures achievable recovery for these particular source groups. Different combinations of extracted paths can form other Xs. [Diagnostic selections](selections.csv) · [Four-arm matching implementation](../../../../tools/wireframe/score_crossing_groups.py)

The next development question is how often four suitably aligned arms are already present around an artist X even when neither extractor calls it an X. Measure that fraction separately from the fraction selected into the finished drawing. It distinguishes the need to extract better paths from the need to recognize promising combinations of existing paths. Three verified source groups on two photographs give little basis for learning a new completion bonus.

The [development branch diagnostic](../arm-opportunity/REPORT.md) measures this fraction on crossings with four traceable nine-pixel branches.

## Verification

The protocol was committed as `e31b2c8` before scoring. Both photo-only selectors used the same folds and path model. Each baseline prediction reproduced the earlier held-out vector file byte for byte. The R package rendered every final drawing and re-detected its junctions. The source labels used to choose diagnostic groups were read for that purpose only after all primary drawings had been predicted. [Reproduction steps](REPRODUCE.md)

A blinded review found that an X export can contain repeated directions. The scorer requires four finite, noncoincident rays on each side, a one-to-one directional match within 25 degrees, and agreement between the incident contours in a 3–9-pixel annulus. A horizontal bar with duplicate directions is a negative test. A second calculation enumerated all 24 arm permutations, checked both sets of contour samples, and used dense one-to-one assignment. It reproduced every four-arm count across all 240 drawings and checked the ink budgets against the rendered pixels. Two runs reproduced all 180 vector exports byte for byte. These are development results from the six Toronto scene categories. [Scoring tests](../../../../tools/wireframe/test_score_crossing_groups.py) · [Selection tests](../../../../tools/wireframe/test_crossing_groups.py)
