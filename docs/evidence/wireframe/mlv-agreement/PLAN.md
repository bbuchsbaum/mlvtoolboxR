# Agreement of MLV statistics with artist drawings

When the MLV contour statistics are computed on an automatic drawing of a photograph, how closely do they reproduce the statistics of the artist's tracing of the same photograph? Do they preserve the scene-category information carried by the artist statistics? This analysis measures both for the existing Toronto drawings. It fits no drawing method and makes no new predictions. Its results set the primary outcomes for later method comparisons.

## Drawings

Use the 180 Toronto photographs that already have drawings from all four methods: the 60 development photos, with out-of-fold selector drawings, and the 120 reserve-test photos, with drawings from the frozen selector. The methods are FDoG, M-LSD, the earlier vector blend and the whole-contour selector. The reference is the artist's polylines from the staged MATLAB files. The 205 unscored reserve photos stay unopened. [Development roster](../curve-rank/development.csv) · [Test roster](../curve-rank/reserve-test.csv)

Most of the analysis takes place in each photo's original 800 × 600 frame. Working-frame coordinates (256 × 192) map to the original frame by inverting the staging transform in `tools/wireframe/toronto_scenes.py`. Each artist drawing is rebuilt from its polylines, and each method drawing from its exported vectors, with `mlv_drawing(..., zero_length = "drop")`. FDoG drawings are the traced vector contours.

## Measurements

Compute `contour_statistics(drawing)` with the package defaults on every drawing. These cover orientation, length, curvature, horizontal/vertical, junction type and junction angle, in eight bins (five types for junction type). Curvature uses the turning estimator and junctions use the exit rule. Artist and method drawings of a photo use the same frame, so the default length range is identical for both. Per-segment lengths and curvatures come from `compute_contour_properties()`.

Two sensitivity analyses repeat every outcome with one change each:

- S1 uses the working 256-pixel frame for both drawings.
- S2 uses the published MLV definitions in the original frame: `feature_spec(estimator = "mlv2023", junction = list(rule = "mlv2023"), normalize = FALSE, length_range = c(2, width + height))`, covering orientation, length, curvature, junction type and junction angle.

## Outcomes

**O1, distribution agreement.** For each photo, method and property, normalize the method and artist histograms to proportions. Measure their Jensen–Shannon distance with base-2 logarithms, which lies between 0 and 1. If one histogram is empty and the other is not, the distance is 1; if both are empty, the pair is excluded and counted. Report the median and mean over photos. As a scale reference, report the median distance between artist drawings of different photos in the same category.

**O2, photo identification.** For each method, property and photo, rank the distance from the method histogram to its own artist histogram among its distances to the artist histograms of all 180 photos. Report the mean percentile rank, where 0.5 is chance and 1 means always closest. Report the top-1 rate as well.

**O3, scalar summaries across photos.** The summaries are total contour length, the number of contours, the density of each junction type per 10,000 pixels of contour length, the length-weighted mean absolute curvature, and the length-weighted fraction of near-horizontal and near-vertical orientation (within 22.5 degrees). For each, report the Spearman correlation across photos between method and artist, and the median method/artist ratio.

**O4, category information.** Describe each drawing with 38 features:

- the proportions of the orientation, length, curvature and junction-angle histograms (32);
- log(1 + density) for each of the five junction types (5);
- log total contour length (1).

Fit multinomial logistic regression (scikit-learn `LogisticRegression`, lbfgs, `C = 1`, `max_iter = 10000`), standardizing with means and scales from the training rows only. Use leave-one-photo-out cross-validation over the 180 photos and report three accuracies:

- artist → artist, the ceiling;
- method → method, the category information present in the method drawings;
- artist → method, which trains on other photos' artist drawings and tests on the held-out photo's method drawing.

Artist → method is the primary measure of whether automatic drawings give MLV statistics on the artist's scale. Chance is 1/6.

## Uncertainty and reporting

Paired differences between the selector and FDoG, and between the selector and the earlier blend, use 2,000 bootstrap resamples of photos within category, seed 20261126. Report pooled results for all 180 photos and separately for the 120 test photos. Report per-category O1 and O3. No outcome gates a decision here. The primary targets for later comparisons will be artist → method accuracy and the mean O1 distance across the six properties. Junction-specific outcomes (O1 junction type, O3 junction densities) will be reported alongside them. Every statistic is computed by one detector, in one frame, on both drawings. The MATLAB junction labels used in earlier point-matching scores are not used here.
