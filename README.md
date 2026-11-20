# mlvtoolbox

`mlvtoolbox` is an R package for measuring contour geometry and constructing
reproducible visual stimuli from line drawings. It is an independently
maintained R port and extension of the Bernhardt-Walther Lab's
[Mid-Level Vision Toolbox](https://github.com/bwlabToronto/MLV_toolbox), with
[pyMLV](https://github.com/bwlabToronto/pyMLV) as an additional reference.

The aim is to make the [original MLV methods](https://doi.org/10.3389/fcomp.2023.1140723)
usable in R with explicit measurement choices, documented bug fixes, tested
native performance improvements, and tools for matching stimulus sets.
Corrected definitions are the defaults; the published behaviour remains
available through explicit `mlv2023` options, and graph/design workflows are
named extensions. Numerical equivalence is qualified
against specific fixtures, not promised for every upstream output.

[Purpose and differences](docs/ATTRIBUTION.md) ·
[Credits and license audit](vendor/LICENSING.md) ·
[Measurement contract](docs/FEATURE_SPEC.md)

## Install locally

A C++17 compiler and Rcpp are required. From this repository:

```r
install.packages("Rcpp")
install.packages(".", repos = NULL, type = "source")
```

Optional packages: xml2 for SVG, R.matlab for v5 MAT files, and png/jpeg for
images. Python, Octave and model weights are not runtime dependencies.

## Drawing to features

Contours are ordered matrices with columns x1, y1, x2, y2. Image coordinates
increase rightward/downward; integer raster pixel centers start at (1,1).
Subpixel coordinates are retained. Zero-length segments require explicit
cleaning with zero_length="drop".

```r
library(mlvtoolbox)
drawing <- mlv_drawing(
  list(rbind(c(5, 5, 35, 5), c(35, 5, 35, 25)),
       matrix(c(20, 1, 20, 30), nrow = 1)),
  width = 40, height = 35
)
drawing <- compute_contour_properties(drawing)
as.data.frame(drawing)             # segment measurements
drawing$junctions                  # positions, types, arms and angles
contour_statistics(drawing)        # long image-level histograms
feature_table(drawing)             # normalized contour-level predictors
plot_properties(drawing)
```

Histogram weighting and normalization are explicit: orientation/curvature
collect segment length; length collects whole-contour length; junctions count
events. A value within 1e-9 × max(1, |edge|) of a bin edge is binned as if on
the edge, so mathematically equal angles share a bin on every platform. Normalized values are per 10,000 contour pixels. Valid zero scores
remain distinct from missing scores.

## How this port differs

- **An R workflow:** drawing objects, feature tables and model compatibility
  metadata; Python and MATLAB serve as references.
- **Explicit methods and corrections:** subpixel geometry, validated
  degeneracies and missingness, vertex-turning curvature by default (each turn
  counted once), an arc-length estimator at a stated scale, and physical
  calibration. `estimator = "mlv2023"` retains the published curvature.
- **Additional stimulus tools:** constrained subset design and contour graph
  descriptors for closure, continuation and enclosure.
- **Measured implementation changes:** compiled geometry and fused histograms,
  with [performance and equality evidence](docs/FUSED_STATISTICS.md).
- **Two modes for curvature, junctions and the medial axis:** corrected defaults, and
  explicit `mlv2023` options that reproduce outputs the original lab computed
  in MATLAB (see [parity evidence](docs/PARITY.md)). No bundled SAM weights or
  MATLAB serialized models.

See the [comparison with the originals](docs/ATTRIBUTION.md) for scope, source
lineage and qualification boundaries.

## Medial-axis analysis

```r
mask <- rasterize(drawing)
result <- analyse_mat(drawing, min_branch_points = 2)
plot(result$drawing, property = "parallelism")
scores <- contour_statistics(result$drawing, properties = "parallelism")
rgb <- render_drawing(result$drawing, property = "parallelism")
```

compute_mat() returns distances, nearest-boundary indices, AOF, skeleton,
ordered branches, and method parameters. compute_mat_properties() adds
parallelism, separation, mirror and optional taper ratings.
map_mat_properties() transfers them back to contours.

By default the flux is pixel-centred and independent of how equidistant
boundary pixels are labelled, the skeleton keeps no removable pixel, and branch
scores do not depend on branch direction or numbering. The published-compatible
recipe reproduces the lab's stored MATLAB results: on 2,757 drawings, 99.1–99.4%
of per-pixel scores agree within 1e-4 (residual: MATLAB's single precision).

```r
compat <- compute_mat(drawing, aof = "mlv2023", morphology = "hmt8",
                      branch_rule = "mlv2023",
                      threshold_angle = asin(0.3 * pi / 2) * 180 / pi)
```

These comparisons qualify computation against stored outputs, not the
perceptual validity of the measures. See [parity evidence](docs/PARITY.md) and
[implementation status](IMPLEMENTATION_STATUS.md).

## Learn more

- `vignette("getting-started", package = "mlvtoolbox")`: from a line drawing to contour features.
- `vignette("junctions", package = "mlvtoolbox")`: how junctions are found and what they measure.
- `vignette("medial-axis", package = "mlvtoolbox")`: medial axis and symmetry scores.
- `vignette("mlv-compatibility", package = "mlvtoolbox")`: reproducing published MLV values.
- `vignette("workflows", package = "mlvtoolbox")`: file input, extraction and feature recipes.
- `vignette("stimulus-design", package = "mlvtoolbox")`: matched conditions across a stimulus set.

## Manipulate and import

```r
parts <- split_contours(drawing, property = "curvature", fraction = 0.5)
rotated <- rotate_drawing(drawing, 30)
clipped <- clip_circle(rotated, radius = 15)
segmented <- segment_at_junctions(drawing)
regions <- split_junction_regions(drawing)

# Optional file adapters:
# drawing <- read_mlv_svg("drawing.svg")
# drawings <- read_mlv_mat("drawings.mat")
# traced <- trace_contours("black_lines_on_white.png")
# write_mlv_png(drawing, "drawing.png")
```

Both ranked groups obey their requested whole-contour length budgets; a
middle contour can remain unassigned. Curvature means and integrated curvature
have separate names. Geometry-changing operations invalidate derived results.

Model scoring accepts an R model/predict method or a function of a feature table.
Bind the exact training table's measurement schema before prediction:

```r
features <- feature_table(drawing)
model <- bind_feature_model(function(z) z$length_1, features)
scores <- predict_contours(drawing, model)
```

The schema records predictor order, bin edges, weights, units and normalization.
Prediction rejects missing or incompatible training metadata before calling the
model. Use fixed `length_range` and `curvature_range` values when comparing
drawings with different canvas sizes. A model can instead receive its training
schema via `expected_spec = feature_schema(features)`. The explicit legacy
option `check_spec = FALSE` warns that compatibility was not checked.

trace_rgb() requires an explicitly supplied provider returning an ink mask.
No SAM model/backend is bundled; arbitrary edge detectors are not silently
substituted for SAM. MATLAB serialized regression models, general SVG paint
effects and MAT v7.3/HDF5 import are outside the current adapters.

For reproducible measurements, use a fixed specification:

```r
spec <- feature_spec(c("length", "curvature", "orientation"),
                     estimator = "arclength", curvature_scale = 2)
measured <- compute_contour_properties(drawing, spec = spec)
features <- feature_table(measured)
```

The default curvature estimator, `turning`, spreads each vertex turn over the
half-segments on either side of it, so every turn is counted once; `arclength`
smooths the turns at a stated scale. Both wrap closed contours. The published
estimator, `estimator = "mlv2023"`, counts the last vertex of every contour
twice and depends on how a curve is sampled; use it to reproduce published
values.
Physical measurement uses `units = "visual_degrees"` with an explicit
`pixels_per_degree`. New branch scores use `branch_estimator = "arclength"`.
See the [measurement contract](docs/FEATURE_SPEC.md) for definitions, migration,
metadata, and the limits of raster resolution invariance.

## References, storage and bugs

- [Upstream bug register](UPSTREAM_BUGS.md): stable IDs, executed reproducers,
  source findings and R regressions. No upstream reports have been submitted.
- [Source provenance](vendor/SOURCES.json): pinned source-only snapshots,
  retained-file checksums and reference-data custody.
- [Cleanup receipt](CLEANUP_RECEIPT.json): one verified shared collection at
  reference-data/mlv; about 2.63 GiB of logical duplicate/history storage removed.
- tools/restore_upstream.py restores full original checkouts to an explicit
  new destination. Large reference data and vendor snapshots are excluded from
  the R package tarball.

## Credits, citation and licenses

The original toolbox and methods are by **Dirk Bernhardt-Walther, Delaram
Farzanfar, Seohee Han, Morteza Rezanejad and the Bernhardt-Walther Lab**.
The R implementation and extensions are by **Bradley Buchsbaum**.
We also credit **Aravind Narayanan** for the Python AOF reference,
**Josep Mones i Teixidor and Carnë Draug** for the Octave morphology tables,
and **Peter Kovesi** for retained edge-linking reference helpers.
[Full source and contribution details](docs/ATTRIBUTION.md) distinguish these roles.

Run `citation("mlvtoolbox")` to cite the R implementation and original methods.
The R implementation is declared **GPL-3.0-or-later**; see [COPYING](COPYING).
[Upstream notices](inst/LICENSES) and [credits](inst/COPYRIGHTS) are also included
in the installed package.

Vendored sources retain their own terms: Python has an MIT root license with
GPL notices in three files, while many MATLAB files have **no explicit license
grant**. [The audit](vendor/LICENSING.md) records this unresolved permission
issue; our GPL declaration does not relicense those files. Dataset and model
redistribution rights are separate, and neither is asserted here.

## Development

Run testthat::test_local() for behavior/numerical checks. The reference generation
and upstream reproducer scripts are in tools/reference and tools/upstream.
Local check results, benchmarks and cross-platform status are recorded separately.

For a bounded source build/check, use `python3 tools/bounded_check.py`.
It stages only package files before invoking R, so retained reference data
are not copied into the build workspace.
