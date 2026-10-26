# mlvtoolbox

An R package for measuring and manipulating contours in vectorized line
drawings. It computes orientation, length, curvature, junctions and AOF-based
medial-axis properties, with explicit image-coordinate and numerical conventions.

The implementation follows the [MLV methods paper](https://doi.org/10.3389/fcomp.2023.1140723).
It corrects demonstrated boundary, geometry and missingness defects; historical
outputs are comparison evidence rather than the definition of correctness.

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
events. Normalized values are per 10,000 contour pixels. Valid zero scores
remain distinct from missing scores.

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

This is a qualified implementation with explicit methods: compiled exact
Euclidean distance transforms, a deterministic nearest-boundary tie rule,
Octave-image Pratt morphology and a pixel-graph tracing method. The Pratt LUTs
and AOF have small executed Octave comparison fixtures. These checks do not
establish equivalence to MATLAB's proprietary morphology or rasterizer, nor
broad natural-image scientific acceptance. Executed comparisons show that
equidistant nearest-boundary choices can change flux and skeleton pixels between
R and Octave, despite identical distances; see [parity evidence](docs/PARITY.md).
Taper uses documented explicit
smoothing. See [implementation status](IMPLEMENTATION_STATUS.md).

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

`mlv2023` remains the default for compatibility. The new `turning` and
`arclength` estimators count each vertex turn once and wrap closed contours.
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

The package and attributed Pratt/AOF-derived code use GPL-3.0-or-later. Original
reference snapshots retain their own notices. Dataset redistribution is not
asserted; package tests use synthetic inputs.

Run testthat::test_local() for behavior/numerical checks. The reference generation
and upstream reproducer scripts are in tools/reference and tools/upstream.
Local check results, benchmarks and cross-platform status are recorded separately.

For a bounded source build/check, use `python3 tools/bounded_check.py`.
It stages only package files before invoking R, so retained reference data
are not copied into the build workspace.
