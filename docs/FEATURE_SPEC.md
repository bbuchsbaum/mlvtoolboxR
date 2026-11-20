# Measurement contract for Release A

`feature_spec()` declares a fixed measurement recipe. Its format version and
method definition versions are checked on use. A prediction schema records the
resolved predictor columns; a feature specification supplies the recipe that
produces them. Binding a model requires its actual training table and never
infers unavailable training history.

The default contour curvature estimator is `turning` (0.3.0): each vertex turn
is spread once over its dual cell, with no tuning parameter. `arclength` gives
curvature at a stated scale, and `mlv2023` reproduces the published estimator,
which counts the last vertex of each contour twice (MLV-M032). Specification-free
calls take the same default and accept `estimator = "mlv2023"`; their feature
schemas record the curvature definition in effect. The default branch estimator
remains `mlv2023`.
Legacy numerical parity fixtures remain unchanged. Metadata fixtures evolve with
explicit metadata assertions. The definition versions describe this port and
its documented deviations, not universal equivalence to MATLAB.

| Settings | Meaning and units |
|---|---|
| properties, bins, ranges | Ordered predictors and fixed histogram definitions; orientation uses axial upper-edge wrapping; other numeric bins clamp overflow; a value within 1e-9 × max(1, \|edge\|) of a bin edge (in the binned coordinate) is binned as if on the edge |
| units | Spatial pixels or degrees of visual angle; tangent and turn angles remain geometric degrees |
| pixels_per_degree | Required positive input calibration for visual degrees; may vary between inputs with otherwise compatible physical definitions |
| normalize | Length/count weights per 10,000 target spatial units; unnormalized outputs retain their natural weight units |
| estimator, curvature_scale | Published, vertex turning, or local arc-length curvature; scale is a physical half-support |
| branch_estimator, branch_window | Published sample windows or new physical arc-length half-windows |
| branch_power, reference_radius, fit_tolerance | Score exponent and explicitly dimensioned reference lengths; published scores retain their historical one-pixel rules |
| junction | Extension, merging, and endpoint radii in target spatial units; relative extension is dimensionless |
| mat | AOF angle in degrees, angular sample count, epsilon in raster pixels, component/branch sizes in raster samples |
| mapping_radius, raster_line_width | Pixel-grid operations, explicitly measured in pixels |
| extraction | Trace method, threshold, thinning, minimum point count, simplification tolerance in pixels, and optional provider plus version |

Raster settings remain pixel-grid definitions. Equal physical contour geometry
does not imply equal raster-derived scores at different raster resolutions.
Those differences belong in sensitivity tests and are retained in the schema.
There is no calibration inference from canvas dimensions. Existing source
provenance is retained separately; declaring an extraction recipe does not prove
that pre-existing contours were produced with that recipe.

Changing a requested measurement definition requires recomputation; input hashes,
timestamps and calibration do not by themselves make converted physical columns
incompatible. Requested extraction settings and observed provenance must remain
distinguishable. Unsupported method versions fail rather than falling back.

## Measurement and prediction

```r
spec <- feature_spec(c("length", "curvature", "orientation"),
  estimator = "arclength", curvature_scale = 0.2,
  units = "visual_degrees", pixels_per_degree = 40)
measured <- compute_contour_properties(drawing, spec = spec)
training <- feature_table(measured)
# Bind the table actually used to train your model:
bound <- bind_feature_model(fitted_model, training)
scores <- predict_contours(other_drawing, bound, spec = spec,
                           pixels_per_degree = 60)
```

`contour_statistics()`, `feature_table()`, `average_property()` and the drawing's
data-frame method use an enriched drawing's saved specification. Supply either
`spec` or individual legacy settings. Explicit `spec = NULL` selects the legacy
path. Geometry transformations retain the recipe but remove derived results.
MAT data retain a geometry snapshot and the relevant measurement settings;
incompatible caches fail with instructions to recompute. `analyse_mat(x, spec =
spec)` returns a list containing the scored `mat` and enriched `drawing`.

`batch_statistics()` requires compatible schemas. Feature tables and long
statistics carry `measurement_metadata`, including the schema, input calibration,
observed extraction provenance and mapped/total pixel coverage for MAT features.
`write_feature_csv()` writes `<file>.schema.rds` alongside the CSV. Keep these
files together. `saveRDS()` preserves complete drawings and bound models.
RGB provider identity/version is declared explicitly; the package cannot verify
the behavior of an arbitrary provider function or infer its model weights.

Arc-length curvature is a source-normalized tent smoothing of signed vertex
turns. Its histogram integrates the affine density within each interval, rather
than assigning an interval's entire length to its mean-value bin. Opposite turns
can cancel at the chosen scale. Use `contour_turning()` for unsmoothed signed and
absolute total turns; use `curvature_profile()` for the local scale-dependent
descriptor. Closed contours wrap and a single elbow contributes one 90-degree
turn. Clipped open kernels preserve signed mass but have boundary bias. Raster
tracing and simplification sensitivity have separate acceptance tolerances.

The new branch separation score is `mean_radius / (mean_radius + reference_radius)`
before exponentiation; it is a different, bounded definition from the published
`1 - 1/R`. Parallelism integrates the radius graph over a physical window, mirror
combines chord efficiency and perpendicular residual, and taper measures the
residual from an affine radius fit. These are geometric descriptors, not newly
validated perceptual scales. All four definitions have explicit version tags.

MAT histograms still count finite mapped raster pixels, including zero scores;
unmapped pixels are excluded and reported in coverage. They are not physical
arc-length histograms and are not resolution invariant. Rendering remains a
pixel-grid visualization. New-spec curvature density maps distribute integrated
absolute turning; numeric histograms are the authoritative quantitative output.

For legacy models, regenerate the original training feature table with its
original settings, then call `bind_feature_model(model, training)`. If those
settings are unavailable, the package cannot establish compatibility. The
warning-emitting `check_spec = FALSE` bypass is explicit and is not a validation.
