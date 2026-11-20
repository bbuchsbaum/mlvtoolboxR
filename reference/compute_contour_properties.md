# Compute contour measurements and junctions

Compute contour measurements and junctions

## Usage

``` r
compute_contour_properties(
  x,
  properties = c("orientation", "length", "curvature", "junctions"),
  ...,
  spec = x$feature_spec,
  pixels_per_degree = NULL,
  estimator = c("turning", "mlv2023")
)
```

## Arguments

- x:

  A drawing.

- properties:

  Requested properties; basic segment measures are cheap and computed
  together. Include "junctions" to detect contour intersections.

- ...:

  Arguments passed to \`detect_junctions()\`.

- spec:

  Optional feature specification, exclusive of properties, dots and
  \`estimator\`.

- pixels_per_degree:

  Optional input calibration override.

- estimator:

  Curvature definition without a specification. Both give one row per
  input segment. \`"turning"\` (default) spreads each vertex turn over
  its dual cell, counting it once; a segment's curvature is its mean
  curvature, the exact integral of that density over the segment divided
  by its length, and \`signed_curvature\` is the signed mean. The
  piecewise density is available from \[curvature_profile()\].
  \`"mlv2023"\` gives the published value: the turn to the next segment
  divided by the segment's own length, with the last segment repeating
  the previous turn. Without a specification the choice is recorded on
  the drawing (\`curvature_estimator\`) and becomes the default of later
  calls on it that give neither \`spec\` nor \`estimator\` (statistics,
  feature tables, segment tables, plots, densities and contour splits);
  explicit arguments still win.

## Value

A drawing with segment and contour tables, and requested junctions.
