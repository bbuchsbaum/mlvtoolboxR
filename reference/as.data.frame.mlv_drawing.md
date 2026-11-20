# Drawing segment table

Drawing segment table

## Usage

``` r
# S3 method for class 'mlv_drawing'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  spec = x$feature_spec,
  pixels_per_degree = NULL,
  estimator = c("turning", "mlv2023")
)
```

## Arguments

- x:

  A drawing.

- row.names, optional, ...:

  Standard data frame method arguments.

- spec:

  Optional feature specification; defaults to an enriched drawing's
  recipe.

- pixels_per_degree:

  Optional input calibration override.

- estimator:

  Curvature definition without a specification; see
  \[compute_contour_properties()\]; exclusive of \`spec\`.

## Value

A data frame with one row per input segment: endpoints, length,
orientation and curvature. Under the turning and arc-length estimators
curvature is the segment's mean curvature (the exact integral of the
estimator's density over the segment divided by its length) and
\`signed_curvature\` its signed mean (NA where an exact reversal leaves
the sign undefined); \[curvature_profile()\] returns the piecewise
density itself. Lengths and curvature use the specification's units.
