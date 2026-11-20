# Wide per-contour feature table

Wide per-contour feature table

## Usage

``` r
feature_table(
  x,
  properties = c("orientation", "length", "curvature", "junction_type"),
  bins = 8,
  normalize = TRUE,
  length_range = c(2, max(3, x$width + x$height)),
  curvature_range = c(0, 90),
  spec = x$feature_spec,
  pixels_per_degree = NULL,
  estimator = c("turning", "mlv2023")
)
```

## Arguments

- x:

  A drawing.

- properties:

  Properties to include.

- bins:

  Number of bins.

- normalize:

  Normalize by contour length.

- length_range, curvature_range:

  Fixed ranges for logarithmic bins.

- spec:

  Optional feature specification, exclusive of individual settings.

- pixels_per_degree:

  Optional input calibration override.

- estimator:

  Curvature definition without a specification; see
  \[contour_statistics()\].

## Value

One row per contour, with named numeric predictors and a
\`feature_schema\` attribute. Inspect it with \`feature_schema()\`.
