# Predict one score per contour

Predict one score per contour

## Usage

``` r
predict_contours(
  x,
  model,
  properties = c("orientation", "length", "curvature", "junction_type"),
  bins = 8,
  ...,
  expected_spec = NULL,
  check_spec = TRUE,
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

- model:

  An R model with a predict method, or a function of a feature table.

- properties, bins:

  Feature-table settings.

- ...:

  Additional model prediction arguments.

- expected_spec:

  Training schema, from \`feature_schema()\`. Required unless the model
  was bound with \`bind_feature_model()\`.

- check_spec:

  Check training/actual measurement compatibility. Setting FALSE is an
  explicit, warning-emitting legacy bypass.

- normalize, length_range, curvature_range:

  Feature-table settings.

- spec:

  Optional feature specification, exclusive of individual settings.

- pixels_per_degree:

  Optional input calibration override.

- estimator:

  Curvature definition without a specification; see
  \[contour_statistics()\].

## Value

A finite numeric vector in contour order.
