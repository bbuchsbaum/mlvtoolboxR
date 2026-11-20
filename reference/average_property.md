# Mean geometric or medial-axis property

Mean geometric or medial-axis property

## Usage

``` r
average_property(
  x,
  property = c("length", "curvature", "orientation", "junctions", "mirror",
    "parallelism", "separation", "taper"),
  spec = x$feature_spec,
  pixels_per_degree = NULL,
  estimator = c("turning", "mlv2023")
)
```

## Arguments

- x:

  A drawing, with medial-axis properties when requested.

- property:

  Property name.

- spec:

  Optional feature specification.

- pixels_per_degree:

  Optional input calibration override.

- estimator:

  Curvature definition without a specification; see
  \[compute_contour_properties()\]. Mean curvature is integrated
  absolute turning per unit length.

## Value

One numeric summary, or NA when undefined.
