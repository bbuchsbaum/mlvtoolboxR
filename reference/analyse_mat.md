# Complete vector drawing to medial-axis feature workflow

Complete vector drawing to medial-axis feature workflow

## Usage

``` r
analyse_mat(
  x,
  properties = c("parallelism", "separation", "mirror"),
  ...,
  spec = x$feature_spec,
  pixels_per_degree = NULL
)
```

## Arguments

- x:

  A drawing.

- properties:

  MAT property names.

- ...:

  Arguments to compute_mat.

- spec:

  Optional feature specification, exclusive of properties and dots.

- pixels_per_degree:

  Optional input calibration override.

## Value

List containing the enriched drawing and medial-axis representation.
