# Exact signed and absolute vertex turning

Exact signed and absolute vertex turning

## Usage

``` r
contour_turning(
  x,
  units = c("pixels", "visual_degrees"),
  pixels_per_degree = NULL
)
```

## Arguments

- x:

  A drawing.

- units:

  Spatial units for arc positions; angles remain geometric degrees.

- pixels_per_degree:

  Required calibration for visual-degree positions.

## Value

A table of canonical vertices, positions and turns. x/y coordinates
remain in input pixels. Collinear forward vertices are removed; exact
retracing is rejected because its signed turn is undefined.
