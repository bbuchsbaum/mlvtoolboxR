# Greedy line fitting for ordered points

Greedy line fitting for ordered points

## Usage

``` r
fit_line_segments(xy, tolerance = 1)
```

## Arguments

- xy:

  Finite two-column point matrix.

- tolerance:

  Maximum mean interior perpendicular distance.

## Value

List with endpoint-index spans and per-point straightness scores.
