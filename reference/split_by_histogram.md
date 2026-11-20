# Split contours using weighted feature histograms

Split contours using weighted feature histograms

## Usage

``` r
split_by_histogram(
  x,
  weights,
  properties = c("orientation", "length", "curvature"),
  bins = 8,
  fraction = 0.5
)
```

## Arguments

- x:

  A drawing.

- weights:

  Named numeric vector matching feature-table columns.

- properties, bins:

  Feature settings.

- fraction:

  Whole-contour length budget.

## Value

The same structure as split_contours.
