# Split contours using weighted feature histograms

Split contours using weighted feature histograms

## Usage

``` r
split_by_histogram(
  x,
  weights,
  properties = c("orientation", "length", "curvature"),
  bins = 8,
  fraction = 0.5,
  spec = NULL
)
```

## Arguments

- x:

  A drawing.

- weights:

  Named numeric vector matching feature-table columns (other than
  \`contour_id\`), each named at most once.

- properties, bins:

  Feature settings, used when \`spec\` is NULL; the histograms are then
  unnormalized.

- fraction:

  Whole-contour length budget.

- spec:

  Optional feature specification, exclusive of \`properties\` and
  \`bins\`; its own histogram settings (including normalization) apply.
  The default ignores any specification the drawing carries.

## Value

The same structure as split_contours.
