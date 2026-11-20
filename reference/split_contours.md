# Split ranked contours under a whole-contour length budget

Each group is a ranked prefix from one end that fits within fraction
times total contour length. A contour that would exceed the budget stops
that prefix. Ties use contour IDs. Curvature ranks by a length-weighted
mean; integrated_curvature explicitly exposes the upstream
length-weighted sum. Medial-axis criteria rank by the mean of a
contour's finite mapped scores. A contour with no finite mapped score
(no defined branch score reaches any of its pixels) cannot be ranked: it
is left unassigned, its score is NA, and it is listed in \`unranked\`.
Ranks are taken among the rankable contours only; the length budget
remains a fraction of the total length of all contours.

## Usage

``` r
split_contours(
  x,
  property = "length",
  fraction = 0.5,
  weights = NULL,
  scores = NULL
)
```

## Arguments

- x:

  A drawing.

- property:

  One or more criterion names.

- fraction:

  Length fraction in \[0,1\] for each group.

- weights:

  Criterion rank weights.

- scores:

  Optional finite per-contour scores, bypassing rank combination.

## Value

List with top, bottom, unassigned drawings, overlap IDs, scores (NA for
unrankable contours) and \`unranked\`, the IDs of contours without a
criterion value (always empty when \`scores\` is supplied).
