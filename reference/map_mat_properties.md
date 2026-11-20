# Map branch ratings back to drawing contours

Uses exact nearest tangent samples per branch, then the maximum across
contributing branches. An explicit mapping mask distinguishes valid zero
ratings from missing samples. Geometry is reused across requested
properties.

## Usage

``` r
map_mat_properties(
  x,
  mat,
  radius = 4,
  spec = mat$feature_spec,
  pixels_per_degree = NULL
)
```

## Arguments

- x:

  A drawing.

- mat:

  An mlv_mat with computed branch properties.

- radius:

  Maximum tangent mapping distance in pixels.

- spec:

  Optional feature specification retained by the scored MAT.

- pixels_per_degree:

  Optional input calibration override.

## Value

A drawing with per-contour MAT scores, score images and mapping masks.

## Details

Under the default \`branch_rule="invariant"\`, each tangent sample of
segment (i, i+1) carries the mean of the defined scores at i and i+1 (NA
if both are undefined). For each branch, an ink pixel takes its nearest
tangent samples of that branch (all samples within 1e-9 px of the
minimum distance) if that distance is below \`radius\`; its value from
the branch is the maximum of the defined scores among those nearest
samples, or NA if all of them are undefined. A farther defined sample of
the same branch is not used. The pixel's score is the maximum over
branches of these values, NA if no branch gives a defined value. Under
\`branch_rule="mlv2023"\` the published first-point, lowest-index and
zeros-as-missing rules apply.
