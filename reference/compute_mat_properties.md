# Compute medial-axis symmetry scores

Compute medial-axis symmetry scores

## Usage

``` r
compute_mat_properties(
  mat,
  properties = c("parallelism", "separation", "mirror"),
  window = 5L,
  power = 10,
  spec = mat$feature_spec,
  pixels_per_degree = NULL
)
```

## Arguments

- mat:

  An mlv_mat object.

- properties:

  Requested branch score names.

- window, power:

  Branch score parameters.

- spec:

  Optional feature specification, exclusive of individual settings.

- pixels_per_degree:

  Optional input calibration override.

## Value

An mlv_mat with branch scores and score images; unmapped pixels are NA.

## Details

Under \`branch_rule="invariant"\` (default) a branch score is NA where
its estimator is undefined (see \[score_branch()\]). The scores of the
two traversal directions are averaged over the defined ones; a score
image pixel is the maximum of the defined scores at that pixel over all
branches (and over both occurrences of a junction pixel closing a loop),
NA if none is defined. For contour mapping see \[map_mat_properties()\].
Under \`branch_rule="mlv2023"\` undefined points carry the published
zeros.
