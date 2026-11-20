# Trace contours from an ink mask or grayscale line drawing

The graph uses eight-neighbor links but suppresses redundant diagonal
shortcuts across an existing orthogonal link. Junctions split paths and
closed loops retain their closure. This is an explicit graph method,
rather than an unqualified emulation of Kovesi edge linking.

## Usage

``` r
trace_contours(
  image,
  threshold = 0.5,
  min_points = 2L,
  thin = TRUE,
  simplify = NULL,
  spec = NULL
)
```

## Arguments

- image:

  Binary ink mask, grayscale matrix, RGB array, or image filename.

- threshold:

  Values below this threshold are ink in numeric images.

- min_points:

  Minimum traced path size.

- thin:

  Thin the mask before tracing (the default \`"minimal"\` thinning, or
  the specification's MAT morphology).

- simplify:

  Optional nonnegative line-fit distance tolerance.

- spec:

  Optional feature/extraction specification, exclusive of individual
  extraction settings. The resulting drawing retains the specification.

## Value

A drawing with extraction provenance.
