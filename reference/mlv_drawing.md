# Construct a vectorized line drawing

Contours are ordered matrices with columns x1, y1, x2, y2. Coordinates
are continuous image coordinates: x increases right, y increases down;
raster pixel centers begin at (1,1). Geometry may extend outside the
image. Zero-length segments require explicit cleaning; a segment whose
endpoint coordinates differ by at most 64 units in the last place of
their own magnitude (rounding residue) counts as zero-length. When such
segments are dropped, the next kept segment starts exactly where the
previous one ends (cyclically for closed contours). Subpixel values are
otherwise retained.

## Usage

``` r
mlv_drawing(
  contours = list(),
  width,
  height,
  ids = NULL,
  source = NULL,
  zero_length = c("error", "drop")
)
```

## Arguments

- contours:

  List of numeric matrices with four columns.

- width, height:

  Positive integer image dimensions.

- ids:

  Unique character contour identifiers.

- source:

  Optional provenance metadata.

- zero_length:

  Reject or drop zero-length segments.

## Value

An S3 \`mlv_drawing\` with geometry and provenance.

## Examples

``` r
d <- mlv_drawing(list(rbind(c(1,1,11,1), c(11,1,11,6))), 20, 20)
as.data.frame(d)
#>   contour_id contour_index segment_index x1 y1 x2 y2 length orientation
#> 1          1             1             1  1  1 11  1     10           0
#> 2          1             1             2 11  1 11  6      5         270
#>   curvature signed_curvature
#> 1         6               -6
#> 2         6               -6
```
