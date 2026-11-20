# Render a drawing into an ink mask

Uses deterministic clipped Bresenham lines for width one. Wider strokes
include pixels whose centers lie within half the stroke width of a
segment. This rasterizer is explicit and is not claimed equivalent to
insertShape.

## Usage

``` r
rasterize(x, width = x$width, height = x$height, line_width = 1)
```

## Arguments

- x:

  A drawing.

- width, height:

  Output dimensions.

- line_width:

  Positive stroke width in output pixels.

## Value

Logical matrix in row/column order, TRUE for ink.
