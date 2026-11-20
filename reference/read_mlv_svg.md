# Import vector contours from SVG

Supports lines, polylines, polygons, rectangles, circles/ellipses,
Bezier paths, elliptical arcs and affine transforms. The viewBox origin
maps to (1,1) in image coordinates. Curves are sampled according to
control-polygon length; curve_step is a sampling parameter, not a
certified Hausdorff bound. Text, embedded images, clipping, masks and
paint effects are not rendered.

## Usage

``` r
read_mlv_svg(file, width = NULL, height = NULL, curve_step = 1, strict = TRUE)
```

## Arguments

- file:

  SVG file.

- width, height:

  Optional output dimensions.

- curve_step:

  Positive curve sampling step.

- strict:

  Error on unsupported visible content; otherwise warn.

## Value

A drawing with import provenance.
