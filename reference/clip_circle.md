# Clip all segment-circle intersections

Handles segments with both endpoints outside but crossing the aperture.
Adjacent retained segments remain in the same contour.

## Usage

``` r
clip_circle(
  x,
  radius = min(x$width, x$height)/2,
  center = c(x$width, x$height)/2
)
```

## Arguments

- x:

  A drawing.

- radius:

  Positive aperture radius.

- center:

  Aperture center.

## Value

A drawing with parent contour IDs and no derived-property cache.
