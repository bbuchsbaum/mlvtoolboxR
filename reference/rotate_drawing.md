# Rotate a drawing in image coordinates

Rotate a drawing in image coordinates

## Usage

``` r
rotate_drawing(x, angle, center = c(x$width, x$height)/2)
```

## Arguments

- x:

  A drawing.

- angle:

  Counterclockwise angle in degrees.

- center:

  Two-coordinate rotation center.

## Value

A drawing with derived properties invalidated.
