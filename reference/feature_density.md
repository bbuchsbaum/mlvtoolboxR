# Feature density on a pixel grid

Feature density on a pixel grid

## Usage

``` r
feature_density(x, property = "length", sigma = 0)
```

## Arguments

- x:

  A drawing.

- property:

  Segment property, junctions, or a mapped MAT property.

- sigma:

  Gaussian smoothing standard deviation in pixels; zero disables.

## Value

Numeric image matrix. Kernel boundaries are zero padded. Curvature
distributes integrated absolute turning (curvature times length of each
measured piece) over raster pixels, under the drawing's specification
or, without one, the vertex turning estimator.
