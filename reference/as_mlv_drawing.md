# Convert an upstream MATLAB drawing representation

Convert an upstream MATLAB drawing representation

## Usage

``` r
as_mlv_drawing(x, zero_length = c("error", "drop"))
```

## Arguments

- x:

  A named list with contours and imsize fields.

- zero_length:

  Reject or explicitly drop zero-length segments.

## Value

An mlv_drawing.
