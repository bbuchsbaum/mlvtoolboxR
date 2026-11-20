# Shift contours, optionally sampling translations within the image

Shift contours, optionally sampling translations within the image

## Usage

``` r
shift_contours(x, offsets = NULL, within = is.null(offsets))
```

## Arguments

- x:

  A drawing.

- offsets:

  A contour-by-two translation matrix, or NULL for random shifts.

- within:

  Require translated contours to fit inside the image.

## Value

A drawing; R's ordinary random seed controls random translations.
