# Trace RGB images with an explicit extraction provider

Trace RGB images with an explicit extraction provider

## Usage

``` r
trace_rgb(image, backend, ..., min_points = 2L, spec = NULL)
```

## Arguments

- image:

  Numeric RGB array or image filename.

- backend:

  Function returning a logical ink mask from the image.

- ...:

  Arguments passed to the provider.

- min_points:

  Minimum traced path size.

- spec:

  Optional specification, including explicit provider name/version.

## Value

A drawing. Model weights and preprocessing belong to the provider.
