# Average outward flux of a signed-distance/nearest-boundary field

Average outward flux of a signed-distance/nearest-boundary field

## Usage

``` r
average_outward_flux(distance, index, samples = 60L, epsilon = 1)
```

## Arguments

- distance:

  Numeric signed-distance matrix.

- index:

  Nearest-boundary indices in R column-major order.

- samples:

  Number of angular samples.

- epsilon:

  Sampling-circle radius.

## Value

Numeric flux matrix, without mutating the inputs.
