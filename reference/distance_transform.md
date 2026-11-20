# Exact Euclidean distance to the nearest feature pixel

Exact Euclidean distance to the nearest feature pixel

## Usage

``` r
distance_transform(features)
```

## Arguments

- features:

  Binary matrix; TRUE marks a feature.

## Value

List with distance and nearest one-based column-major index matrices.
Empty feature sets produce Inf distances and NA indices. Ties select the
smallest linear index.
