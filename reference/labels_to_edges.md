# Convert segmentation labels to an ink boundary mask

Convert segmentation labels to an ink boundary mask

## Usage

``` r
labels_to_edges(labels)
```

## Arguments

- labels:

  Finite integer-valued label matrix.

## Value

Logical matrix marking both sides of horizontal/vertical label
transitions. Outer image borders are not added automatically.
