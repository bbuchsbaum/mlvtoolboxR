# Export image or contour statistics

Export image or contour statistics

## Usage

``` r
write_feature_csv(x, file, ...)
```

## Arguments

- x:

  A drawing.

- file:

  CSV destination.

- ...:

  Arguments to contour_statistics.

## Value

The written data frame invisibly. A companion file named
\`\<file\>.schema.rds\` preserves the measurement schema, calibration
and provenance.
