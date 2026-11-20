# Inspect the measurement schema of a feature table

A schema records ordered predictor names, actual bin boundaries,
weights, units, normalization and the measurement definition version.
Canvas size is represented only through the resulting bin boundaries,
not as a separate compatibility condition. Row subsets retain the
schema; renamed, reordered or removed predictor columns are rejected.

## Usage

``` r
feature_schema(x)
```

## Arguments

- x:

  A feature table returned by \`feature_table()\`.

## Value

An \`mlv_feature_schema\`, suitable for binding or checked prediction.
