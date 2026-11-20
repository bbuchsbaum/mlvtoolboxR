# Bind a prediction model to its training feature schema

Bind using the exact feature table used for training. The wrapper does
not train the model or infer its training history. It retains the schema
through R serialization without relying on a model to preserve
data-frame attributes.

## Usage

``` r
bind_feature_model(model, features)
```

## Arguments

- model:

  An R prediction model or function accepting a predictor table.

- features:

  The training feature table, or its \`feature_schema()\`.

## Value

An \`mlv_feature_model\` containing the model and training schema.
