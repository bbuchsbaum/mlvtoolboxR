# Plot several property panels

Plot several property panels

## Usage

``` r
plot_properties(
  x,
  properties = c("orientation", "length", "curvature", "junctions"),
  columns = 2L,
  ...
)
```

## Arguments

- x:

  A drawing.

- properties:

  Properties passed to plot.

- columns:

  Panel column count.

- ...:

  Additional drawing plot arguments.

## Value

The drawing invisibly; graphics parameters are restored.
