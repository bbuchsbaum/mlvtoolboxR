# Render a drawing with property or model-score colors

Render a drawing with property or model-score colors

## Usage

``` r
render_drawing(
  x,
  property = NULL,
  values = NULL,
  width = x$width,
  height = x$height,
  line_width = 1,
  palette = grDevices::hcl.colors(256, "Viridis")
)
```

## Arguments

- x:

  A drawing.

- property:

  Optional segment or mapped medial-axis property name.

- values:

  Optional numeric score for each contour, for example model
  predictions; exclusive of \`property\`.

- width, height:

  Output dimensions.

- line_width:

  Stroke width.

- palette:

  Color vector used from low to high values.

## Value

RGB array with white background. Missing scores are gray.
