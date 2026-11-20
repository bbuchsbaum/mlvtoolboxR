# Plot vector contours or properties

Draws the drawing on its canvas in image coordinates (x rightward, y
downward, equal aspect ratio). The canvas spans the raster extent, 0.5
to \`width + 0.5\` and 0.5 to \`height + 0.5\`, because pixel centres
lie at integers starting at 1; the plot coordinates are those of
\[plot.mlv_mat()\], so MAT stages and drawings overlay exactly. With a
\`property\`, segments are coloured by that quantity and, by default, a
colour bar naming the quantity and its unit is drawn to the right of the
canvas, inside the plot region, so later low-level calls (\`points()\`,
\`text()\`) use the same coordinates.

## Usage

``` r
# S3 method for class 'mlv_drawing'
plot(
  x,
  property = NULL,
  col = "black",
  lwd = 1,
  ...,
  legend = TRUE,
  palette = NULL,
  zlim = NULL,
  arms = FALSE,
  arm_length = NULL,
  junction_col = "#D55E00",
  frame = TRUE
)
```

## Arguments

- x:

  A drawing.

- property:

  Optional segment property (\`"orientation"\`, \`"length"\`,
  \`"curvature"\`), \`"junctions"\`, or a mapped MAT property.

- col:

  Contour colour when no segment property is shown.

- lwd:

  Line width.

- ...:

  Arguments to the initial plot, such as \`main\`, \`xlim\` or \`axes =
  FALSE\`.

- legend:

  Draw a colour bar when colour encodes a quantity.

- palette:

  Optional palette: a function of the number of colours, or a character
  vector of colours. Defaults are described above.

- zlim:

  Value range mapped onto the palette for length and curvature; values
  outside it take the end colours. Defaults to the finite data range.
  Ignored for orientation (always 0 to 180 degrees) and MAT scores (0 to
  1).

- arms:

  Draw junction arms.

- arm_length:

  Drawn arm length in pixels; defaults to 8 canvas side.

- junction_col:

  Colour of junction markers and arms.

- frame:

  Outline the canvas in light grey.

## Value

The drawing invisibly.

## Details

Orientation is axial: a segment at 10 degrees and one at 190 degrees
have the same colour, and the default palette is a cyclic hue wheel over
0 to 180 degrees, so 0 and 180 share a colour. Length and curvature use
the sequential Viridis palette over \`zlim\`. Mapped medial-axis scores
are drawn as filled mapped pixels on a 0 to 1 Viridis scale over grey
contours. With \`property = "junctions"\`, or \`arms = TRUE\`, junction
positions are marked; \`arms = TRUE\` also draws each arm as an arrow
leaving the junction along its orientation (the \`arms\` column of
\[detect_junctions()\]). Junctions come from the drawing's stored
\`junctions\` (for example from \[compute_contour_properties()\] with
\`rule = "mlv2023"\`), its feature specification, or default
\[detect_junctions()\].

## Examples

``` r
d <- mlv_drawing(list(rbind(c(5, 10, 35, 10)), rbind(c(20, 10, 20, 35))), 40, 40)
plot(d, property = "orientation")

plot(d, property = "junctions", arms = TRUE, axes = FALSE)
```
