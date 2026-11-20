# Plot a stage of a medial-axis representation

Draws one stage of the pipeline stored in an \`mlv_mat\`, in image
coordinates: x increases to the right, y increases downward, pixel
centres lie at integers and the aspect ratio is one. Quantities are
colour-coded with a legend that names the quantity and its unit:

## Usage

``` r
# S3 method for class 'mlv_mat'
plot(
  x,
  what = "flux",
  drawing = NULL,
  legend = TRUE,
  zlim = NULL,
  col = NULL,
  main = NULL,
  axes = FALSE,
  overlay_col = "black",
  overlay_lwd = 0.6,
  ...
)
```

## Arguments

- x:

  An \`mlv_mat\` from \[compute_mat()\].

- what:

  Stage to draw; see Description.

- drawing:

  Optional \`mlv_drawing\` whose segments are drawn on top, for example
  the drawing the MAT was computed from. Its dimensions must match.

- legend:

  Draw a colour bar or key to the right of the image.

- zlim:

  Optional colour limits for \`"distance"\`, \`"flux"\` and scores. The
  default flux limits are symmetric about zero.

- col:

  Optional vector of colours replacing the default palette for a
  continuous stage.

- main:

  Panel title; the default names the stage. Titles use \`cex.main = 1\`
  and \`font.main = 1\` unless set in \`...\`.

- axes:

  Draw light pixel axes.

- overlay_col, overlay_lwd:

  Colour and width of the \`drawing\` overlay.

- ...:

  Further graphical parameters passed to \[graphics::plot.window()\] and
  \[graphics::title()\], such as \`cex.main\`.

## Value

\`x\`, invisibly.

## Details

\* \`"ink"\`: the rasterized ink mask (black on white). \*
\`"distance"\`: the signed distance map in pixels (sequential palette).
Ink pixels have nonpositive values. \* \`"flux"\`: the average outward
flux, stored as the sum of \`samples\` dot products (diverging palette
with limits symmetric about zero). The legend marks the skeleton
threshold, \`threshold \* samples\`. \* \`"skeleton"\`: the pixels whose
flux exceeds the threshold (light) and the thinned skeleton (dark), with
the ink in grey. \* \`"branches"\`: traced branches, one colour per
branch, each labelled with its index in \`x\$branches\`. \* a branch
score name present in \`x\$score_images\` (for example
\`"parallelism"\`; see \[compute_mat_properties()\]): the score at each
skeleton pixel on a 0-1 sequential palette. Pixels without a score are
not drawn.

The function changes no persistent graphics parameter: it reserves room
for the legend inside the plotting region, so later calls such as
\[graphics::points()\] use the same image coordinates.

## Examples

``` r
d <- mlv_drawing(list(matrix(c(5, 10, 45, 10), 1),
                      matrix(c(5, 22, 45, 22), 1)), 50, 32)
m <- compute_mat(d)
plot(m, "flux", drawing = d)

plot(m, "branches")
```
