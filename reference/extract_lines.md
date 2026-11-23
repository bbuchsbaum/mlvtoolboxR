# Extract a line drawing from an image

Turns a photograph or an image of a painting into a binary ink mask
whose lines follow the image's edges, ready for \[trace_contours()\] or,
through \[xdog_lines()\], for \[trace_rgb()\]. The only method,
\`"fdog"\`, is a flow-based difference of Gaussians with an XDoG
threshold, in four stages:

## Usage

``` r
extract_lines(
  image,
  method = "fdog",
  sigma = 2,
  k = 1.6,
  tau = 0.98,
  eps = 0,
  phi = Inf,
  sigma_c = 2,
  sigma_m = 3,
  smooth = 5,
  colour = FALSE,
  min_component = 20L,
  junction_gap = 0
)
```

## Arguments

- image:

  Greyscale matrix or RGB(A) array with values in \\0, 1\\
  (sRGB-encoded, as returned by \[read_mlv_image()\]), or a PNG/JPEG
  file name. Transparent pixels are composited over white.

- method:

  Extraction method; only \`"fdog"\` (matched exactly).

- sigma:

  Scale of the narrow Gaussian across the flow, in pixels (0.1 to 50);
  sets the scale of the edges that are drawn.

- k:

  Ratio of the wide to the narrow Gaussian scale (greater than 1).

- tau:

  Weight of the wide Gaussian, in \\0, 1\\; lower values draw only
  higher-contrast edges.

- eps:

  XDoG threshold on the response, in lightness units (\\-1, 1\\).

- phi:

  XDoG threshold sharpness (positive; \`Inf\` is a hard threshold).

- sigma_c:

  Gaussian scale of the structure tensor, in pixels (0.1 to 50).

- sigma_m:

  Scale of the integration along the flow, in pixels: 0 skips it,
  otherwise 0.1 to 50.

- smooth:

  Radius of the anisotropic Kuwahara pre-filter in pixels: 0 skips it,
  otherwise 0.5 to 50.

- colour:

  Use CIE \`L\*a\*b\*\` colour for the pre-filter and the flow.

- min_component:

  Smallest eight-connected ink component kept, in pixels. Narrow light
  gaps between dark marks are drawn as short centre lines, so lowering
  it below the default adds many short marks in textured areas.

- junction_gap:

  Optional maximum bridge length in pixels (0 to 50). Zero disables
  repair. The repair considers gaps longer than 1.5 pixels.

## Value

A logical ink mask (\`TRUE\` is ink) with the dimensions of the image
and an \`"extraction"\` attribute giving the method and validated
parameters. The original definition, \`"fdog-xdog-1"\`, omits
\`junction_gap\`. With repair enabled, \`"fdog-xdog-repair-1"\` records
it.

## Details

1\. \*\*Lightness and pre-filter.\*\* The image is converted from sRGB
to CIE \`L\*\` (D65 white), divided by 100. An anisotropic Kuwahara
filter (Kyprianidis, Kang and Doellner 2009) removes brush texture while
keeping form boundaries: each pixel's neighbourhood is an ellipse along
the local edge, of radius \`smooth\` (stretched to \`smooth \* (1 + A)\`
along the edge and shrunk to \`smooth / (1 + A)\` across it, \`A\` the
anisotropy), split into 8 overlapping polynomial sectors (Kyprianidis
2010) with a fall-off \`1 - 0.5426 r^2\` (\`r\` the radius in the unit
disc, close to the reference Gaussian). The output averages the sector
means with weights \`1 / (1 + (v / 0.01^2)^4)\`, \`v\` a sector's
lightness variance, so the most uniform sectors dominate. \`smooth = 0\`
skips the filter. Like any Kuwahara filter it rounds sharp corners
slightly (at the default radius the traced corner of a square lies up to
about 1.2 pixels inside the true corner, against 0.85 pixels without the
filter), and it removes dark strokes narrower than about 2 pixels unless
they lie along a pixel row or column; light strokes survive from 1
pixel. For thin line art, use \`smooth = 0\`. 2. \*\*Edge tangent
flow.\*\* The structure tensor of the smoothed image (Sobel gradients,
Gaussian scale \`sigma_c\`) gives the local edge tangent as its minor
eigenvector (Kyprianidis and Doellner 2008). 3. \*\*Flow-based DoG.\*\*
At every pixel a 1-D difference of Gaussians \`G_sigma - tau \* G\_(k
sigma)\` is taken across the flow (along the gradient, unit steps,
bilinear samples), then integrated along the tangent streamline with a
Gaussian of scale \`sigma_m\`, which joins edge fragments into coherent
strokes (Kang, Lee and Chui 2007). Streamlines take unit midpoint steps
for up to \`ceil(2 \* sigma_m)\` steps each way and stop early once they
have turned more than 20 degrees from their starting direction, so
corners are not rounded off. \`sigma_m = 0\` skips the integration. 4.
\*\*Threshold, then this package's placement and cleanup.\*\* The XDoG
soft threshold (Winnemoeller, Kyprianidis and Olsen 2012), \`T = 1\`
where the response \`D \>= eps\` and \`1 + tanh(phi \* (D - eps))\`
elsewhere, marks an ink band where \`T \< 1/2\`; with \`phi = Inf\` (the
default) that is \`D \< eps\`, and a finite \`phi\` moves the threshold
to \`eps - atanh(1/2) / phi\`. The band lies on the dark side of an
edge, about \`1.2 \* sigma\` deep, so its skeleton would sit beside the
edge. The remaining steps are this package's own, not part of the
published methods. Lines are placed between the band and its brighter
side: for each band pixel and each 4-neighbour outside the band that is
brighter (after the across-flow \`G\_(k sigma)\` blur, integrated along
the flow), the mask keeps whichever of the two lies nearer the zero
crossing of the pure difference of Gaussians (\`tau = 1\`, whose zero
crossing lies on an isolated step edge), or nearer the threshold
crossing where that does not change sign between them. Thin strokes are
drawn once, along their middle: sampling the response bilinearly across
the flow, a band at most \`3 \* sigma\` pixels wide with brighter pixels
on both sides is a dark stroke, and a gap of the pure difference of
Gaussians at most \`3 \* sigma\` wide with darker band on both sides
(not itself a dark stroke) is a light stroke. Each is drawn as its
centre line (if at least \`min_component\` pixels long), and edge lines
beside it are dropped. Parallel dark strokes 2 to 4 pixels apart
(hatching) merge into one band and are drawn as that band's two outer
outlines, about 2 pixels outside the outer strokes, with no line between
them; from about 6 pixels apart each gets its own centre line. Enclosed
holes of one or two pixels are filled, and eight-connected components
smaller than \`min_component\` pixels are removed.

The decision is a sign test of a filter of lightness, so with \`eps =
0\` and \`phi = Inf\` it compares edge contrast with the local
lightness: with \`tau \< 1\`, an edge is drawn when its lightness step
exceeds roughly \`(1 - tau) / 0.1\` of the darker side's lightness (20
percent at the default \`tau\`). Multiplying the lightness by a positive
constant therefore leaves the mask unchanged when the pre-filter is off
(\`smooth = 0\`), but the method is not invariant to scaling image
intensities: the sRGB to \`L\*\` conversion is non-linear, and the
Kuwahara weights use an absolute variance scale. Raising \`eps\` widens
the band; once it covers the whole image there is no band boundary left
and the mask is empty (XDoG's stylised output would be all ink).

FDoG can leave a gap at an occlusion T junction in the thresholded band.
A positive \`junction_gap\` enables a repair after mask cleanup. It fits
straight segments to the thinned contours and considers extending a free
endpoint to a nearby contour that crosses its projected path. It adds a
bridge when contrast across that path in the original CIE lightness
image has the same sign and at least 40 percent of the magnitude of the
median contrast along the visible stem. The stem must have at least 3.5
L\* units of median contrast. The new pixels must touch two distinct
patches of existing ink, one at each endpoint. The thresholded band
stays as computed by FDoG. Validate junction measurements against
annotations on representative images. For line art whose junctions
follow stroke centre lines, pass a centre-line mask to
[`trace_contours`](https://bbuchsbaum.github.io/mlvtoolboxR/reference/trace_contours.md).

The extraction is deterministic. With \`smooth = 0\` and \`junction_gap
= 0\`, the arithmetic is arranged so that quarter turns and reflections
of the image transform every stage, and so the mask, exactly. The
Kuwahara filter accumulates its sums in an order that these transforms
change, so with the pre-filter on the intermediate images agree only to
rounding (about 1e-14); the masks have agreed exactly in every test, but
a pixel within rounding of the threshold could differ.

With \`colour = TRUE\` the pre-filter and both structure tensors use all
three CIE \`L\*a\*b\*\` channels (each divided by 100), so colour
boundaries steer the flow and the smoothing; the difference of Gaussians
is still taken on lightness, so a boundary between two colours of equal
lightness is not drawn.

The defaults were chosen by a visual review of four public-domain
paintings 640 pixels wide (see the [Line
extraction](https://bbuchsbaum.github.io/mlvtoolboxR/articles/line-extraction.html)
article). \`sigma\`, \`sigma_c\`, \`sigma_m\` and \`smooth\` are in
pixels, so images at other resolutions want them scaled with the image.

## References

Kang, H., Lee, S. and Chui, C. K. (2007). Coherent line drawing.
\*Proceedings of the 5th International Symposium on Non-Photorealistic
Animation and Rendering (NPAR)\*, 43-50.
[doi:10.1145/1274871.1274878](https://doi.org/10.1145/1274871.1274878)

Kyprianidis, J. E. and Doellner, J. (2008). Image abstraction by
structure adaptive filtering. \*Proceedings of EG UK Theory and Practice
of Computer Graphics\*, 51-58.

Kyprianidis, J. E., Kang, H. and Doellner, J. (2009). Image and video
abstraction by anisotropic Kuwahara filtering. \*Computer Graphics
Forum\*, 28(7), 1955-1963.
[doi:10.1111/j.1467-8659.2009.01574.x](https://doi.org/10.1111/j.1467-8659.2009.01574.x)

Kyprianidis, J. E. (2010). Anisotropic Kuwahara filtering with
polynomial weighting functions. \*Proceedings of EG UK Theory and
Practice of Computer Graphics\*, 25-30.

Winnemoeller, H., Kyprianidis, J. E. and Olsen, S. C. (2012). XDoG: an
eXtended difference-of-Gaussians compendium including advanced image
stylization. \*Computers & Graphics\*, 36(6), 740-753.
[doi:10.1016/j.cag.2012.03.004](https://doi.org/10.1016/j.cag.2012.03.004)

## See also

\[line_extraction_stages()\] for the intermediate images;
\[trace_contours()\] and \[trace_rgb()\] to vectorise the mask.

## Examples

``` r
# A dark disc on light paper.
xy <- expand.grid(y = 1:48, x = 1:48)
img <- matrix(ifelse((xy$x - 24.5)^2 + (xy$y - 24.5)^2 < 15^2, 0.2, 0.9), 48, 48)
ink <- extract_lines(img)
sum(ink)
#> [1] 84
attr(ink, "extraction")$method
#> [1] "fdog-xdog-1"
d <- trace_contours(ink, simplify = 1)
length(d$contours)
#> [1] 1
```
