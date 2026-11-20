# Intermediate images of the line extractor

Runs \[extract_lines()\] and returns every stage, for figures and
teaching.

## Usage

``` r
line_extraction_stages(image, ...)
```

## Arguments

- image:

  Greyscale matrix or RGB(A) array with values in \\0, 1\\
  (sRGB-encoded, as returned by \[read_mlv_image()\]), or a PNG/JPEG
  file name. Transparent pixels are composited over white.

- ...:

  Parameters of \[extract_lines()\], each named in full; \`method\`, if
  given, must be \`"fdog"\`.

## Value

A list of \`h x w\` matrices unless noted: \`lightness\` (CIE \`L\*\`
divided by 100), \`smoothed\` (after the Kuwahara pre-filter), \`flow\`
(a list of the unit edge-tangent components \`x\` and \`y\`, zero where
the tensor is isotropic, and the \`anisotropy\` in \\0, 1\\), \`dog\`
(the across-flow difference of Gaussians), \`response\` (after
integration along the flow), \`xdog\` (the XDoG soft threshold in \\0,
1\\), \`band\` (the logical ink band), \`mask\` (the result of
\[extract_lines()\]) and \`extraction\` (its definition and parameters).

## Examples

``` r
img <- matrix(0.9, 32, 32); img[8:24, 8:24] <- 0.2
s <- line_extraction_stages(img, min_component = 5)
names(s)
#> [1] "lightness"  "smoothed"   "flow"       "dog"        "response"  
#> [6] "xdog"       "band"       "mask"       "extraction"
identical(s$mask, extract_lines(img, min_component = 5))
#> [1] TRUE
```
