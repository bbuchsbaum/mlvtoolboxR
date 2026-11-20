# Line extraction as a trace_rgb() backend

\`xdog_lines()\` is \[extract_lines()\] with \`method = "fdog"\`, shaped
as a \[trace_rgb()\] backend: \`trace_rgb(image, backend = xdog_lines,
...)\` records the provider \`"xdog_lines"\`, its version
\`"fdog-xdog-1"\` and every parameter in the drawing's \`source\`.

## Usage

``` r
xdog_lines(image, ...)
```

## Arguments

- image:

  Greyscale matrix, RGB(A) array or image file name; see
  \[extract_lines()\].

- ...:

  Parameters of \[extract_lines()\], each named in full; \`method\`, if
  given, must be \`"fdog"\`.

## Value

The logical ink mask of \[extract_lines()\].

## Examples

``` r
rgb <- array(0.9, c(40, 40, 3)); rgb[10:30, 12:28, ] <- 0.25
d <- trace_rgb(rgb, backend = xdog_lines, min_component = 10)
d$source$provider_version
#> [1] "fdog-xdog-1"
```
