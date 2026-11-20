# Trace RGB images with an explicit extraction provider

Trace RGB images with an explicit extraction provider

## Usage

``` r
trace_rgb(image, backend, ..., min_points = 2L, spec = NULL)
```

## Arguments

- image:

  Numeric RGB array or image filename.

- backend:

  Function returning a logical ink mask from the image.

- ...:

  Arguments passed to the provider.

- min_points:

  Minimum traced path size.

- spec:

  Optional specification, including explicit provider name/version.

## Value

A drawing. Model weights and preprocessing belong to the provider.

## Details

A backend may attach an \`"extraction"\` attribute to its mask,
\`list(method = \<version string\>, parameters = \<list\>)\`, as
\[extract_lines()\] does. Without a specification the drawing then
records that version as \`source\$provider_version\` and the parameters
as \`source\$extraction_parameters\`. A specification whose extraction
settings carry line-extraction parameters (\`extraction\$lines\`, see
\[feature_spec()\]) passes them to the backend, which must not also
receive them through \`...\`, and the mask's attribute must then match
the specification's version and parameters exactly. With any
specification, a backend whose attribute names a version other than the
specification's \`provider_version\` is an error.
