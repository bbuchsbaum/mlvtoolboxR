# Read MATLAB v5 vector drawing structures

Read MATLAB v5 vector drawing structures

## Usage

``` r
read_mlv_mat(file, variable = "vecLD", zero_length = c("error", "drop"))
```

## Arguments

- file:

  MAT file path.

- variable:

  Name of the drawing structure array.

- zero_length:

  Zero-length input policy.

## Value

One drawing for a scalar structure, otherwise a list of drawings.
