# Thin a binary image

\`"minimal"\` (default) is \`"hmt8"\` thinning followed by sequential
removal (column-major scan, repeated to a fixed point) of every 8-simple
pixel that has at least two neighbours. No simple pixel with two or more
neighbours remains; a 2x2 block can remain only where every one of its
pixels is non-simple (removing any of them would change the 8-connected
topology). The result has the same 8-connected topology as the input,
keeps branch end points and is idempotent. \`"hmt8"\` is sequential
eight-direction hit-or-miss thinning with edge and corner structuring
elements, applied with background to the S, SW, E, SE, W, NW, NE and N
in turn. It preserves topology but is not always one pixel wide (a 4x4
square leaves a 2x2 block; 4-connected staircase corners remain). That
order reproduces MATLAB \`bwmorph(...,"skel",Inf)\` exactly on the
Forest_DBW reference (flux skeleton and ink thinning) and is the
published-compatible choice; MATLAB's own algorithm is undocumented, so
wider equivalence is not claimed. \`"pratt"\` reproduces the Octave
image 2.20.0 \`bwmorph(...,"skel",Inf)\` lookup tables; its skeletons
are not one pixel wide (strokes of even width keep two-pixel ribbons and
2x2 blocks). All three are sequential or direction-ordered, so none is
equivariant under lattice rotations/reflections for general inputs.

## Usage

``` r
skeletonize(mask, min_component = 1L, method = c("minimal", "hmt8", "pratt"))
```

## Arguments

- mask:

  Binary ink matrix.

- min_component:

  Minimum eight-connected component size.

- method:

  Thinning algorithm: \`"minimal"\` (default), \`"hmt8"\` or
  \`"pratt"\`.

## Value

Logical skeleton, with backend metadata described in the package docs.
