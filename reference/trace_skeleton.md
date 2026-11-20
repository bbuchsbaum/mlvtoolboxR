# Trace ordered medial-axis branches

Skeleton pixels with a nonpositive signed distance lie on ink and are
not medial points; they are excluded before tracing on every rule, and
their count is returned as attribute \`skeleton_ink_excluded\`.
\`rule="invariant"\` decomposes the pixel graph (diagonal steps only
where no orthogonal connection exists) into maximal paths between
end/junction pixels, then prunes spurs and merges to a fixed point. A
\*spur\* is a path from an end pixel to a junction; it is short when it
has fewer than \`min_points\` pixels (junction included). Each round
removes the shortest short spurs (all tied ones together) at every
junction that keeps at least two paths afterwards; a junction where
removing its tied shortest spurs would leave fewer than two paths is
protected and keeps all its paths. Paths are then joined through every
pixel where exactly two remain. Without ties this is sequential
shortest-first spur pruning; with ties (the four equal arms of a square)
every tied spur is kept, so the result does not depend on traversal
order, start pixel or lattice orientation. Only spurs are pruned: a path
between two junctions, or a component that is a single path, is kept
whatever its length, so no skeleton component loses its last path. A
closed path whose start pixel joins no other path is an isolated loop.
Junction pixels are shared by the paths that meet there but never
repeated within a path. \`rule="mlv2023"\` is the published MLV tracing:
Kovesi \`edgelink\` with \`cleanedgelist(minlength=min_points)\`
(Euclidean length), ported line by line, including its raster-scan
traversal direction and the duplicated point at every merged junction.

## Usage

``` r
trace_skeleton(mat, min_points = 10L, rule = NULL)
```

## Arguments

- mat:

  An mlv_mat object.

- min_points:

  Minimum spur size: pixels for \`"invariant"\` (spurs only), Euclidean
  length for \`"mlv2023"\` (cleanedgelist spurs and isolated edges).

- rule:

  Branch rule; defaults to the rule recorded in \`mat\`.

## Value

A list of branch data frames with row, column, x, y, radius and AOF.
