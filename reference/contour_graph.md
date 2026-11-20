# Construct a contour arrangement and relationship graph

Exact crossings, including self-crossings, split source spans. Nearby
endpoints/intersections merge transitively; collinear overlaps error.
Bends remain degree-two nodes. Enclosure applies only to disjoint simple
closed polygons. Continuation links free endpoints; parallel relations
use shared projected extent and mean normal separation. These geometric
relations do not establish border ownership, global symmetry or
perceptual validity.

## Usage

``` r
contour_graph(
  x,
  spec = NULL,
  tolerance = 1e-06,
  gap_max = 3,
  angle_tolerance = 5,
  parallel_distance = 6,
  min_shared_extent = 1,
  pixels_per_degree = NULL
)
```

## Arguments

- x:

  A drawing.

- spec:

  Optional feature specification; uses the drawing's saved recipe or a
  default pixel recipe. Distances use its target spatial units.

- tolerance:

  Positive node-clustering tolerance and enclosure clearance.

- gap_max:

  Maximum continuation gap.

- angle_tolerance:

  Tangent/axial angular tolerance in geometric degrees, nonnegative and
  strictly less than 90.

- parallel_distance:

  Maximum mean normal separation of parallel spans.

- min_shared_extent:

  Minimum projected overlap for parallel spans.

- pixels_per_degree:

  Optional input calibration override.

## Value

An \`mlv_contour_graph\` with nodes, source-span edges, continuation,
parallel and enclosure tables, method specification and source geometry.
At most 1,000 source segments and 20,000 candidate split points are
accepted. Node coordinates are in target spatial units; edge lengths
retain source geometry even when clustering moves the representative
node centroid.
