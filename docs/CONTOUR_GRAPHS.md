# Contour relationships beyond marginal histograms

```r
graph <- contour_graph(drawing, spec = spec, tolerance = .001,
  gap_max = 3, angle_tolerance = 5, parallel_distance = 6,
  min_shared_extent = 1)
graph_features(graph)
graph$nodes
as.data.frame(graph)        # source contour/segment spans
graph$continuation
graph$parallel
graph$enclosure
```

The graph explicitly splits crossings, including self-crossings. It keeps bends
as degree-two nodes and records each source segment's parameter span. Node
clustering is transitive within the declared tolerance. Its representative
coordinate is a centroid; edge lengths retain the source geometry. Collinear
overlaps are rejected because their multiplicity needs a separate convention.
The existing published junction detector is unchanged.

`graph_features()` reports connected components, cycle rank, free endpoints,
branch nodes, source-contour closure, continuation pairs, parallel contour pairs,
shared projected extent and strict enclosure. Cycle rank is a graph invariant,
not a count of filled objects. Extra degree-two vertices change raw vertex/edge
counts but not cycle rank. A short span collapsed by tolerance is reported.

Good continuation requires two free source endpoints whose outward tangents
face across the gap. Parallel relations are span-level and report angle, shared
extent and mean normal separation. Their summed extent is pairwise support,
not union coverage. Enclosure is restricted to disjoint simple closed polygons
whose boundaries are separated by more than the node tolerance. Self-crossing,
touching and overlapping boundaries do not receive a strict enclosure relation.

Distance parameters and graph coordinates use pixels or degrees of visual
angle according to the feature specification. Angles remain geometric degrees.
`pixels_per_degree` can override input calibration. No canvas-based calibration
is inferred. `graph_features(graph, x = drawing, spec = spec)` verifies that
the requested source geometry and recipe still match the cached graph. Recompute
after a geometry or parameter change; serialization preserves the original
graph and its provenance, not a live link to a mutable drawing.

The initial implementation accepts at most 1,000 source segments and 20,000
candidate split points. It is intended for bounded contour drawings; dense
intersection arrangements can exceed the latter limit. Mean parallel separation
is NA when no qualifying relation exists.

The [frozen tests](GRAPH_ACCEPTANCE.md) compare four pairs with exactly matched
orientation/length/curvature histograms but different spatial organization. They
also establish topology by hand enumeration, test perturbation thresholds and
transformation invariances, and reject stale geometry/specifications. These
checks establish the declared geometric descriptors. They do not establish
border ownership, global symmetry, natural-image predictive value or perceptual
validity.
