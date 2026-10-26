# Frozen contour graph acceptance v1

Frozen before graph implementation. The deterministic generator is
`inst/fixtures/graph-v1.R`; no random seed or learned provider is used.
Each paired comparison must have matching raw orientation, contour-length and
arc-length-curvature histograms under the shared fixed specification. No claim
is made about other marginals unless explicitly tested. Graph relationships must
distinguish the pairs without reading their marginal histogram differences.

## Meaning of the graph

Nodes are source endpoints, polyline bends, and exact segment intersections,
including self-crossings. Segments split at every intersection; edges retain
source contour/segment identity and span parameters. Collinear overlaps have
ambiguous multiplicity and are rejected explicitly; shared endpoints and
consecutive collinear segments are allowed. Legacy junction behavior is untouched.
Degree-two subdivision vertices are retained as provenance; they do not alter
connected components, cycle rank, endpoint count or branch-node count.

Points within the positive spatial tolerance are merged by transitive connected
components and represented by their centroid. This can collapse short spans;
report the count. The tolerance applies to existing endpoints/intersections,
not an implicit endpoint-to-segment snapping rule. Exact crossings split both
segments. All distances/lengths use the specification's target spatial unit;
calibration is explicit. Node positions describe tolerance clusters, while edge
lengths retain original source-span lengths.

Topological cycle rank is `edges - vertices + connected_components`, counting
multiedges and loops. It is not a count of independent objects or holes in a
filled image. A source contour is closed only when its final endpoint equals
its initial point. Geometric enclosure is limited to disjoint simple closed
polygon contours, separated by more than the node tolerance. Intersecting,
self-crossing or touching boundaries do not establish enclosure. T-junctions do
not establish border ownership.

Good-continuation relations connect two free degree-one source endpoints across
a gap: both outward tangents face the connecting line within the angular
tolerance, and gap length is within the maximum. They need not be from distinct
components. Closed source contours supply no free endpoints.

Parallel relations concern spans of distinct source contours, within the axial
angle tolerance, with at least the minimum common projected extent and bounded
mean normal separation. Intersecting spans are excluded. Shared extent is
measured on the normalized mean axis; separation is averaged at its overlap
endpoints. Report segment-level relations and unique contour-pair counts.
Summed pairwise shared extent is not union coverage. Straightness or local
parallelism does not prove reflection symmetry.

## Frozen cases and expectations

| Pair | Matched marginals | Required distinction |
|---|---|---|
| Four box sides vs separated sides | Orientation, length, curvature | One component/one cycle vs four components/no cycle |
| Collinear gap vs staggered pair | Same three marginals | One continuation vs none |
| Parallel lines 5 vs 20 pixels apart | Same three marginals | One parallel contour pair vs none at distance limit 6; shared extent 10 and separation 5 |
| Nested vs separated simple squares | Same three marginals | One strict enclosure vs none; both have two graph cycles |

The generator also fixes exact vertex/edge/component/cycle/degree counts for a
cross, a T, a self-crossing bowtie, empty input, and endpoint gaps 0.0005/0.002
under tolerance 0.001. Bowtie cycle rank is two and it is not a simple polygon.
Overlapping collinear segments must error. Topology expectations follow hand
enumeration in the generator, not production graph output.

Require invariance of relational/topological summaries under reversal, rotation,
translation and physical-coordinate scaling with calibration. Raw node/edge
counts may change with collinear subdivision; component/cycle/degree summaries
and total parallel extent must not. Check near-threshold and negative controls.
Serialization and tabular output preserve source identity; querying a graph
against changed geometry or a different spec must reject the stale result.

Numerical comparisons use absolute tolerance 1e-8 except spatial node matching,
which uses the declared 0.001. Bound the initial implementation to 1,000 source
segments and 20,000 candidate split points; exceeding these is an explicit input
limit, not an optimality or scientific result. Record runtime on a 10-by-10 line
grid and keep the small-fixture suite below 10 seconds on the local qualification
machine. No natural-image predictive claim follows from this gate.
