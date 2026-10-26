# Whole-contour subset design

`design_subsets()` assigns each eligible contour to A, B, or unassigned. The
primary mode matches an experiment across its complete drawing set. It returns
both drawing sets, original source/contour IDs, retained input geometry and
provenance, an explicit feature recipe, every residual, and a solver status.

```r
spec <- feature_spec(c("length", "orientation", "curvature"),
  normalize = FALSE, estimator = "arclength", curvature_scale = 0.5)
result <- design_subsets(drawings, spec,
  target = "integrated_curvature", length_bounds = c(40, 50),
  count_bounds = c(2, 6), nuisance = c(length = 2, orientation = 10),
  count_imbalance = 1, time_limit = 5, max_nodes = 100000)
result$status
result$diagnostics
as.data.frame(result)
saveRDS(result, "design.rds")
```

`result$conditions$A` and `$B` contain one drawing per original source, including
empty drawings. `unassigned` preserves remaining contours. All selection is
whole-contour and disjoint. Subset first to exclude ineligible contours. No
transformations occur. Geometry changes invalidate derived properties.

The target is an additive difference of integrated curvature, length, or
length-integrated horizontal-minus-vertical orientation. Raw histogram mass
constraints have a linear formulation and preserve their physical units. They
do not compare separate normalized distributions. Use the same fixed bin
specification across inputs, with `normalize = FALSE`. For physical inputs,
`pixels_per_degree` can supply a different calibration per drawing. The
curvature target is scale-dependent; opposite turns can cancel within support.

`mode = "per_drawing"` applies all bounds separately within each drawing. A
drawing with one contour cannot form two nonempty conditions. In pooled mode,
per-drawing diagnostic rows still show all imbalances; their local feasibility
flags do not decide pooled feasibility. No tolerances are silently relaxed.

The dependency-free exact search accepts at most 64 contours. Its node and
wall-time limits bound search effort, while preprocessing and independent
reconstruction checks have separate timings. Difficulty grows combinatorially;
the input cap is not a promise that 64-contour problems finish. Deterministic
traversal resolves tied solutions. Numeric feasibility tolerance is 1e-8 in the
declared constraint units. Objective certification uses the same tolerance.

| Status | Meaning |
|---|---|
| optimal | Search completed with a verified best solution, within numeric tolerance |
| infeasible | Search completed without any feasible assignment |
| feasible | Budget stopped search; a verified solution and conservative upper bound are returned |
| budget_exhausted | Budget stopped search before any feasible assignment was found |

For the last two states, the root upper bound may be loose. An absent solution
does not prove infeasibility. No-solution results have `conditions = NULL` and
unassigned condition labels (`NA`), rather than a fabricated empty design.
`solver$working_data_bytes` is an object-size estimate, not peak RSS.

Junctions and medial-axis scores change when contours are removed. Request
`diagnostics = "junctions"` or named MAT properties to recompute them on each
selected drawing. Results are under `relational_diagnostics`; they never become
optimization constraints. Physical MAT diagnostics require an arc-length branch
specification. Retaining an original relational score would be misleading.

The [frozen benchmark](DESIGNER_ACCEPTANCE.md) contains deliberately constructed
small combinatorial cases, not natural-image or perceptual validation. Raw
evaluation, budget exhaustion, infeasible cases and ties are retained together.
