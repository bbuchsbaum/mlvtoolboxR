# mlvtoolbox 0.3.0

## Qualified against stored MATLAB outputs

* The package was compared with outputs the original lab computed in MATLAB
  for 3,017 line drawings (2,757 with medial-axis scores). Published-compatible
  options now reproduce them: `detect_junctions(rule = "mlv2023")` matches
  395,125 of 395,134 junctions exactly, and `compute_mat(aof = "mlv2023",
  morphology = "hmt8", branch_rule = "mlv2023")` reproduces Forest_DBW's distance map, flux and
  skeleton exactly and 99.1-99.4% of all stored per-pixel symmetry scores within
  1e-4. See `docs/PARITY.md`.

## Behaviour changes (corrected defaults)

Defaults now correct demonstrated MLV defects (MLV-M018-M031 in
`UPSTREAM_BUGS.md`). Feature specifications record the new method definitions,
so specifications and cached MATs from 0.2.x are rejected and must be rebuilt.

* `detect_junctions()`: arms are the directions in which each contour leaves the
  junction neighbourhood (`rule = "exit"`), independent of contour direction,
  numbering and closed-contour start vertex; the absolute gap is Euclidean
  (`gap_metric = "euclidean"`), so junctions are rotation-invariant.
* `compute_mat()`: pixel-centred, tie-invariant average outward flux
  (`aof = "centered"`, which requires `samples` divisible by 4); a skeleton
  with no removable pixel left (`morphology = "minimal"`: MATLAB-exact
  hit-or-miss thinning, `"hmt8"`, followed by removal of every remaining simple
  pixel with at least two neighbours; 2x2 blocks remain only where no pixel is
  simple); branch tracing, scoring and mapping that do not depend on branch
  direction, numbering or loop start (`branch_rule = "invariant"`). Isolated
  closed loops are scored with circular windows.
  `skeletonize()` and `trace_contours()` use the same thinning.
* Skeleton pixels lying on ink are excluded before branch tracing.
* `mlv_drawing()` treats segments whose endpoints differ only by rounding residue
  as zero-length.

## Bug fixes

* Junction arm angles no longer depend on the merged junction position, and
  gapped T-junctions no longer gain a phantom arm (MLV-R005-R008, R014).
* Sweep and exhaustive junction search agree exactly, including gaps exactly at
  the extension limit, independent of segment direction (MLV-R015).
* Dropping rounding-residue segments keeps contours exactly connected and
  closed (MLV-R016).
* Width-1 rasterization is independent of segment direction (MLV-R010).
* Coincident branch points no longer produce NaN tangents (MLV-R011).

# mlvtoolbox 0.2.1

* Orientation, length and curvature histograms are computed in a fused native
  pass directly from contours. Measured `feature_table()` workloads are
  3.8-9.4x faster than 0.2.0; all 122 baseline public outputs and their
  metadata are unchanged (see `docs/FUSED_STATISTICS.md`).
* Added contour graph descriptors: source-linked topology, continuation,
  parallel relations and strict enclosure, with calibrated cache validation
  and rejection of stale geometry, parameters and calibration.
* Credited the original MLV, Python AOF and GNU Octave image contributors in
  `Authors@R`, added software and methods citations, and installed the full
  applicable license texts and original notices.

# mlvtoolbox 0.2.0

## Subsequent candidate: whole-contour stimulus design

* Added a bounded exact subset designer with set-level matching, optional
  per-drawing constraints, raw histogram L1 tolerances and explicit resource
  limits. Optimal, feasible, infeasible and exhausted outcomes are distinct.
* Results retain original contour identities, input geometry/calibration,
  reconstructed conditions, independently recomputed constraints and optional
  fresh junction/MAT diagnostics. A complete stimulus-design vignette is included.

## Feature contracts and estimators

* Prediction now requires the actual training feature schema and rejects bin,
  predictor, unit, and method mismatches before invoking a model. An explicit
  warning-emitting legacy bypass remains available.
* Added versioned feature specifications with fixed bins, pixels-per-degree
  calibration, physical curvature scales and branch windows, and extraction
  provenance. Batch statistics and CSV sidecars preserve measurement metadata.
* Added vertex-turning and arc-length tent curvature estimators, closed-contour
  wraparound, and exact integration of affine curvature into histograms.
* Added arc-length branch scores with independent analytical and invariance
  tests. These are explicitly versioned alternative geometric definitions.
* Published `mlv2023` contour and branch behavior remains the default for one
  release. Existing numeric parity fixtures remain unchanged.
* Geometry changes invalidate derived measurements; incompatible MAT caches
  are rejected. Mapping coverage distinguishes missing scores from valid zero.
