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
