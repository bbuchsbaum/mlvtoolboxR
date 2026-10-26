# mlvtoolbox 0.2.0

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
