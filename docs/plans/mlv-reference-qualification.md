# Plan: qualify mlvtoolbox against real MATLAB MLV outputs

Status: completed 2026-11-20 (approved 2026-11-19, with the instruction to strengthen the package and to adopt the analytically supported behaviour wherever it differs from MLV). Results: docs/PARITY.md and docs/evidence/mlvref/.

## Principle

Agreement with MLV is evidence, not the goal. Every divergence is adjudicated
against an independent oracle (analytical construction, exhaustive search,
hand-enumerated geometry or an invariance argument):

- Oracle disagrees with R → fix R, add a regression test, record as MLV-R00x.
- Oracle disagrees with MLV → keep the R behavior as the default, add a
  minimal synthetic reproducer with the analytical answer, record the
  difference as MLV-M0xx.
- Underdetermined (ties, unspecified smoothing, rasterizer) → choose a
  principled definition, document it, and quantify its effect on natural data.

Never loosen a tolerance or reclassify a difference to make a comparison pass.

## Reference sources

1. **Stored MATLAB outputs** in `reference-data/mlv` (local, not redistributable):
   30 dataset files, 3,017 drawings with MATLAB orientations, lengths,
   curvatures and junctions; 2,757 also carry per-contour parallelism,
   separation and mirror scores, `_allX/_allY/_allScores` samples and
   histograms. `examples/Forest_DBW.mat` also holds MATLAB `MAT.skeleton`,
   `distance_map` and `AOF` (768×1024). The producing MLV revision is unknown.
2. **Pinned MATLAB source under Octave 11.1 + image 2.20.0** (commit
   f12418a…): stage-wise execution, and a check that the stored outputs
   correspond to the pinned code (version drift).
3. **pyMLV** only as a secondary diagnostic.

Natural data never enter the package or Git. Committed artifacts are harness
code, small synthetic fixtures, and aggregate receipts (counts, quantiles,
hashes of inputs).

## Phases

### P0 Export (shared infrastructure)
`tools/mlvref/export.m`: Octave reads each dataset and writes, per dataset,
gzip CSVs under `.agent-work/mlvref/`: segments (drawing, contour, segment,
x1..y2, orientation, length, curvature), junctions (drawing, junction, x, y,
type, contour/segment IDs, angles), MAT scores (drawing, contour, segment,
parallelism/separation/mirror) and `_all` samples, plus a manifest with input
SHA-256s. Check free disk before export (target < 1 GB).

### P1 Version drift
On a stratified sample (each dataset), run pinned code under Octave on the
stored contours and compare with the stored fields. Classify each field as
"stored = pinned" or "drift". Only fields without drift support claims about
the pinned code; drifted fields are still MATLAB ground truth for "published
MLV values".

### P2 Geometry and junctions (track A)
- R vs stored MATLAB on all 3,017 drawings: segment geometry; junction
  matching (position within 1e-6, type, contour sets, angles).
- Synthetic junction suite (L, T, X, Y, arrow, ψ, near-miss, collinear
  continuation, endpoint-on-interior, multi-way, negative coords, merge chains)
  with analytical expectations; run R and Octave.
- Adjudicate every mismatch class; fix R or record MLV differences.

### P3 Medial axis (track B)
- Forest_DBW: recover MATLAB's binary raster from `distance_map` (ink ⇔ d ≤ 0,
  verify), feed it to R → compare distance, skeleton, AOF with MATLAB's;
  separately compare R rasterization with MATLAB's raster.
- Tie policy: quantify default-label effects on Forest and on stored
  per-contour scores. Evaluate a tie-invariant AOF (average boundary direction
  over all equidistant nearest pixels) against the lowest-index rule; adopt
  the principled one if it is defensible, keep `mlv2023`-compatible option only
  if it reproduces published values.
- End-to-end per-contour parallelism/separation/mirror vs stored MATLAB on all
  2,757 drawings; attribute residuals to rasterization, ties, smoothing.
- Adjudicate and fix as in P2.

### P4 Integration
Update PARITY.md, DIFFERENCES.md, IMPLEMENTATION_STATUS.md, tests and NEWS;
fresh-context review subagent; full test suite; CI green on all platforms.

## Constraints

- Disk: < 10 GB free at start. No R CMD check --as-cran locally; use CI.
- Octave lacks `insertShape` and MATLAB `smoothdata` heuristics; stage-wise
  inputs replace them, and that substitution is recorded.
- Commit per milestone; push only green states.
