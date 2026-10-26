# Reliable measurements and controlled stimulus construction

Mote epic `mlv-reliable-stimuli` has all fourteen child outcomes completed with
separate qualification gates. The implementation keeps `mlv2023` defaults for
compatibility and introduces explicit new measurement and construction APIs.

| Delivered capability | Entry points and evidence |
|---|---|
| Training/measurement compatibility | `feature_schema()`, `bind_feature_model()`, checked `predict_contours()`; missing/mismatched metadata fails before model invocation |
| Reusable units and methods | `feature_spec()`, pixels-per-degree calibration, estimator versions, fixed bins, metadata/export sidecars; [contract](FEATURE_SPEC.md) |
| Corrected optional estimators | `contour_turning()`, `curvature_profile()`, physical `score_branch()` windows; [frozen acceptance](ESTIMATOR_ACCEPTANCE.md), [Release A](QUALIFICATION_0.2.0.md) |
| Matched whole-contour subsets | `design_subsets()`, set-level/per-drawing constraints, exact bounded search, explicit feasibility/budget status; [guide](SUBSET_DESIGNER.md), [Release B](QUALIFICATION_DESIGNER.md) |
| Spatial relationships | `contour_graph()`, `graph_features()`, source spans, connectivity/continuation/parallel/enclosure relations; [guide](CONTOUR_GRAPHS.md), [qualification](QUALIFICATION_GRAPHS.md) |
| Backend research decision | [Boost segment-Voronoi spike](VORONOI_SPIKE.md): narrow-gap and correspondence evidence supports a separate implementation study; no production backend/default/dependency change |

Release A was qualified before the designer was implemented. Designer fixtures,
objectives, tolerances and compute budgets were frozen before optimization code.
Both required synthetic comparisons achieved target separation 270 versus ranked
100 at identical constraints; feasible controls tied. Infeasible inputs and
budget exhaustion were retained. Graph fixtures were likewise frozen before
implementation and match marginal histograms exactly while differing in the
intended spatial organization.

The new estimators and graphs have explicit geometric definitions. They do not
establish perceptual validity, border ownership, global symmetry, or universal
natural-image superiority. Raster extraction and mapped-pixel MAT histograms
remain resolution-sensitive. The bounded designer is not a guarantee that large
combinatorial experiments will finish. Curved/unbounded Voronoi edges remain
outside the research probe.

## Final integrated candidate

The final candidate includes the concurrent fused-statistics work and its last
large-index offset safety correction. The 394 package files match
[the exact source manifest](evidence/final/package-source-manifest.json), with
zero mismatches. The later local source checkpoint is bound in
[checkpoint evidence](evidence/final/checkpoint.json). Earlier phase receipts
retain their own source identities and are not silently relabeled as this tree.

| Check | Result |
|---|---|
| Installed tests | 1,332 assertions, zero failures/warnings/skips |
| R CMD check --as-cran --no-manual | Zero errors/warnings; two triaged notes |
| R and C++ line coverage | 95.0669% |
| Examples and both vignette rebuilds | Passed |
| Published numeric parity and frozen feature/designer/graph gates | Passed |
| Final package-source verification | All 394 files match |

The notes are “New submission” and inability to verify current time. The
installed testthat binary's R 4.5.2 build notice is distinct from test warnings
under the local R 4.5.1 runtime. Build/install/check took about 9.2/6.6/134.5
seconds and coverage 18.5 seconds. Sampled owned-workspace peaks were 14,920,658
bytes during check and 19,467,117 during coverage, below the existing 200 MiB cap.
These workspace measurements are not peak process memory.

[Final source archive](../mlvtoolbox_0.2.1-final.tar.gz): 352,747 bytes, SHA256
`7be461b5ae8d2447a6ff37b22dc8d88114fc4e5c9b25d9f084130c7fb52a12fe`.
It is ignored by Git and retained locally with the earlier candidate archives.
[Qualification receipt](evidence/final/qualification.json),
[full check](evidence/final/R-CMD-check.log),
[installed tests](evidence/final/testthat.Rout),
[coverage](evidence/final/coverage.json),
[environment](evidence/final/sessionInfo.txt).

The package qualification ran on macOS arm64. Actual Linux/Windows, MATLAB,
hosted CI, PDF manual generation and natural-image predictive validation remain
unexecuted. Local Git checkpoints preserve the work; no push, tag, package
submission or deployment is part of this implementation.

## Reproduction

Use `python3 tools/bounded_check.py --work-name verification` for an isolated,
bounded current-source build/check. It includes newly added non-ignored package
files and records content hashes. The frozen designer comparator/oracle are in
`tools/designer/`; graph generators are in `inst/fixtures/graph-v1.R`. Prototype
source, preprocessing and fixed fixtures are under `tools/voronoi/`, with native
input files and raw output in its evidence directory. Re-run experiments in an
isolated checkout, initialize a fresh explicit timebox, and build the package
library before invoking the prototype driver. Historical receipts are retained.
