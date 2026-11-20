# Implementation status

The R package is version 0.3.0. This file supersedes the original
assessment's statements that no implementation exists.

## Current qualification (0.3.0)

- Qualified against outputs the original lab computed in MATLAB: 3,017 drawings
  with geometry and junctions, 2,757 with per-pixel medial-axis scores, and the
  complete Forest_DBW MAT. Pinned source under Octave reproduces the stored
  fields (no version drift). See [PARITY](docs/PARITY.md).
- Published-compatible options reproduce MATLAB: all segment geometry within
  1e-9; 395,125/395,134 junctions exactly; Forest distance map, flux and
  skeleton exactly; 99.1–99.4% of stored per-pixel symmetry scores within 1e-4,
  with the residual traced to MLV's single-precision arithmetic.
- Defaults use alternative definitions for 14 newly documented differences from
  MLV (MLV-M018–M031) and are
  tested for invariance to contour/branch direction and numbering and to lattice
  symmetries. Thirteen R-port revisions were made (MLV-R005–R017), four of them by an independent review.
- The installed suite passes locally (macOS arm64, R 4.5.1); hosted CI covers
  macOS, Windows and Ubuntu (R-devel, release, oldrel-1).
- 77 retained upstream MATLAB files have no explicit license grant and are not
  in the public repository ([licensing audit](vendor/LICENSING.md)).
- Not qualified: actual MATLAB execution (stored outputs are the reference),
  photograph-to-drawing extraction, and MATLAB rasterization ties on
  non-integer endpoints.

The sections below record the 0.1.1 port qualification and remain the basis
for parity claims.

## Implemented

- Drawing validation and stable contour IDs; subpixel geometry; explicit cleaning.
- Orientation, length, curvature, axial mean orientation and typed statistics.
- Exact sweep-filtered junctions, consistent angles and transitive merging.
- Conservative whole-contour splitting, model/histogram scores, rotation,
  shifting, stable circle clipping, junction segmentation and arc-length regions.
- Deterministic rasterization, exact EDT/nearest indices, compiled AOF, attributed
  Octave-image Pratt morphology, graph tracing, branch scores and tangent mapping.
- Valid-zero score masks, local score rendering, feature densities and batch tables.
- SVG paths/arcs/transforms/aspect ratios, v5 MAT structures, PNG/JPEG adapters,
  and explicit RGB extraction providers.
- All 62 original top-level functions mapped in [FUNCTION_MAP](docs/FUNCTION_MAP.md).
  Bundled SAM/legacy-forest models are not part of the package.

## Executed numerical evidence

- The installed suite has 583 passing assertions, including independent geometry,
  exhaustive EDT/nearest-search, metamorphic and public workflow tests.
- Eighteen synthetic cases execute pinned Octave/Python references across
  contours, histogram bins, rotations, MAT masks, AOF radii and branch shapes.
  All 289 comparisons are retained: 226 match and 63 have documented explanations.
- A retained forest drawing (562 contours / 2,072 segments after removing one
  zero-length contour) matches reference geometry and histograms within explicit
  tolerances. Python length histograms require a labeled shape adaptation.
- Six MAT distance maps agree. Default skeletons differ on four masks at exact
  nearest-boundary ties. Keeping R distances and substituting Octave labels
  makes all six skeletons identical, with AOF differences below 1.5e-14.
- Analytical boundary and distance tests led to two R-port revisions:
  odd-bin axial histogram boundaries and isolated-pixel signed distances.
- Exact old-output witnesses cover grouped statistics and tangent samples;
  filtered mapping agrees with exhaustive nearest-sample search.

See [PARITY](docs/PARITY.md) for the stage matrix, plots, tolerances, adapters,
natural drawing and the reference-implementation differences it reproduces. Actual MATLAB was not executed.

## Qualification (0.1.1)

Source d83ed952fdc9f639b98e849b0d8697259894de46 passes R CMD check
--as-cran --no-manual with 0 errors, 0 warnings and one informational
New submission NOTE. R and instrumented C++ coverage is 92.8989%
(previously 92.6117%). Examples and vignette rebuilds pass.

Paired local measurements improve the 1,000-contour feature table from 1.609 s
to 30.875 ms (52.1x), the long statistics table 45.6x, and full MAT on four
256x256 elliptical rings 1.91-1.92x. All five representative old/new workload
outputs are identical. Forest histograms take 6.656 ms, versus 619 ms in R 0.1.0,
83.292 ms in Octave and 7.471 ms in shape-adapted Python. These are local
workload-specific observations, not MATLAB JIT or complete image-pipeline claims.

See [qualification](docs/QUALIFICATION_0.1.1.md),
[performance](docs/PERFORMANCE_0.1.1.md) and receipts in
docs/evidence/parity-perf. The local installable archive is
mlvtoolbox_0.1.1.tar.gz; it is ignored by Git and can be rebuilt with
tools/bounded_check.py. The 0.1.0 archive and historical evidence are preserved.

Prepared CI covers macOS, Windows and Linux. No CRAN/R-universe submission,
site deployment or upstream issue submission has occurred.

The full Python image pipeline and actual MATLAB have not run. Natural-image
scientific acceptance, peak RSS, parallel scaling, general SVG paint effects,
MAT v7.3/HDF5 and bundled SAM extraction remain outside this qualification.
Rasterization, graph tracing, EDT ties, taper smoothing and axial mean orientation
have explicit method choices documented in DIFFERENCES.md.

## Storage and provenance

Cleanup verified 159 source files and 39 retained data files, removed about
2.63 GiB of logical content, and preserved one local reference-data collection.
See CLEANUP_RECEIPT.json, vendor/SOURCES.json and tools/restore_upstream.py.

The R implementation and attributed derived kernels are GPL-3.0-or-later.
Original reference notices and mixed rights remain recorded; this local port
does not assert dataset redistribution rights or blanket upstream relicensing.
