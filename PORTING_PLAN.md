# Current execution status

The cleanup and R implementation are complete, with a qualified 0.1.1
performance/parity follow-up. See
[IMPLEMENTATION_STATUS](https://bbuchsbaum.github.io/mlvtoolboxR/IMPLEMENTATION_STATUS.md),
[qualification](https://bbuchsbaum.github.io/mlvtoolboxR/docs/QUALIFICATION_0.1.1.md),
[performance](https://bbuchsbaum.github.io/mlvtoolboxR/docs/PERFORMANCE_0.1.1.md),
[parity](https://bbuchsbaum.github.io/mlvtoolboxR/docs/PARITY.md) and
the [differences
register](https://bbuchsbaum.github.io/mlvtoolboxR/DIFFERENCES.md).
Current evidence: 583 assertions, 92.8989% R/C++ coverage, check 0
errors / 0 warnings, 18 synthetic reference cases and a retained forest
drawing. Default MAT tie differences are explicitly retained. The
assessment below is the historical design baseline; its
no-implementation statements are superseded.

# MLV Toolbox to R: assessment and porting plan

Assessment date: 2026-10-25; revised after local Octave and Python
probes. Proposed package name: `mlvtoolbox`.

## Recommendation

Proceed with an analysis-first R package: vector drawings and contour
statistics, followed by a qualified medial-axis pipeline, then stimulus
manipulation and input adapters. Keep photograph-to-contour extraction
optional. This order gives a useful package early while testing the
hardest scientific dependency before spending time on peripheral
plotting and import features.

The intended outcome is a more robust, faster implementation that, by
default, uses an alternative definition wherever a demonstrated
difference from the stated method exists. Preserve the scientific
quantities and make all method changes explicit. Upstream outputs help
explain existing behavior; independent analytical and numerical evidence
decides the R behavior when an upstream implementation departs from that
evidence.

This is a feasible medium-sized scientific software project. Source
translation is the smaller part of the work. Correctness, robustness,
computational cost, and reproducibility are acceptance criteria for each
delivered capability. Record compatibility differences instead of
reproducing documented edge-case behavior to obtain matching output.
Claims of equivalence or speedup require executed evidence for the
stated inputs, source revisions, runtimes and backends.

This document is a plan, not an implementation or a numerical
qualification report. The local workspace was empty at assessment time.

## Evidence and scope

Upstream:
[bwlabToronto/MLV_toolbox](https://github.com/bwlabToronto/MLV_toolbox).
Pinned source:
[`f12418a5e03c2f38e222fed1dffccf42b227cd29`](https://github.com/bwlabToronto/MLV_toolbox/commit/f12418a5e03c2f38e222fed1dffccf42b227cd29),
the `main` head returned by GitHub during this assessment, committed
2025-04-30.

The complete Git tree was inspected, and all 62 top-level `MLVcode/*.m`
files were fetched. Detailed review focused on the main analysis
pipelines, statistics, drawing representation, junctions, rasterization,
tracing, SVG handling, and manipulation. Counts below exclude vendored
code unless stated:

| Component | Observed size | Porting implication |
|----|---:|----|
| Top-level MATLAB core | 62 files; 171,451 bytes; 5,045 lines | Small enough for a direct, selective port |
| Nonblank lines not beginning with a MATLAB comment | 2,510 | A source-size indicator, not a complexity estimate |
| Auxiliary contour-tracing directory | 8 MATLAB files | Assess reachability; it is not the principal tracing path |
| Vendored structured-edge stack | Most of the remaining `MLVcode` tree | Do not translate the entire dependency tree |
| Dataset directory | 36 files; 682,150,255 bytes, about 651 MiB | Keep outside the package and routine CI |
| First-party automated tests / CI | None identified in the inspected tree | Build characterization and independent tests first |

The published [2023 methods
paper](https://www.frontiersin.org/journals/computer-science/articles/10.3389/fcomp.2023.1140723/full)
provides scientific definitions and example workflows. Use it alongside
the pinned code: the paper predates the current SAM extraction path.

The user supplied local reference checkouts. MATLAB source was observed
at `MLV_toolbox/`, and the Python port at `vendor/pyMLV/`; refresh their
locations as the user finishes the vendoring layout. The MATLAB HEAD
matches the pinned SHA above. Python
[pyMLV](https://github.com/bwlabToronto/pyMLV) was pinned to
[`14b3183b91146b1021eb7ef97db75ad703b06832`](https://github.com/bwlabToronto/pyMLV/commit/14b3183b91146b1021eb7ef97db75ad703b06832).
Both tracked-file status checks were clean after the probes. Treat the
reference checkouts as read-only and put diagnostic adapters outside
them.

Small R arithmetic checks, upstream contour calls under Octave, one
adapted Octave MAT smoke, and a Python subpixel-orientation regression
were executed. See
[ASSESSMENT_EVIDENCE.md](https://bbuchsbaum.github.io/mlvtoolboxR/ASSESSMENT_EVIDENCE.md)
for inputs, outputs, source revisions, adaptations and qualification
boundaries. No R implementation, end-to-end cross-language parity,
speedup or image-quality result is claimed.

## Octave and Python reference strategy

Use local Octave actively for the runnable MATLAB source paths. Octave
11.1.0 and image 2.20.0 are installed. Loading the image package
supplies `bwdist` with nearest-feature indices, `bwmorph`, `bwareaopen`,
and `imbinarize`. The current environment lacks `smoothdata`,
`knnsearch`, `insertShape`, `xmlread`, and `imsegsam`. Record a
per-operation capability matrix instead of treating either complete
MATLAB compatibility or complete unavailability as the baseline.

The unmodified MAT path hits Octave’s parser at two `jj++m_Neighbors8`
expressions in `computeGradientVectorField.m`. A temporary copy changed
only these expressions to explicit addition. With that adaptation,
`octave-cli` ran a 25-by-35 synthetic MAT case, producing finite
distance/AOF matrices and 44 skeleton pixels. This is an executability
probe, not a correctness or performance result. Retain a SHA-bound
adapter manifest for future reference runs; never silently patch the
reference checkout.

Octave’s installed `bwmorph` explicitly documents differences from
MATLAB’s `skel`. Distinguish analytical correctness, Octave agreement,
Python agreement, and actual MATLAB agreement. MATLAB-specific evidence
remains useful where proprietary rasterization or morphology affects a
claim, but MATLAB access does not block building and scientifically
validating the R implementation.

The Python tree has 60 top-level code modules, with 52 same-named MATLAB
core counterparts and eight additional helpers. No package build
metadata, automated test directory or GitHub workflow was identified in
the pinned tree. Several MATLAB capabilities have no same-named Python
module, including SVG import, raster rendering and RGB/edge-map
extraction. Its README also reports that the medial-axis demo currently
needs MATLAB-rendered input.

Review useful NumPy/SciPy/OpenCV choices and compare selected runnable
functions. The Python source is a second implementation to examine, with
shared lineage; agreement between the ports is not independent
scientific validation.

Demonstrated and source-level Python concerns include:

- An executed subpixel segment test returns orientation 0 degrees rather
  than the analytical 315 degrees because endpoints are converted to
  int32.
- Branch scoring/tracing imports missing local modules: `smoothData`,
  `diff`, `fitLineSegments`, and `GetConSeg`.
- The gradient-field code loads a fixed `outerBoundaryOriginal.mat` and
  opens plots; numerical computation must not depend on external
  diagnostic state.
- Skeleton thresholding aliases and overwrites the flux array, and
  component cleanup replaces its calculated threshold with 100.
- Branch scoring contains derivative/window differences and an
  uninitialized taper intermediate. These source findings require
  targeted execution before relying on the corresponding outputs.

The Python root declares MIT, while numerical files retain GPL notices.
Include both repositories in the per-file rights inventory; the root
Python license does not settle the terms for all included components.

## Capability assessment

| Area | Main source anchors | Difficulty / recommended approach |
|----|----|----|
| Representation and contour measures | `computeOrientation`, `computeLength`, `computeCurvature` | Low: vectorized base R |
| Histograms and summaries | `get*Stats`, `averageProperty` | Low–medium: preserve bin boundaries, weighting, normalization, and empty cases explicitly |
| Junction detection and types | `detectJunctions`, `cleanupJunctions`, `computeJunctionAnglesTypes` | Medium: numerical geometry plus heuristic thresholds; start with a simple reference implementation |
| Rendering | `renderLinedrawing`, `MATpropertiesToContours` | Medium–high: rasterization changes downstream scientific results |
| Signed distance and average outward flux | `computeGradientVectorField`, `computeAOF` | High: exact Euclidean distance, nearest-boundary indices, boundary conditions, and nested loops |
| Skeletons and branches | `extract2DSkeletonFromBinaryImage`, `traceSkeleton`, Kovesi `edgelink` | High: thinning and graph topology are part of the method |
| MAT scores and contour projection | `computeMATpropertyPerBranch`, `mapMATtoContour` | High: smoothing, derivatives, tangent geometry, nearest-neighbor ties, and missing scores |
| Manipulation and statistical models | `splitLD*`, `rotateLinedrawing`, `predictContoursByStatsModel` | Medium: retain scientific selection rules; use R model prediction interfaces |
| SVG and raster drawing input | `importSVG`, `traceLinedrawingFromEdgeMap` | Medium–high: SVG support is broad but experimental; raster tracing shares skeleton dependencies |
| RGB input | `traceLineDrawingFromRGB` | Separate optional integration: SAM or structured-edge model plus postprocessing |

The main MAT code is MATLAB source, not an opaque native binary. The
compiled MEX dependency stack is concentrated in the legacy
structured-edge backend. That makes a native R/C++ analysis
implementation practical without porting the whole photograph extraction
stack.

## Scientific definitions and explicit changes

The definitions below establish the reference method. They do not
require preserving behavior demonstrated to depart from the stated
definition. Use the stated definition by default, document the changed
results, and test the invariant. A substantive new estimator, skeleton
method or approximation must have an explicit method identity and its
own validation.

1.  **Drawing geometry.** A contour is an ordered numeric matrix with
    columns `x1, y1, x2, y2`, not merely an unordered collection of
    edges. Preserve contour identity, segment order, and pen lifts.
    Dimensions upstream are `[width, height]`, whereas matrices index
    `[row, column]`.

2.  **Coordinates and orientation.** Image y increases downward. The
    upstream orientation is `atan2(y1-y2, x2-x1)` in degrees modulo 360.
    Orientation histograms use modulo 180, while curvature and junction
    calculations retain directed orientation. Skeleton fields called `X`
    and `Y` actually feed MATLAB row and column indexing; adapters must
    make this distinction explicit.

3.  **Curvature.** Preserve the piecewise definition: the wrapped
    absolute angle change divided by the current segment length, with
    the last segment compared to the previous one. A single-segment
    contour has curvature zero. This is not interchangeable with fitting
    a smooth curve and taking its differential curvature.

4.  **Statistics.** Orientation and curvature bins collect segment
    length; length bins collect whole-contour length. Junction type and
    angle histograms count events, with per-contour participation
    tracked separately. Several normalized statistics are scaled per
    10,000 pixels of contour length. Specify bin-edge and overflow
    behavior explicitly, characterize upstream differences, and test
    boundaries explicitly instead of inheriting the default settings of
    R’s [`hist()`](https://rdrr.io/r/graphics/hist.html).

5.  **MAT construction.** Preserve image polarity, the two distance
    transforms, nearest-boundary map, 60 angular samples, epsilon 1,
    flux thresholding, thinning, and component removal. A generic
    morphological skeleton is not automatically equivalent to this
    AOF-based transform.

6.  **MAT scores.** Preserve branch sampling, radius smoothing,
    derivative padding, window K=5, final moving mean of width 3, and
    the tenth-power score transformation. Taper is available but absent
    from the default set of three properties. Its unparameterized
    `smoothdata()` calls need explicit reference characterization.

7.  **Projection and missingness.** Preserve the maximum across
    contributing skeletal branches, tangent-point mapping and radius
    cutoff. Upstream converts zero-valued projected samples to missing
    values; determine where this means “unmapped” versus a valid zero
    score. The R representation should distinguish these cases.

Source anchors: [contour
measures](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/computeContourProperties.m),
[orientation](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/computeOrientation.m),
[curvature](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/computeCurvature.m),
[MAT
extraction](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/extract2DSkeletonFromBinaryImage.m),
[branch
scores](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/computeMATpropertyPerBranch.m),
[projection](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/mapMATtoContour.m).

## Issues to resolve before claiming compatibility

Horizontal/vertical histogram loss, the empty junction-angle dimensions,
and split-budget overshoot were reproduced by running the upstream
MATLAB functions in Octave. Other entries below remain source-level
findings or arithmetic checks. The Python subpixel-orientation
difference was executed separately.

| Finding | Evidence | Proposed treatment |
|----|----|----|
| Exactly horizontal segments disappear from the horizontal/vertical histogram | The Octave call on a length-10 horizontal segment returned total 0; the strict comparison excludes thisHV=1 | Include the upper endpoint and require histogram mass conservation |
| Empty junction-angle output has inconsistent default dimensions | The Octave call returned 5 histogram entries and 8 centers | Return a consistent zero vector of the requested length |
| Omitted MAT threshold differs from explicit 28 degrees | `computeMAT` uses 0.3 when omitted, versus `2/pi*sin(28*pi/180)=0.2988749` | Make omitted and explicit threshold_angle=28 use the same conversion; retain raw-threshold control for replication |
| A split’s top group can exceed its requested length fraction | The Octave length-ranked split of lengths 4,4,2 at fraction 0.5 retained length 8/10 at the top after sorting | Enforce the documented whole-contour length budget, with deterministic ties and explicit unassigned contours |
| Curvature ranking description and calculation differ | `splitLDbyProperties` describes an average but computes a length-weighted sum without dividing by contour length | Preserve the actual criterion under an explicit name; add a separately named mean if useful |
| Zero scores and missing projections are conflated | `MATpropertiesToContours` marks `thisProp~=0` as valid | Carry a separate mapping mask and test the effects of retaining valid zeros |
| RGB function help and actual default differ | Help mentions structured-edge detection; the arguments block sets SAM | Make the backend explicit and store its identity and settings |

The R definition is the default. For each demonstrated difference,
retain the minimal reproducing input, upstream behavior, independent
expected result and regression test. Maintain a concise table of
intentional differences.

Do not add a package-wide “legacy mode” prematurely. Add a narrow
replication option only when a real analysis requires the historical
behavior. Ambiguous scientific choices require a stated definition and
evidence, rather than silently changing the quantity being measured.

Direct source: [horizontal/vertical
bins](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/getHorizontalVerticalStats.m#L53),
[empty junction
output](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/getJunctionStats.m#L76),
[threshold
default](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/computeMAT.m#L31),
[splitting](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/splitLDbyProperties.m#L126),
[score
missingness](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/MATpropertiesToContours.m#L44),
[RGB
extraction](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/traceLineDrawingFromRGB.m).

## Proposed R design

Use a small S3 API. Do not export 62 one-for-one wrappers.

- `mlv_drawing`: ordered list of numeric segment matrices, width,
  height, stable contour IDs, original-image metadata, and provenance.
- Derived segment and contour properties: explicit tables attached to
  the drawing or returned as tables; geometry-changing operations
  invalidate or recompute them. Do not retain stale cached properties
  after rotation or segmentation.
- `mlv_mat`: signed-distance matrix, nearest-boundary indices, AOF
  matrix, logical skeleton, branch records, parameters, and provenance.
- Raster data: ordinary `[height, width]` matrices, with TRUE meaning
  ink in a logical drawing mask. Convert polarity at upstream-compatible
  boundaries. Raster-aligned pixel centers are (1,1) at the top left.
  SVG coordinate transforms and clipping rules must be explicit.
- Statistics: data frames carrying property, bin center/edges, value,
  units, normalization, and contour ID where relevant.

Proposed public verbs, not yet implemented:

``` r

drawing <- mlv_drawing(contours, width = 640, height = 480)
drawing <- compute_contour_properties(drawing)
stats <- contour_statistics(drawing)
plot(drawing, property = "curvature")

mask <- rasterize(drawing)
mat <- compute_mat(mask, threshold_angle = 28)
mat <- compute_mat_properties(mat, properties = c("mirror", "parallelism", "separation"))
drawing <- map_mat_properties(drawing, mat)

parts <- split_contours(drawing, property = "length", fraction = 0.5)
```

Provide
[`read_mlv_mat()`](https://bbuchsbaum.github.io/mlvtoolboxR/reference/read_mlv_mat.md),
[`read_mlv_svg()`](https://bbuchsbaum.github.io/mlvtoolboxR/reference/read_mlv_svg.md),
tabular export and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html). For
model-based selection, accept an R prediction function or supported
[`predict()`](https://rdrr.io/r/stats/predict.html) method with a
specified feature schema. Importing serialized MATLAB
regression/ensemble objects is outside the initial scope.

Use base R for geometry, statistics, object handling and initial plots.
Use Rcpp for measured hotspots and kernels requiring precise raster
semantics, especially AOF and potentially distance transforms/thinning.
Keep small, transparent R reference implementations for numerical
kernels where feasible.

Make `xml2`, `png`, `jpeg`, and MATLAB-file readers optional adapters.
Detect MAT-file formats before deciding how to read them. Evaluate an
existing image library against fixtures before accepting its distance or
morphology routines as a backend; availability of a similarly named
operation is insufficient.

Avoid imposing Python, deep-learning weights, a GPU, or a large spatial
stack on users who only need contour statistics. Optional RGB extraction
should accept precomputed masks/labels, or use a separately configured
local backend. A Canny alternative would be a different extraction
method, with its own name and provenance.

Suggested package layout: `R/` for objects and workflows, `src/` for
qualified kernels, `tests/testthat/` for behavior and oracle
comparisons, `inst/fixtures/` for tiny reusable inputs,
`tools/reference/` for Octave/MATLAB reference runners and adapters, and
`vignettes/` for reproducible workflows.

## Robustness acceptance criteria

- Validate dimensions, finite coordinates, column meanings and segment
  order at construction/import boundaries. Normalize MATLAB-file nesting
  once at import. Preserve subpixel doubles and use explicit index
  conversions at native/image boundaries.
- Define empty results with stable types and shapes. Reject or
  explicitly clean zero-length segments while retaining ID mappings.
  Undefined quantities must have documented missing/error behavior
  rather than arbitrary epsilon denominators that turn them into
  plausible scores.
- Exercise nearly parallel/collinear segments, duplicate points, tiny
  branches, zero radii, extreme coordinate magnitudes,
  all-ink/all-background masks and border-touching shapes. Use
  scale-aware numerical predicates where needed; preserve separate
  absolute and relative gap criteria.
- Keep zero-valued scores separate from missing mappings. Protect input
  arrays and AOF results from aliasing and in-place threshold side
  effects.
- Derived properties must be invalidated on geometry changes. Keep IDs
  and tie rules stable, and expose the selected method, parameters and
  backend in results so analyses can be reproduced.
- Numerical functions must run without plot windows,
  working-directory-specific files, implicit model downloads,
  unsolicited output files or network calls. Visualization and export
  are explicit operations.
- Keep memory use bounded and check dimensions/index capacity before
  allocating. Support interruption of long native kernels and produce
  actionable errors.

## Delivery sequence and acceptance gates

### Phase 0 — Baseline, rights inventory, and reference capture

- Freeze the source SHA, source manifest, function dependency map, and
  licenses.
- Inventory rights per reused core/helper file and per candidate
  dataset. The README explicitly defers licensing to individual modules;
  several core files have GPL-3.0-or-later notices, while many others
  only have copyright notices. Do not assume one blanket upstream
  license. Resolve the package license after the inventory and any
  necessary maintainer clarification.
- Establish Octave as the local reference runner, load image explicitly,
  and record runtime/package versions, platform, source SHA, parameters
  and seeds. Enumerate supported functions, unavailable operations and
  diagnostic adapters. Prefer `octave-cli` for headless runs.
- Compare selected runnable Python paths at their pinned SHA. Record
  independent analytical expectations alongside the two upstream
  outputs; behavior the two share must not become the expected R result
  when it departs from the analytical expectation.
- Capture MATLAB-specific reference evidence when available and useful
  for proprietary raster/morphology behavior. Its absence is a
  qualification boundary, not a blocker for the R package.
- Generate synthetic fixtures in the drawing representation and as
  binary matrices. Export intermediate arrays, event records, bin
  boundaries, and branch data, not only screenshots or final image
  summaries.
- Preserve NaN/Inf and dimensions in the fixture format; write a
  metadata manifest and convert to compact R fixtures reproducibly.
- Use `demos/testPatterns.svg` as a candidate workflow fixture after
  rights review. Treat `demos/Forest_DBW.mat` as a historical example,
  not a newly executed oracle: its generating version is not
  established.

**Gate:** reproducible analytical fixtures, an Octave/Python capability
matrix and executable reference probes, attribution/rights inventory for
the first reused scope, and a table of intentional differences.
Reference adapters must be visible and version-bound. Do not label an
Octave run as MATLAB evidence.

### Phase 1 — Useful contour-analysis package

- Create package metadata, S3 constructors/validators, and stable IDs.
- Implement orientation, length, curvature, junction
  detection/cleanup/types, histograms, summaries, and simple vector
  plotting.
- Define empty-drawing, zero-length, nonfinite coordinate and
  disconnected segment behavior. Reject invalid inputs or provide
  explicit cleaning; never permit unnoticed division by zero.
- Retain a simple exhaustive pairwise junction implementation as the
  small-case correctness witness. Profile representative sparse and
  dense drawings, then use exact spatial candidate filtering for larger
  inputs when it improves cost. Candidate rejection must account for the
  existing gap tolerances.

**Gate:** independent analytical fixtures and qualified upstream
comparisons for supported cases, regression tests for each documented
difference, and a runnable “drawing to feature table” example. This is
the first useful release scope even if the MAT pipeline is still under
development.

### Phase 2 — Raster and MAT feasibility gate

- Implement deterministic rasterization separately from display
  plotting.
- Implement or qualify exact Euclidean distance plus nearest-feature
  indices, including tie rules; preserve the signed-distance
  construction.
- Port AOF with the specified sampling and boundaries.
- Reproduce or explicitly qualify thinning, component cleanup, edge
  linking, and branch construction. The primary tracing route is
  `GetConSeg` -\> `edgelink(image,10)` -\> cleanup, with its
  minimum-length and topology effects.
- Profile the pilot on modest and larger synthetic images before
  expanding native-code scope.

**Gate:** compare each intermediate array on rectangles, disks, rings,
parallel lines, border-touching shapes and junctions. Test raster masks,
topology and geometric accuracy, not merely visual plausibility. Label
the morphology backend. When a deliberate algorithm change alters
topology, establish its scientific behavior with independent evidence
and document the method difference; do not hide it by loosening
comparison tolerances.

### Phase 3 — Full medial-axis analysis

- Implement parallelism, separation, mirror and taper scores with
  explicit smoothing/window semantics.
- Implement tangent projection, nearest-neighbor selection, mapping
  masks, per-contour scores and MAT histograms.
- Port the full drawing -\> raster -\> MAT -\> branch scores -\> contour
  scores -\> feature-table workflow.

**Gate:** separate fixtures for branch formulas and geometric
reconstruction, then end-to-end comparisons. Compare branch structure,
score locations, missingness and histograms in addition to aggregate
means.

### Phase 4 — Complete the practical workflows

- Add SVG and MATLAB drawing import, raster edge-map tracing, rotation,
  circular apertures, shifting, segmentation and property-based
  splitting.
- Preserve the existing SVG feature set deliberately: lines, polylines,
  polygons, rectangles, circles/ellipses, Bezier paths, arcs and
  transforms. Stage this support; unsupported constructs must be
  reported rather than silently dropped. General SVG rendering is
  outside scope.
- Implement model feature-table construction and R prediction adapters.
- Add property plots, image exports, feature-density maps and batch
  examples. Consolidate MATLAB draw/render variants into a small R
  surface.
- Make stochastic operations reproducible with R seeds and recorded
  settings; cross-language random sequences are not an initial
  equivalence requirement.

**Gate:** contour IDs/order remain traceable; transformations invalidate
derived properties correctly; split lengths, ties and model schemas are
tested; imports have a documented support matrix.

### Phase 5 — Optional RGB backend and distribution qualification

These are separate work streams; RGB support does not gate an analysis
release.

- RGB: establish a backend contract for segmentation/edge maps, then
  qualify model identity, weights, preprocessing, thresholds, output
  masks and the shared tracing/postprocessing path. Review pyMLV as it
  develops; its inspected tree does not yet supply a complete extraction
  backend. A Python SAM bridge remains a candidate requiring its own
  qualification against the stated method.
- Distribution: run package checks on macOS, Windows and Linux, native
  checks if C++ is included, meaningful coverage, and build the
  vignettes from scratch.
- Publish documentation for contour analysis, medial-axis analysis,
  stimulus manipulation and optional extraction, with a feature/status
  matrix.

**Gate:** executed platform checks and reproducible examples for the
advertised scope. Publication, registry submission and deployment remain
separate actions.

## Verification design

Use both upstream comparisons and independent evidence; reproducing
upstream behavior on a documented edge case is not scientific
validation.

- Geometry fixtures: horizontal/vertical/diagonal and subpixel segments,
  directed reversal, one-segment contours, corners, closed loops,
  zero-length segments, empty drawings and large/small coordinate
  magnitudes.
- Junction fixtures: T, Y, X, Arrow, Star, near misses around both
  epsilon thresholds, collinear/overlapping segments, nearby junction
  merges and same-contour self-crossings (excluded by the upstream
  method).
- Raster/MAT fixtures: exact small masks, unequal image dimensions, tiny
  images, blank/full masks, borders, symmetric nearest-boundary ties,
  short branches, loops and adjacent branch points.
- Score fixtures: constant-radius straight branches, changing radius,
  bent branches, endpoints and smoothing windows shorter than the data;
  test the final exponent explicitly.
- Statistics: exact bin boundaries, mass conservation for applicable
  histograms, zero totals, nondefault bin counts and normalization
  units.
- Metamorphic checks: translation preserves vector geometry statistics;
  rigid rotation preserves lengths and curvature; contour reversal
  preserves axial orientation histograms. Scaling preserves orientation
  and scales length/curvature inversely as appropriate. Do not assert
  these invariances indiscriminately for pixel-threshold junctions or
  raster-based scores.
- Compare topology and identifiers exactly after a documented
  canonicalization of branch direction/order. Use explicit per-output
  floating-point tolerances selected from reference precision and
  sensitivity tests, not one broad tolerance that can hide skeleton
  changes.
- Measure core geometry, junction search, rasterization, distance/AOF,
  branch tracing and mapping separately. Report input dimensions,
  contour counts, runtime versions, cold versus warm behavior and memory
  use. Do not promise a speedup before measurement.

## Performance work and acceptance criteria

A faster implementation is an explicit deliverable. Start with correct
small reference kernels and measured upstream baselines. Improve
algorithms and data movement before introducing parallelism.

| Hot path | Candidate improvement | Correctness witness |
|----|----|----|
| Junction detection | Expanded bounding boxes plus a grid/sweep/spatial index, then the exact intersection predicate | Exhaustive pairs on small random and adversarial drawings; sparse and dense cases |
| Distance transforms | Linear-time exact EDT with nearest-feature coordinates; bounded working arrays | Brute-force nearest-boundary search on small masks, including ties |
| AOF | Compiled loops, precomputed samples/normal data, efficient index decoding and reusable buffers | Slow scalar evaluation on tiny signed-distance/nearest-boundary inputs |
| Branch scores | Prefix sums for fixed/local windows and shared derivative computations | Direct window sums and analytical branches, including endpoints |
| Contour projection | Exact nearest-neighbor indexing with deterministic ties; reuse geometry across properties | Brute-force nearest-neighbor/tangent mapping on small drawings |
| Raster and feature tables | Avoid drawing each full image repeatedly; reuse masks and preallocate tables | Pixel fixtures, ID tracing and histogram conservation |

These are candidates to measure, not established speedups. Preserve
exactness where the contract requires it. Approximate search, lower
image resolution, fewer AOF samples and lower numerical precision
require an explicitly named method/quality tradeoff; they cannot be
hidden inside a speed claim.

Plan benchmark families around 100, 1,000 and 10,000 segments (including
sparse and dense junction layouts), and 128, 512 and 1,024 pixel square
masks plus rectangular cases. Bound the expensive reference runs and
record any cases that cannot finish within the benchmark budget.

Compare identical completed work and supported parameters against the
pinned MATLAB source under Octave, runnable Python paths where
available, and the simple R reference implementation. Record host,
versions, threads, setup, cold/first/warm timings, allocations/peak
memory, and per-stage as well as end-to-end time. Missing or failing
upstream functions are capability gaps, not speedups. Separate
extraction/model loading from contour analysis.

**Gate:** on declared representative workloads, demonstrate an
end-to-end improvement and explain per-stage regressions or tradeoffs,
while passing the correctness suite and memory bounds. Report measured
ratios and uncertainty; set numerical performance targets after baseline
profiling. Do not claim the package is faster merely because it uses
compiled code.

Once a package exists, the minimum local checks will be the relevant
`testthat` subset followed by `R CMD build` and `R CMD check --as-cran`
on the produced tarball. Exact package commands and platform jobs will
be recorded against the implemented revision. These checks have not been
run because no R package has been implemented.

## Data and dependency boundaries

Keep the user-supplied reference checkouts available locally and
unchanged. Exclude full vendor checkouts, the approximately 651 MiB
dataset directory and trained-model trees from the R package tarball and
normal CI. Configure package build exclusions once the layout settles;
this does not authorize deleting or moving the user’s reference data.

Use synthetic fixtures by default; add small real-image examples only
after checking redistribution rights. Availability in an upstream
repository does not establish permission for a new package distribution.

Local assessment environment: R 4.5.1; Rcpp 1.1.2, cpp11 0.5.5, testthat
3.3.2, xml2 1.6.0, png 0.1-9, jpeg 0.1-11, R.matlab 3.7.0, imager 1.0.8
and hdf5r 1.3.12 were reported installed. Package presence does not
establish backend suitability. R startup emitted locale warnings. Octave
11.1.0 with image 2.20.0 is available and was exercised as described
above. Python 3.14.7 with NumPy 2.4.3 ran the orientation probe; cv2 and
skimage were absent from that interpreter. No MATLAB executable was
found on PATH or in the standard app location. Full MATLAB equivalence
remains untested.

Initial sizing: contour analysis is the smallest work package;
rasterization, skeleton topology and MAT projection are the largest
uncertainty. Set a calendar estimate after the Phase 0 reference capture
and a Phase 2 pilot rather than estimating solely from source line
counts. RGB backend work has its own model and runtime dependencies and
should be estimated separately.

## First implementation task

Start with manifests for both reference checkouts, the Octave
runner/adaptation specification, the regression cases for known
differences and the `mlv_drawing` contract. Implement constructors and
orientation/length/curvature with subpixel and degenerate-input tests,
then define histogram boundaries and junction output shapes as those
APIs are added. Establish the correctness and timing harness before
optimizing the expensive kernels.

This assessment updated planning/evidence documents and ran bounded
synthetic probes. The user supplied the upstream checkouts; no R
package, model download or publication was created.

## Complete top-level function inventory

Phase assignments describe implementation work packages, not a promise
of a one-to-one exported R function. Helpers can remain internal;
presentation wrappers should be consolidated.

### 1 — Contour geometry and statistics (17 files)

`averageProperty`, `cleanupJunctions`, `computeContourProperties`,
`computeCurvature`, `computeJunctionAnglesTypes`, `computeJunctions`,
`computeLength`, `computeOrientation`, `detectJunctions`,
`getContourPropertiesStats`, `getCurvatureStats`,
`getHorizontalVerticalStats`, `getJunctionStats`, `getLengthStats`,
`getOrientationStats`, `lineIntersection`, `removeZeroLengthContours`.

### 2 — Raster and medial-axis feasibility (7 files)

`InitializeNeighborhoods`, `computeAOF`, `computeGradientVectorField`,
`computeMAT`, `extract2DSkeletonFromBinaryImage`, `renderLinedrawing`,
`traceSkeleton`.

### 3 — Medial-axis scores and projection (9 files)

`MATpropertiesToContours`, `computeAllMATfromVecLD`,
`computeAllMATproperties`, `computeMATproperty`,
`computeMATpropertyPerBranch`, `fitLineSegments`,
`getDistanceFromLineSegment`, `getMATpropertyStats`, `mapMATtoContour`.

### 4 — Import, manipulation, models, and presentation (28 files)

`allLDHistogramsToTable`, `applyCircularAperture`, `computeColorIndex`,
`drawAllProperties`, `drawContoursByStatsModel`, `drawJunctions`,
`drawLinedrawing`, `drawLinedrawingProperty`, `drawMATproperty`,
`generateFeatureDensityMap`, `histogramToTable`, `importSVG`,
`labelMatrix2edges`, `mergeLineSegments`, `predictContoursByStatsModel`,
`randomlyShiftContours`, `renderContoursByStatsModel`,
`renderJunctions`, `renderLinedrawingProperty`, `renderMATproperty`,
`rotateLinedrawing`, `saveSceneLDsToTable`,
`segmentContoursAtJunctions`, `splitLDbyHistogramWeights`,
`splitLDbyProperties`, `splitLDbyStatsModel`,
`splitLDmiddleSegmentsVsJunctions`, `traceLinedrawingFromEdgeMap`.

### 5 — Optional RGB extraction (1 files)

`traceLineDrawingFromRGB`.

Additional reachable dependencies to inventory and qualify include
`3rdParty/contour_tracing/main/GetConSeg.m`,
`3rdParty/contour_tracing/edge_linking/edgelink.m`, `cleanedgelist.m`
and `findendsjunctions.m`. The vendored edge-linking source contains a
permissive notice requiring retention of its copyright and permission
notice; review each reused helper individually.

## Primary source links

- [Pinned Python
  port](https://github.com/bwlabToronto/pyMLV/tree/14b3183b91146b1021eb7ef97db75ad703b06832)

- [Python rendering
  limitation](https://github.com/bwlabToronto/pyMLV/blob/14b3183b91146b1021eb7ef97db75ad703b06832/README.md)

- [Python orientation
  implementation](https://github.com/bwlabToronto/pyMLV/blob/14b3183b91146b1021eb7ef97db75ad703b06832/MLVcode/computeOrientation.py)

- [Pinned README and
  requirements](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/README.md)

- [Pinned MATLAB
  core](https://github.com/bwlabToronto/MLV_toolbox/tree/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode)

- [Published
  methods](https://www.frontiersin.org/journals/computer-science/articles/10.3389/fcomp.2023.1140723/full)

- [Distance and boundary
  indices](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/computeGradientVectorField.m)

- [Average outward
  flux](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/computeAOF.m)

- [Deterministic upstream rasterization
  settings](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/renderLinedrawing.m)

- [Principal branch-tracing
  wrapper](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/3rdParty/contour_tracing/main/GetConSeg.m)

- [Edge-linking implementation and
  notice](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/3rdParty/contour_tracing/edge_linking/edgelink.m)

- [Experimental SVG
  importer](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/MLVcode/importSVG.m)

- [Paper figure
  workflow](https://github.com/bwlabToronto/MLV_toolbox/blob/f12418a5e03c2f38e222fed1dffccf42b227cd29/demos/Frontiers2023_Figures.m)
