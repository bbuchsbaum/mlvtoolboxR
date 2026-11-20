# Differences from the reference implementations

This package ports and extends the Mid-Level Vision (MLV) Toolbox
developed in Dirk Bernhardt-Walther’s lab (MATLAB) and its Python port,
pyMLV. This register records, under stable IDs, the observed differences
between this package and the pinned reference sources (what each
implementation does, and why), the package’s own method choices, and
revisions made to the R port. The `mlv2023` options reproduce the
reference behaviour except where an entry says otherwise.

Pinned references: MATLAB
[f12418a5e03c2f38e222fed1dffccf42b227cd29](https://github.com/bwlabToronto/MLV_toolbox/commit/f12418a5e03c2f38e222fed1dffccf42b227cd29);
Python
[14b3183b91146b1021eb7ef97db75ad703b06832](https://github.com/bwlabToronto/pyMLV/commit/14b3183b91146b1021eb7ef97db75ad703b06832).

None of these entries has been submitted upstream. The notes distinguish
executed reproducers, source-confirmed differences, portability issues
and open questions about the intended scientific definition. Agreement
between the MATLAB and Python ports reflects their shared lineage; it
does not by itself settle which definition is intended. All reproducible
inputs here are synthetic.

## Executed MATLAB-source edge cases

Run `octave-cli --quiet --no-history tools/upstream/reproduce_matlab.m`.
Executed output is retained in `docs/evidence/upstream-matlab.log`. The
interpreter was Octave 11.1.0; none of these contour/geometry
reproductions requires a MATLAB proprietary toolbox.

| ID | Source | Minimal trigger and result in the pinned source | Behaviour of this package / R regression |
|----|----|----|----|
| MLV-M001 | getHorizontalVerticalStats.m | One horizontal segment of length 10: histogram total 0 | R includes the upper endpoint (+1) and conserves the 10 pixels. test-geometry.R |
| MLV-M002 | getJunctionStats.m | No junctions, eight angle bins: value vector length 5, centers length 8 | R returns eight zeros, also for nondefault requested bin counts. test-geometry.R |
| MLV-M003 | splitLDbyProperties.m and related splitters | Lengths 4,4,2, fraction .5: top set length 8/10 | In R both ranked prefixes respect the five-pixel budget. In the pinned source the same cumulative comparison can also retain a last contour at fraction zero. test-geometry.R |
| MLV-M005 | applyCircularAperture.m | Segment (0,5)-(10,5), radius 2 at (5,5): all geometry discarded | R keeps (3,5)-(7,5), length 4: a segment can cross the circle while both of its endpoints lie outside it. test-geometry.R |
| MLV-M006 | cleanupJunctions.m | Positions (0,0),(4,0),(2,0), merge radius 2: two groups | R forms one transitive connected component. The pinned source traverses only the upper-triangular adjacency, so backward connections are not followed. test-geometry.R |
| MLV-M007 | computeJunctionAnglesTypes.m | Diagonal bar (0,0)-(10,10), perpendicular stem (5,5)-(10,0): T minimum angle 0 degrees | R reports a minimum angle of 90 degrees. In the pinned source, interior arms use a different y convention from endpoint arms. test-geometry.R |
| MLV-M008 | lineIntersection.m | Crossing negative-coordinate segments around (-3,-3): empty result | R returns the intersection (-3,-3); it does not require geometry to lie inside the image. test-geometry.R |
| MLV-M009 | fitLineSegments.m | Two-point straight line: the pinned source stops with an undefined-lineSegs error | R returns a single span with defined scores. test-kernels.R |

## Other MATLAB findings

| ID | Status | Source and evidence | R treatment |
|----|----|----|----|
| MLV-M004 | Source + executed arithmetic; output sensitivity not established broadly | computeMAT.m uses .3 when threshold is omitted, but explicit 28 converts to .29887487943380647 | Omitted and explicit threshold_angle=28 use the same conversion; retain parameters |
| MLV-M010 | Executed Octave portability issue; not shown to affect MATLAB | computeGradientVectorField.m has two jj++m_Neighbors8 expressions rejected by Octave | Reference harness uses a temporary, SHA-bound explicit-addition adaptation outside vendor |
| MLV-M011 | Source-confirmed missingness convention | MATpropertiesToContours.m marks every score equal to zero as invalid and replaces it with NaN | Carry a mapping mask; retain valid zeros; test-kernels.R |
| MLV-M012 | Documentation and ranking differ; requires scientific clarification | splitLDbyProperties.m describes mean curvature but ranks by the unnormalized length-weighted sum | Expose curvature (weighted mean) and integrated_curvature as distinct criteria |
| MLV-M013 | Source-confirmed global-state side effect | lineIntersection.m disables MATLAB:divideByZero warnings without restoring their state | Native predicates handle parallelism without changing R warning options |
| MLV-M014 | Documentation discrepancy | traceLineDrawingFromRGB.m help says structured-edge default; arguments block selects SAM | RGB extraction provider is explicit and recorded; no implicit model download |
| MLV-M015 | Source-confirmed degenerate-input case; needs a valid-domain policy | computeCurvature.m divides turning angle by zero segment length on a multi-segment contour | Constructor rejects zero-length geometry or explicitly drops it, retaining contour IDs |

## Python findings

Run
`PYTHONDONTWRITEBYTECODE=1 python3 tools/upstream/reproduce_python.py`.
Executed output is retained in `docs/evidence/upstream-python.log`.
Python 3.14.7 and NumPy 2.4.3 were used. cv2/skimage were unavailable in
that interpreter; the full image pipeline was not executed.

| ID | Status | Source and minimal evidence | R treatment |
|----|----|----|----|
| MLV-P001 | Executed function | computeOrientation.py truncates endpoints (1.1,1.1)-(1.9,1.9) to int32 and returns 0 degrees | Preserve subpixel doubles; analytical result 315 degrees |
| MLV-P002 | Executed KeyError + further source differences | getHorizontalVerticalStats.py reads the field orientation where the upstream schema provides orientations; it also applies radian trig to degree values and builds one fewer bin edge than its loop expects | Use one consistent schema, degree conversion, explicit complete boundaries and conservation tests |
| MLV-P003 | Source/module inventory | Branch functions import smoothData, diff, fitLineSegments and GetConSeg modules that are not present in the pinned repository. traceSkeleton.py also returns R/F fields while downstream expects Radius/AOF | Native graph and explicit branch schema; no missing helper dependency |
| MLV-P004 | Source-confirmed alias + executed NumPy semantics witness, not full-pipeline execution | extract2DSkeletonFromBinaryImage.py sets skeletonImage=fluxImage then thresholds in place, so returned AOF is overwritten | Separate flux/skeleton buffers; test input/result immutability |
| MLV-P005 | Source-confirmed file-based external state and side effects | computeGradientVectorField.py replaces computed boundaries with outerBoundaryOriginal.mat and opens plots; computeMAT.py unconditionally writes MAT.mat | Numerical APIs require no working-directory files, plot windows or unsolicited output files |
| MLV-P006 | Source-confirmed: taper not complete; parallelism formula differs from MATLAB | computeMATpropertyPerBranch.py uses newddR before assignment; parallelism uses cumulative sums rather than the MATLAB local derivative magnitudes; taper uses dR where a second derivative is described | Explicit branch formulas tested on constant-radius and tapered branches |
| MLV-P007 | Source-confirmed array handling | computeMAT.py unconditionally selects binaryImage\[:,:,0\], rejecting a two-dimensional input described by its API | Binary matrices are an explicit supported input |
| MLV-P008 | Source-confirmed: projection step not runnable as written | MATpropertiesToContours.py initializes an ndarray then calls append; loop uses numContours inconsistently, and coordinate bookkeeping stores only row indices before unravel_index | Typed contour/score structures, exact indices and explicit mapping masks |
| MLV-P009 | Source-confirmed substituted cleanup parameter | extract2DSkeletonFromBinaryImage.py replaces calculated area_threshold with 100; the thinned array and returned skeleton are also handled inconsistently | Named, recorded component threshold and an explicitly thinned returned skeleton |
| MLV-P010 | Source-confirmed padding/smoothing discrepancy | getBranchDerivative.py leaves dR shorter than dX/dY for nontrivial branches and uses a radius moving-mean width of three versus MATLAB’s two | Equal-sized padded derivatives; radius window two, with explicit endpoint convention |
| MLV-P011 | Shared empty-result edge case plus source count semantics | getJunctionStats.py shares the five-versus-eight empty-histogram shape of MLV-M002, and passes centers as np.histogram edges on nonempty input | Explicit event histogram edges and stable requested dimensions |

## Intentional method choices

These are design choices of this package, not observations about the
reference sources.

- Thinning (0.3.0): the default `morphology = "minimal"` is the
  MATLAB-exact sequential hit-or-miss thinning (`"hmt8"`, identified on
  Forest_DBW) followed by removal of remaining simple pixels; `"hmt8"`
  alone is the published-compatible choice and `"pratt"` is the
  Octave-image backend used before 0.3.0.
- Exact EDT ties select the lowest R linear index. Other equally near
  choices can change AOF values. The reference report measures such
  differences.
- Contour tracing uses an explicit pixel graph without diagonal
  shortcuts across an existing orthogonal connection. It is not a claim
  to reproduce Kovesi’s spur cleanup and contour ordering.
- Width-one rasterization uses clipped Bresenham lines. MATLAB
  insertShape equivalence is not claimed.
- Taper uses explicit moving means instead of unspecified smoothdata
  defaults.
- average_property(…,“orientation”) is an explicitly documented axial
  circular mean, with NA for an isotropic resultant; it replaces the
  original concatenate-and-flip vector summary.
- Per-contour junction participation counts a junction once per contour.
- SVG curve_step controls sampling, not a certified geometric error
  bound.

Before raising any entry with the MLV maintainers, rerun the minimal
reproducer at their current revision and, where the observation is
Octave-specific, confirm whether it also occurs in MATLAB. Attach source
SHA, runtime, minimal input, observed and expected results, and the
scientific consequence. Do not present source-only observations as
executed full-pipeline results.

## Findings from the expanded parity suite (0.1.1)

The fixtures in `inst/fixtures/parity` and `tools/parity/generate.py`
execute the pinned sources on curved, integer and subpixel contours,
several bin counts, varying-radius branches and six image topologies.
Temporary adapters are recorded in the fixture manifest.
`tools/parity/compare.R` retains differences explicitly.

| ID | Status and trigger | Finding and R treatment |
|----|----|----|
| MLV-M016 | Executed terminal-bend branch; points (1,0) through (20,0), then (21,100) | fitLineSegments assigns the last one or two points to the final span even when their mean perpendicular distance exceeds its threshold of 1. The returned all-point span has mean interior distance about 9.81. This also changes mirror scores on curved branches. R ends the span before the point that exceeds the threshold and tests every returned chord against an independent projection calculation. |
| MLV-M017 | Executed 13-bin axial histogram, orientations 90 and 270 degrees | getOrientationStats shifts, wraps and subtracts a half-bin; rounding puts axial-equivalent vertical directions into different bins. R scales the axial angle directly and assigns an exact boundary to its upper bin. Agreement is retained on ordinary 5/8-bin fixtures; the 13-bin boundary divergence is intentional and tested analytically. |
| MLV-P012 | Executed AOF on spatially varying nearest labels, radii 1 and 2 | computeAOF indexes IDX at int(px-1), int(py-1) after constructing zero-based sampling coordinates; this selects the preceding grid cell. It also omits epsilon from the sampling displacement. With the same distance/label fields, flux values differ from the MATLAB source by more than 45 at radius 1. R agrees with the MATLAB source within 1e-10 on both radii. |
| MLV-P013 | Executed 37-degree rotation | rotateLinedrawing rounds sine and cosine to five decimal places. On synthetic input, coordinates differ from the full-precision rotation by up to approximately 0.000895 pixels, and the transform is not exactly length-preserving. R uses full precision; its counterclockwise image-coordinate angle corresponds to the negative upstream angle. |
| MLV-P014 | Executed raw getLengthStats after computeLength; NumPy 2.4.3 | computeLength produces contourLengths with shape (n,1), but getLengthStats adds its one-element arrays into scalar histogram cells. Current NumPy raises ValueError: setting an array element with a sequence. The harness records the raw exception and separately labels an adapted comparison that reshapes only this field to (n,). |
| MLV-P015 | Executed 13-bin histograms; shared formula with MLV-M017 | Python shares the axial boundary sensitivity of the MATLAB half-bin shift/subtraction formula. R’s explicit upper-bin rule applies consistently after reversing segment direction. |

The retained Forest_DBW demonstration contains one zero-length contour.
The optional natural-drawing harness removes that one segment/contour
explicitly and identically for all three implementations, retaining 562
contours and 2,072 segments. This is a real-data witness for the already
recorded MLV-M015 degenerate-input case; no source snapshot or canonical
data file is changed.

## R-port revisions in 0.1.1

| ID | Earlier R behaviour (executed) | Revision and regression |
|----|----|----|
| MLV-R001 | Sparse/isolated masks: the original R port assigned distance 0 instead of -1 to isolated interior ink pixels. Its outer-boundary predicate checked only for a white neighbor. | Require both an ink and white neighbor, matching the signed-distance construction. Six diverse masks now agree on distances within the single-precision reference tolerance; isolated pixels have analytical distance -1. |
| MLV-R002 | The original R port placed an exact vertical 13-bin boundary in the lower bin because its half-bin shift rounded downward. | Use the axial angle scaled by the bin count before rounding to a centered bin. Tests cover reversed vertical directions and odd bin counts 1, 3, 5, 7, 13 and 31. |

## Observed method differences, not numerical parity

Equidistant nearest-boundary pixels need not have the same label in R
and Octave. The R EDT selects the lowest column-major index. All
differing labels in the six-mask comparison have the same exact squared
distance after the MLV-R001 revision, but their directions can change
AOF and the final skeleton. Tests therefore establish full distance
parity, AOF parity with identical distance/label fields, and morphology
parity with identical flux fields. They do not label the differing
default end-to-end skeletons as matches. See [the complete parity
report](https://bbuchsbaum.github.io/mlvtoolboxR/docs/PARITY.md) for
counts and images. Keeping R distances but substituting Octave nearest
labels makes all six skeletons identical, with maximum AOF difference
below 1.5e-14. Actual MATLAB’s distance/morphology tie behavior remains
unexecuted.

## Pre-release native integration review (0.2.1)

These are revisions to the new R implementation made during review; they
are not observations about either vendored library.

| ID | Evidence and boundary | Revision |
|----|----|----|
| MLV-R003 | Source review: native contour-column offsets used an integer row stride; sufficiently large matrices could overflow the multiplied offset. No huge-memory or crash reproduction was executed. | Use R_xlen_t for pointer strides and offsets. Ordinary integer/double contour and scale cases execute in the native regression suite; the huge allocation remains unexecuted. |
| MLV-R004 | Source review: malformed manually edited S3 contour geometry could produce nonfinite lengths/orientations before native bin indexing. Normal constructors already reject these inputs. No crash was executed before the revision. | Reject zero, NA and infinite segment lengths before bin indexing. Executed public-API regression tests confirm all three edited inputs fail with an explicit diagnostic. |

See [the integrated performance
report](https://bbuchsbaum.github.io/mlvtoolboxR/docs/FUSED_STATISTICS.md)
for the exact qualified archive, numerical witnesses and check evidence.

## Findings from qualification against stored MATLAB outputs (0.3.0)

The lab’s MLV datasets in the local reference collection carry outputs
computed in MATLAB: 3,017 drawings with orientations, lengths,
curvatures and junctions, 2,757 of them with per-pixel parallelism,
separation and mirror scores, and the full Forest_DBW distance map, AOF
and skeleton. Re-executing the pinned source under Octave reproduces the
stored fields, so these findings apply to the pinned revision. Each
entry has a synthetic reproducer with an independently derived answer;
full evidence, counts and oracles are in [the junction
report](https://bbuchsbaum.github.io/mlvtoolboxR/docs/evidence/mlvref/junctions/REPORT.md)
and [the medial-axis
report](https://bbuchsbaum.github.io/mlvtoolboxR/docs/evidence/mlvref/mat/REPORT.md).
Natural data are not redistributed.

### MLV findings

| ID | Finding in the pinned source | Natural effect | Behaviour of this package |
|----|----|----|----|
| MLV-M018 | Contract clarification: lineIntersection.m converts the 1 px gap allowance to a parameter fraction with max(abs(dx),abs(dy)), so the absolute gap is Chebyshev (up to 1.41 px on a diagonal) | 7,745 junctions change under a Euclidean gap | Default `gap_metric = "euclidean"` (rotation-invariant); `"chebyshev"` available |
| MLV-M019 | Intersection parameters within 1e-4 of 0 or 1 snap to the segment endpoint, a displacement of up to 1e-4 × length | 29 junctions | Exact intersection |
| MLV-M020 | A zero-length neighbour adds an arm of orientation atan2d(0,0) = 0 | 8 junctions | Zero-length segments are rejected at construction or dropped |
| MLV-M021 | Interior arms are measured from the merged centroid rather than along the leaving segment, so straight-arm angles shift with the centroid offset | 960 junctions | Arm direction from the segment through which the contour leaves the junction |
| MLV-M022 | Segments attach at their start when it is within 2 px, even when the end is nearer, so reversing a contour changes the junction | part of 7,484 | Direction-free exit rule |
| MLV-M023 | Splitting one polyline into contours with a short middle piece creates a junction | part of 7,484; 376 junctions at contour split points | Contours ending inside the neighbourhood contribute no arm |
| MLV-M024 | Every segment of a contour inside the neighbourhood is an arm (duplicate arms; overshoots read as X) | 1,255 exact duplicates; most of 5,659 Star→X | One arm per contour exit |
| MLV-M025 | A member segment with both endpoints beyond 2 px always gets two arms, even when it lies wholly beyond the junction | part of 7,484 | Only the arm leading away |
| MLV-M026 | Closed contours are treated as ending at their first vertex; the junction depends on the start vertex | 68 arm-count changes; 164 junctions MLV does not report | Closed contours wrap |
| MLV-M027 | computeAOF centres its sampling circle on the pixel corner and reads labels from floor() cells, a half-pixel offset between the two grids; the result is not lattice-equivariant, and two parallel lines give a two-row ridge | 3,285 Forest skeleton pixels | Default `aof = "centered"` (pixel-centred, mean over equidistant boundary pixels) |
| MLV-M028 | Branch scores depend on traversal direction (forward differences, movmean(2) radius, pair (i,i+1) scored by i) | Reversing branches drops stored-value agreement from 0.999 to 0.044 | Default scores are direction-symmetric |
| MLV-M029 | mapMATtoContour reads a last-writer-wins score image and resolves exact nearest-tangent ties by enumeration order | Forest agreement 0.848 → 0.999 when emulated | Each branch keeps its own score; maximum over exact ties |
| MLV-M030 | Kovesi cleanedgelist duplicates junction points when merging edges | 272 duplicates on Forest | Default tracing never repeats a junction pixel within a branch |
| MLV-M031 | separation = (1 - 1/R)^10 is unbounded: radius 0 gives +Inf, radii below 0.5 exceed 1 | 11 stored +Inf samples | NA below radius 1 |
| MLV-M032 | computeCurvature.m gives each segment the whole turn to the next segment divided by its own length; the last segment repeats the previous turn. The last vertex of every contour is counted twice (a two-segment roof with a 90 degree apex integrates to 180 degrees), a closed contour’s closing vertex is not counted (a right triangle integrates to 405 degrees), and a vertex’s value depends on one neighbouring length, so resampling changes it (a circle of radius 10 with alternating 1:2 vertex spacing falls into two 8-bin curvature bins) | Median ratio of integrated curvature, published / each-turn-once, 1.22 over 3,017 drawings; median L1 between normalized 8-bin curvature histograms 1,650 per 10,000 | Default `estimator = "turning"` (each turn spread once over its dual cell); `"mlv2023"` reproduces the published values; test-default-estimator-snapping.R, test-invariance-contract.R |
| MLV-M033 | computeMATpropertyPerBranch scores undefined points 0: the two end points of parallelism and taper (zero-width derivative window), every point of a branch with fewer than three points, and mirror on three or fewer points; the width-3 moving mean then averages that 0 into the neighbouring point (a straight constant-radius branch scores 0.5^10 and (2/3)^10 at its first two points instead of 1) | Every branch end on the 2,757 stored drawings; applying this change alone lowers per-contour parallelism r between default and published values from 0.926 to 0.811 | NA where undefined; smoothing averages defined values only (`branch_rule = "invariant"`) |

MLV-M006, M007, M008 and M015 were re-verified on the natural data (462,
24,759, 0 and 552 affected items respectively). `rule = "mlv2023"` and
`aof = "mlv2023"`, `branch_rule = "mlv2023"` reproduce all of these
behaviours for published-compatible results.

### R-port revisions

| ID | Earlier R behaviour | Revision |
|----|----|----|
| MLV-R005 | Junction arm direction was the chord from the merged position, biasing straight-arm angles (28,975 natural junctions) | Direction of the leaving segment |
| MLV-R006 | Junction arm membership depended on contour direction; contour pieces inside the neighbourhood counted as arms | Arms are exits from the open disc of radius `endpoint_radius` |
| MLV-R007 | Closed contours did not wrap; junction type depended on the start vertex (247 natural junctions) | Indices wrap for closed contours |
| MLV-R008 | A member segment that never reaches the disc got an additional arm back across the junction (172 natural junctions) | Only the arm leading away |
| MLV-R009 | Rounding-residue segments (length ~1e-13) were accepted, producing curvatures ~1e13 deg/px (19 natural segments) | Endpoints within 64 ulp of their magnitude are the same point |
| MLV-R010 | Width-1 rasterization depended on segment direction (1,984 Forest pixels changed under reversal) | Major-axis DDA with MATLAB’s half-pixel tie rule; reversal-invariant, 0 pixels missing against MATLAB’s raster |
| MLV-R011 | Coincident consecutive branch points produced NaN tangent subscripts and stopped | Zero-length pairs contribute no tangent samples |
| MLV-R012 | Skeleton pixels on ink (distance \<= 0) reached branch scoring and aborted it (7 of 580 drawings, Pratt path) | Excluded before tracing on every path, with a recorded count |
| MLV-R013 | The Octave-derived Pratt thinning is neither thin nor MATLAB’s skeleton (18,035 redundant pixels vs 327 on Forest) | `morphology = "hmt8"` reproduces MATLAB’s Forest skeleton exactly and is the published-compatible choice; the default `"minimal"` adds removal of remaining simple pixels; Pratt kept as an option |
| MLV-R014 | A merged junction’s neighbourhood was one disc around the centroid; distant merged crossings lost or gained arms | Union of discs around every merged intersection |
| MLV-R015 | Sweep prefilter was not conservative and exact-limit extensions depended on segment direction | Canonical segment frame, shared acceptance arithmetic; sweep and exhaustive bitwise identical |
| MLV-R016 | Dropping a rounding-residue segment left a gap that opened closed contours | Kept segments are re-joined exactly, cyclically for closed contours |
| MLV-R017 | `rule = "mlv2023"` visited merge partners in sweep order, not MATLAB’s ascending index | Ascending order; iterative traversal |
