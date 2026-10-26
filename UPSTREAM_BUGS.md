# Upstream bug register

Pinned references: MATLAB [f12418a5e03c2f38e222fed1dffccf42b227cd29](https://github.com/bwlabToronto/MLV_toolbox/commit/f12418a5e03c2f38e222fed1dffccf42b227cd29);
Python [14b3183b91146b1021eb7ef97db75ad703b06832](https://github.com/bwlabToronto/pyMLV/commit/14b3183b91146b1021eb7ef97db75ad703b06832).

No upstream issues have been submitted. These notes distinguish executed
reproducers, source-confirmed defects, portability issues and unresolved
scientific contracts. Shared lineage and agreement between ports do not
establish correctness. All reproducible inputs here are synthetic.

## Executed MATLAB-source defects

Run `octave-cli --quiet --no-history tools/upstream/reproduce_matlab.m`.
Executed output is retained in `docs/evidence/upstream-matlab.log`.
The interpreter was Octave 11.1.0; none of these contour/geometry reproductions
requires a MATLAB proprietary toolbox.

| ID | Source | Minimal trigger and observed result | Expected behavior / R regression |
|---|---|---|---|
| MLV-M001 | getHorizontalVerticalStats.m | One horizontal segment of length 10: histogram total 0 | Include the upper endpoint (+1); conserve 10 pixels. test-geometry.R |
| MLV-M002 | getJunctionStats.m | No junctions, eight angle bins: value vector length 5, centers length 8 | Eight zeros, including nondefault requested bin counts. test-geometry.R |
| MLV-M003 | splitLDbyProperties.m and related splitters | Lengths 4,4,2, fraction .5: top set length 8/10 | Both ranked prefixes obey a five-pixel budget. The same cumulative comparison can retain a last contour at fraction zero. test-geometry.R |
| MLV-M005 | applyCircularAperture.m | Segment (0,5)-(10,5), radius 2 at (5,5): all geometry discarded | Keep (3,5)-(7,5), length 4. Both endpoints may be outside a circle while the segment crosses it. test-geometry.R |
| MLV-M006 | cleanupJunctions.m | Positions (0,0),(4,0),(2,0), merge radius 2: two groups | One transitive connected component. Traversing only the upper-triangular adjacency drops backward connections. test-geometry.R |
| MLV-M007 | computeJunctionAnglesTypes.m | Diagonal bar (0,0)-(10,10), perpendicular stem (5,5)-(10,0): T minimum angle 0 degrees | Minimum angle 90 degrees. Interior arms use a different y convention from endpoint arms. test-geometry.R |
| MLV-M008 | lineIntersection.m | Crossing negative-coordinate segments around (-3,-3): empty result | Return the valid intersection (-3,-3); geometry need not be confined to the image. test-geometry.R |
| MLV-M009 | fitLineSegments.m | Two-point straight line | Returns an undefined-lineSegs error; a single span and defined scores are required. test-kernels.R |

## Other MATLAB findings

| ID | Status | Source and evidence | R treatment |
|---|---|---|---|
| MLV-M004 | Source + executed arithmetic; output sensitivity not established broadly | computeMAT.m uses .3 when threshold is omitted, but explicit 28 converts to .29887487943380647 | Omitted and explicit threshold_angle=28 use the same conversion; retain parameters |
| MLV-M010 | Executed Octave portability failure, not a demonstrated MATLAB parser bug | computeGradientVectorField.m has two jj++m_Neighbors8 expressions rejected by Octave | Reference harness uses a temporary, SHA-bound explicit-addition adaptation outside vendor |
| MLV-M011 | Source-confirmed missingness defect | MATpropertiesToContours.m marks every score equal to zero as invalid and replaces it with NaN | Carry a mapping mask; retain valid zeros; test-kernels.R |
| MLV-M012 | Contract discrepancy requiring scientific clarification | splitLDbyProperties.m describes mean curvature but ranks by the unnormalized length-weighted sum | Expose curvature (weighted mean) and integrated_curvature as distinct criteria |
| MLV-M013 | Source-confirmed global-state side effect | lineIntersection.m disables MATLAB:divideByZero warnings without restoring their state | Native predicates handle parallelism without changing R warning options |
| MLV-M014 | Documentation discrepancy | traceLineDrawingFromRGB.m help says structured-edge default; arguments block selects SAM | RGB extraction provider is explicit and recorded; no implicit model download |
| MLV-M015 | Source-confirmed undefined degenerate behavior; valid-domain policy needed | computeCurvature.m divides turning angle by zero segment length on a multi-segment contour | Constructor rejects zero-length geometry or explicitly drops it, retaining contour IDs |

## Python findings

Run `PYTHONDONTWRITEBYTECODE=1 python3 tools/upstream/reproduce_python.py`.
Executed output is retained in `docs/evidence/upstream-python.log`.
Python 3.14.7 and NumPy 2.4.3 were used. cv2/skimage were unavailable in that
interpreter; the full image pipeline was not executed.

| ID | Status | Source and minimal evidence | R treatment |
|---|---|---|---|
| MLV-P001 | Executed function | computeOrientation.py truncates endpoints (1.1,1.1)-(1.9,1.9) to int32 and returns 0 degrees | Preserve subpixel doubles; analytical result 315 degrees |
| MLV-P002 | Executed KeyError + additional source defects | getHorizontalVerticalStats.py reads orientation instead of orientations; it also applies radian trig to degree values and builds one fewer bin edge than its loop expects | Use one consistent schema, degree conversion, explicit complete boundaries and conservation tests |
| MLV-P003 | Source/module inventory | Branch functions import absent smoothData, diff, fitLineSegments, GetConSeg modules. traceSkeleton.py also returns R/F fields while downstream expects Radius/AOF | Native graph and explicit branch schema; no missing helper dependency |
| MLV-P004 | Source-confirmed alias + executed NumPy semantics witness, not full-pipeline execution | extract2DSkeletonFromBinaryImage.py sets skeletonImage=fluxImage then thresholds in place, so returned AOF is overwritten | Separate flux/skeleton buffers; test input/result immutability |
| MLV-P005 | Source-confirmed hidden external state and side effects | computeGradientVectorField.py replaces computed boundaries with outerBoundaryOriginal.mat and opens plots; computeMAT.py unconditionally writes MAT.mat | Numerical APIs require no working-directory files, plot windows or unsolicited output files |
| MLV-P006 | Source-confirmed incomplete taper and changed formula | computeMATpropertyPerBranch.py uses newddR before assignment; parallelism uses cumulative sums rather than the MATLAB local derivative magnitudes; taper uses dR where a second derivative is described | Explicit branch formulas tested on constant-radius and tapered branches |
| MLV-P007 | Source-confirmed array handling | computeMAT.py unconditionally selects binaryImage[:,:,0], rejecting a two-dimensional input described by its API | Binary matrices are an explicit supported input |
| MLV-P008 | Source-confirmed unusable projection state | MATpropertiesToContours.py initializes an ndarray then calls append; loop uses numContours inconsistently, and coordinate bookkeeping stores only row indices before unravel_index | Typed contour/score structures, exact indices and explicit mapping masks |
| MLV-P009 | Source-confirmed substituted cleanup parameter | extract2DSkeletonFromBinaryImage.py replaces calculated area_threshold with 100; the thinned array and returned skeleton are also handled inconsistently | Named, recorded component threshold and an explicitly thinned returned skeleton |
| MLV-P010 | Source-confirmed padding/smoothing discrepancy | getBranchDerivative.py leaves dR shorter than dX/dY for nontrivial branches and uses a radius moving-mean width of three versus MATLAB's two | Equal-sized padded derivatives; radius window two, with explicit endpoint convention |
| MLV-P011 | Shared empty-result defect plus source count semantics | getJunctionStats.py repeats the five-versus-eight empty histogram mismatch, and passes centers as np.histogram edges on nonempty input | Explicit event histogram edges and stable requested dimensions |

## Intentional changes that are not reported as upstream bugs

- R uses an explicitly named Octave-image Pratt morphology backend. Its LUT
  output matches the executed Octave fixtures; MATLAB skeletonization equivalence
  is not claimed.
- Exact EDT ties select the lowest R linear index. Other equally near choices
  can change AOF values. The reference report measures such differences.
- Contour tracing uses an explicit pixel graph without diagonal shortcuts across
  an existing orthogonal connection. It is not a claim to reproduce Kovesi's
  spur cleanup and contour ordering.
- Width-one rasterization uses clipped Bresenham lines. MATLAB insertShape
  equivalence is not claimed.
- Taper uses explicit moving means instead of unspecified smoothdata defaults.
- average_property(...,"orientation") is an explicitly documented axial circular
  mean, with NA for an isotropic resultant; it is not the original ad hoc
  concatenate-and-flip vector summary.
- Per-contour junction participation counts a junction once per contour.
- SVG curve_step controls sampling, not a certified geometric error bound.

Before filing a candidate upstream, rerun the minimal reproducer at the
maintainer's current revision and, where the observation is Octave-specific,
confirm whether it also occurs in MATLAB. Attach source SHA, runtime,
minimal input, observed/expected result, and the scientific consequence.
Do not turn source-only concerns into claims of executed full-pipeline failures.

## Findings from the expanded parity suite (0.1.1)

The fixtures in `inst/fixtures/parity` and `tools/parity/generate.py` execute
the pinned sources on curved, integer and subpixel contours, several bin counts,
varying-radius branches and six image topologies. Temporary adapters are recorded
in the fixture manifest. `tools/parity/compare.R` retains differences explicitly.

| ID | Status and trigger | Finding and R treatment |
|---|---|---|
| MLV-M016 | Executed terminal-bend branch; points (1,0) through (20,0), then (21,100) | fitLineSegments forces the last one/two points into the final span even when their mean perpendicular error exceeds its threshold of 1. The returned all-point span has mean interior error about 9.81. This also changes mirror scores on curved branches. R stops before the failing endpoint and tests every returned chord against an independent projection calculation. |
| MLV-M017 | Executed 13-bin axial histogram, orientations 90 and 270 degrees | getOrientationStats shifts, wraps and subtracts a half-bin; rounding puts axial-equivalent vertical directions into different bins. R scales the axial angle directly and assigns an exact boundary to its upper bin. Agreement is retained on ordinary 5/8-bin fixtures; the 13-bin boundary divergence is intentional and tested analytically. |
| MLV-P012 | Executed AOF on spatially varying nearest labels, radii 1 and 2 | computeAOF indexes IDX at int(px-1), int(py-1) after constructing zero-based sampling coordinates; this selects the preceding grid cell. It also omits epsilon from the sampling displacement. With the same distance/label fields, flux errors exceed 45 at radius 1. R agrees with the MATLAB source within 1e-10 on both radii. |
| MLV-P013 | Executed 37-degree rotation | rotateLinedrawing rounds sine and cosine to five decimal places. Synthetic coordinate errors reach approximately 0.000895 pixels and the transform is no longer exactly length-preserving. R uses full precision; its counterclockwise image-coordinate angle corresponds to the negative upstream angle. |
| MLV-P014 | Executed raw getLengthStats after computeLength; NumPy 2.4.3 | computeLength produces contourLengths with shape (n,1), but getLengthStats adds its one-element arrays into scalar histogram cells. Current NumPy raises ValueError: setting an array element with a sequence. The harness records the raw failure and separately labels an adapted comparison that reshapes only this field to (n,). |
| MLV-P015 | Executed 13-bin histograms; shared formula with MLV-M017 | Python retains the axial boundary instability of the MATLAB half-bin shift/subtraction formula. R's explicit upper-bin rule applies consistently after reversing segment direction. |

The retained Forest_DBW demonstration contains one zero-length contour. The
optional natural-drawing harness removes that one segment/contour explicitly
and identically for all three implementations, retaining 562 contours and
2,072 segments. This is a real-data witness for the already recorded MLV-M015
degenerate-input issue; no source snapshot or canonical data file is changed.

## R-port defects found and fixed in 0.1.1

| ID | Executed failure | Fix and regression |
|---|---|---|
| MLV-R001 | Sparse/isolated masks: the original R port assigned distance 0 instead of -1 to isolated interior ink pixels. Its outer-boundary predicate checked only for a white neighbor. | Require both an ink and white neighbor, matching the signed-distance construction. Six diverse masks now agree on distances within the single-precision reference tolerance; isolated pixels have analytical distance -1. |
| MLV-R002 | The original R port placed an exact vertical 13-bin boundary in the lower bin because its half-bin shift rounded downward. | Use the axial angle scaled by the bin count before rounding to a centered bin. Tests cover reversed vertical directions and odd bin counts 1, 3, 5, 7, 13 and 31. |

## Observed method differences, not numerical parity

Equidistant nearest-boundary pixels need not have the same label in R and Octave.
The R EDT selects the lowest column-major index. All differing labels in the
six-mask comparison have the same exact squared distance after MLV-R001 is fixed,
but their directions can change AOF and the final skeleton. Tests therefore
establish full distance parity, AOF parity with identical distance/label fields,
and morphology parity with identical flux fields. They do not label the differing
default end-to-end skeletons as matches. See the complete parity report for counts
and images. Actual MATLAB's distance/morphology tie behavior remains unexecuted.
