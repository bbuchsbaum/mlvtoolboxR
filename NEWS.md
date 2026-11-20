# mlvtoolbox 0.3.0

## Qualified against stored MATLAB outputs

* The package was compared with outputs the original lab computed in MATLAB
  for 3,017 line drawings (2,757 with medial-axis scores). Published-compatible
  options now reproduce them: `detect_junctions(rule = "mlv2023")` matches
  395,125 of 395,134 junctions exactly, and `compute_mat(aof = "mlv2023",
  morphology = "hmt8", branch_rule = "mlv2023", threshold_angle =
  asin(0.3 * pi / 2) * 180 / pi)` reproduces Forest_DBW's distance map, flux and
  skeleton exactly and 99.1-99.4% of all stored per-pixel symmetry scores within
  1e-4; with `estimator = "mlv2023"` all 2,284,239 stored segment orientations,
  lengths and curvatures agree within 1e-9. See `docs/PARITY.md`.

## Line drawings from images

* New `extract_lines()` turns a photograph or an image of a painting into a
  binary ink mask: CIE L* lightness, an anisotropic Kuwahara pre-filter
  against brush texture, the structure-tensor edge tangent flow, a flow-based
  difference of Gaussians integrated along the flow, and the XDoG threshold
  (Kyprianidis, Kang and Döllner 2009; Kyprianidis and Döllner 2008; Kang, Lee
  and Chui 2007; Winnemöller, Kyprianidis and Olsen 2012). Two further steps
  are this package's own: lines are placed on edges rather than in the
  dark-side XDoG band, and thin dark or light strokes are drawn once, as their
  centre lines. The default pre-filter removes dark strokes thinner than about
  2 pixels that do not lie along a pixel row or column, so use `smooth = 0`
  for thin line art. It is deterministic and needs no model weights or new
  dependencies. Without the pre-filter, quarter turns and reflections of the
  image transform every stage exactly; with it, the intermediate images agree
  to rounding and the masks have agreed exactly in every test. The mask
  records its definition, `"fdog-xdog-1"`, and every parameter in an
  `"extraction"` attribute.
* `xdog_lines()` is the same extractor shaped as a `trace_rgb()` backend, and
  `line_extraction_stages()` returns the intermediate images.
* `trace_rgb()` records a backend's `"extraction"` attribute, when present, as
  `source$provider_version` and `source$extraction_parameters`.
* `feature_spec(extraction = list(lines = list(...)))` fixes the line
  extractor's parameters; the specification records `definitions["lines"]`,
  `trace_rgb()` passes the parameters to the backend and checks its output,
  and specifications naming another version of the extractor are rejected. With
  any specification, `trace_rgb()` refuses a backend whose `"extraction"`
  attribute names a version other than the specification's `provider_version`.
  Specifications without line settings are unchanged.
* New vignette `line-extraction` ("From an image to a line drawing"), with a
  public-domain (CC0) painting from The Cleveland Museum of Art bundled in
  `inst/extdata`.

## Behaviour changes (new default methods)

The defaults now follow the alternative definitions recorded as MLV-M018-M032
in `DIFFERENCES.md` (for example, junctions and branch scores that do not
depend on contour direction, and curvature that counts each vertex turn once);
the `mlv2023` options keep the published behaviour. Feature specifications record the new method definitions,
so specifications and cached MATs from 0.2.x are rejected and must be rebuilt.

* Curvature: the default estimator is now vertex turning (`estimator =
  "turning"`) in `feature_spec()` and in every specification-free call
  (`as.data.frame()`, `compute_contour_properties()`, `contour_statistics()`,
  `feature_table()`, `predict_contours()`, `average_property()`, plotting and
  `feature_density()`). Each vertex turn is spread once over its dual cell, the
  half-segments on either side of the vertex. Segment tables keep one row per
  input segment for every estimator; with `"turning"` and `"arclength"` a
  segment's curvature is its mean curvature (the exact integral of the density
  over the segment divided by its length), with a signed mean in
  `signed_curvature`; `curvature_profile()` returns the piecewise density, and
  histograms integrate it exactly. The published estimator, now `estimator = "mlv2023"`, gives each
  segment the whole turn to the next segment divided by its own length and
  repeats the previous turn on a contour's last segment: the last vertex of every
  contour is counted twice (a 90 degree apex integrates to 180 degrees), the
  closing vertex of a closed contour is not counted, and values depend on
  sampling (MLV-M032). Turning has no tuning parameter and gives sampling-stable
  histograms: a circle of radius 10 drawn with 16, 32, 64 or 128 vertices,
  regularly or with alternating 1 : 2 spacing, puts all its length in one
  curvature bin (4.43-8.54 deg/px of 8), whereas `mlv2023` splits the
  alternating circles between 2.09-4.43 and 8.54-15.76. `"arclength"` remains
  for curvature at a stated spatial scale. Specifications state their
  estimator, so a specification built when `mlv2023` was the default keeps
  computing `mlv2023` values; specification-free feature schemas now record the
  curvature definition in effect (`vertex-dual-1` or `mlv2023-port-0.1.1`), so
  models bound to 0.2.x specification-free tables are rejected unless
  `estimator = "mlv2023"` is given.
  `compute_contour_properties(estimator = )` records its choice on the drawing
  (`curvature_estimator`), and later calls on that drawing that give neither
  `spec` nor `estimator` use it; geometry operations keep the choice and drop
  derived measurements. The turning estimator accepts exact
  retracing (an unsigned 180 degree turn, signed value `NA`), which occurs in 178
  of the 3,017 stored MLV drawings; `contour_turning()` and the arc-length
  estimator still reject it. Specification-free turning histograms are computed
  in the fused native pass (12,000 segments: 2.1 ms vs 1.9 ms for
  `"mlv2023"`).
* Histograms: a value within 1e-9 × max(1, |edge|) of a bin edge, in the
  coordinate in which the bins are equally spaced (degrees, log10(1 + x) or the
  raw value), is binned as if exactly on the edge, so the upper-bin rule applies
  and mathematically equal values share a bin on every platform (for example a
  computed junction gap of 59.999999999999943 and an exact 60 degrees). R and the
  fused native kernel apply the same rule with the same arithmetic. Rounding
  residue is about 1e-13 relative or less; the tolerance is far below any real
  measurement difference. Snapping is a numerical robustness fix that acts
  only on values within that tolerance of an edge; it does not change the
  histogram definition (`definitions["histogram"]` stays `"1"`) or feature
  schema compatibility. On the stored MLV drawings this changes the bin of
  1,109 of 389,974 junction smallest gaps (8 bins) and of no segment orientation
  or curvature value.

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
* Default branch tracing (`branch_rule = "invariant"`, definition
  `direction-symmetric-tiemax-graph-2`) no longer deletes small shapes.
  `min_branch_points` now removes only spurs (paths from an end pixel to a
  junction), shortest first, joining paths through every pixel where nothing
  else branches; tied spurs are kept together when removing them would leave
  their junction with fewer than two paths, and no skeleton component loses
  its last path. A 20 px square previously had no branch and no mapped score;
  its four arms are now kept and all 80 contour pixels are mapped. The result
  does not depend on scan order, start pixel or lattice orientation.
  Specifications and cached MATs with the `-graph-1` definition are rejected.
* Default branch scores are NA where their estimator is undefined, not 0
  (`branch_rule = "invariant"`, definition
  `direction-symmetric-tiemax-graph-defined-3`; MLV-M033). The published
  estimator scores the end points of parallelism and taper (zero-width
  window), branches of fewer than three points, and mirror on three or fewer
  points as 0, and its width-3 moving mean averages that 0 into the next
  point. The default marks them NA and averages only defined values, so a
  branch end takes the value of its defined neighbour (before smoothing);
  traversal directions, score images and contour mapping combine defined
  values only. The default also sums each smoothing window directly: the
  published cumulative-sum differences can round a radius pair (1, 1) below 1
  and turn an exact separation of 0 into NA. A raw score of exactly 0 is a
  defined value and is never replaced by NA. `score_branch()` and
  `branch_rule = "mlv2023"` keep the published values. Compared with the
  published estimator, scores change at undefined points and the points
  whose smoothing window contains one, and, for separation, around radius
  pairs that the published sums round below 1; elsewhere they agree to
  rounding (at most 3.1e-13 on twelve sampled Forest drawings). On the 2,757
  stored drawings per-contour parallelism correlates 0.811 with the published
  values (0.926 before), whose branch ends are lowered by those zeros;
  separation and mirror are unchanged (0.982, 0.864). Specifications and
  cached MATs with the `-graph-2` definition are rejected.
* `split_contours()` no longer fails when a contour has no finite medial-axis
  score (more common now that undefined scores are NA): such contours cannot
  be ranked, so they are left unassigned, scored NA and listed in the new
  `unranked` element. Ranks are taken among rankable contours; the length
  budget is still a fraction of the total length of all contours.
* Skeleton pixels lying on ink are excluded before branch tracing.
* `mlv_drawing()` treats segments whose endpoints differ only by rounding residue
  as zero-length.

## Documentation and plotting

* New vignettes: `getting-started` (from a line drawing to contour features),
  `junctions`, `medial-axis` and `mlv-compatibility`; `workflows` now covers file
  input, extraction and feature recipes.
* New `plot.mlv_mat()` shows each medial-axis stage (ink, distance, flux,
  skeleton, branches, scores) in image coordinates with labelled legends.
* `plot.mlv_drawing()` draws colour bars with units, uses a cyclic palette for
  orientation, can draw junction arms, frames the canvas on the raster extent
  so drawings and MAT stages overlay exactly, and draws mapped scores as
  filled pixels.
* `detect_junctions()` returns an `arms` list-column with each junction's arm
  directions; `angles` are the gaps between consecutive arms.

## Revisions in this port

* Junction arm angles no longer depend on the merged junction position, and
  gapped T-junctions no longer gain an additional arm (MLV-R005-R008, R014).
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
