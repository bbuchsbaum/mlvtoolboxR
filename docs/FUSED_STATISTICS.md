# Integrated histogram performance, 0.2.1

The earlier Rcpp prototype is now integrated into the package. A fused native
pass computes the published geometric histograms directly from contours,
avoiding a full segment data frame and grouped R lists. On this machine the
measured geometric workloads are **3.83–9.35× faster than the qualified 0.2.0
implementation**, and **1.33–1.89× faster than the earlier grouping/length
prototype adapted to the same feature-schema API**.

This qualification covers an immutable 388-file source snapshot, including
the completed feature-schema, estimator and designer additions. Subsequent
graph work in the shared checkout is outside this snapshot. The check ran
before the combined Git checkpoint existed; the
[source manifest](evidence/fused-statistics/package-source-manifest.json) binds
the actual tested bytes. All 388 package files at local source commit
`bf9f0796010c18c58168f518cdbbc35ca37a5a1c` match that manifest exactly
([checkpoint receipt](evidence/fused-statistics/checkpoint.json)). The retained archive is
[mlvtoolbox_0.2.1-fused.tar.gz](../mlvtoolbox_0.2.1-fused.tar.gz), 344,689 bytes,
SHA-256 `36adfe06494da8a2e6a07ffc6c325007374ab14225356bd6889695e19885fee1`.
The separate designer archive and its earlier evidence remain unchanged.

## Measurements

Milliseconds per public `feature_table()` call, including output and feature
metadata construction. The default geometric cases request orientation,
length and curvature, eight bins, normalized per contour.

| Workload | R 0.2.0 | Earlier prototype on 0.2.0 | Integrated 0.2.1 | Gain over 0.2.0 | Gain over prototype |
|---|---:|---:|---:|---:|---:|
| Forest: 562 contours / 2,072 segments | 5.875 | 1.852 | 1.391 | 4.22× | 1.33× |
| 1,000 contours × 20 segments | 18.125 | 3.656 | 1.938 | 9.35× | 1.89× |
| One contour × 20,000 segments | 8.000 | 2.625 | 1.539 | 5.20× | 1.71× |
| Forest with physical units and fixed feature specification | 6.219 | 2.391 | 1.625 | 3.83× | 1.47× |
| Arc-length curvature: 8 contours × 20 segments | 8.938 | 8.562 | 8.875 | 1.01× | 0.96× |

The table uses the retained second baseline/prototype run and final qualified
candidate. Each entry is the median of five warmed batch averages; adaptive
batch size targets at least 0.12 seconds. Processes ran sequentially with BLAS
and OpenMP thread variables set to one. This was not a globally isolated
machine. Earlier repeated runs and the superseded pre-stride archive are retained,
including less favorable timings;
the gains are not inferred from a single call. Reading inputs, process startup,
compilation, rendering and the image/MAT pipeline are outside these timings.

Fresh adapted Python forest timings have a median of **5.241 ms**
(three round medians: 5.239–5.406 ms), or **3.77×** the final R time.
Python's raw length histogram still fails under NumPy 2.4.3; the timed path
reshapes contour lengths from `(n,1)` to `(n,)` before timing. Its maximum
difference from the exported R histogram is **5.46e-12**, within the existing
1e-8 tolerance. The Python comparison used an intermediate candidate whose
complete output is exactly identical to the final one.

The older isolated 0.1.1 prototype reported 5.44× Python. That prototype did
not construct the new feature schema and measurement metadata, so its absolute
timing is not a comparison of the current public API. The prototype column above
preserves the 0.2.0 metadata work and provides the relevant comparison for this
integration. Arc-length curvature retains its exact affine integration path;
these measurements show essentially no improvement for it.

[Summary CSV](evidence/fused-statistics/performance-summary.csv),
[all retained R runs](evidence/fused-statistics/benchmarks),
[Python timings](evidence/fused-statistics/python-timings.csv),
[Python adapter and validation](evidence/fused-statistics/python-validation.json).

## Allocation changes

Separate `Rprofmem()` measurements average five calls. These count bytes
allocated through R, not peak memory or native `std::vector` allocations.

| Workload | Baseline R bytes/call | Final R bytes/call | Reduction |
|---|---:|---:|---:|
| Forest | 1,710,584 | 627,144 | 63.3% |
| 1,000 contours | 9,640,208 | 1,110,696 | 88.5% |
| One long contour | 7,650,960 | 3,696 | 99.95% |
| Physical forest | 1,725,608 | 631,328 | 63.4% |
| Arc-length | 1,573,376 | 1,417,256 | 9.9% |

The native path keeps scratch lengths/orientations for one contour plus the
histogram accumulators. Thus the long-contour row does **not** mean its total
memory use fell to 3.7 KB. No peak-RSS or native-allocation claim is made.

## Implementation and numerical contract

[Native kernels](../src/statistics.cpp) provide grouped accumulation, contour
lengths and fused published geometric histograms.
[R dispatch](../R/statistics.R) retains bin definitions, normalization,
calibration, schema construction and missingness behavior. Turning and
arc-length estimators continue through their existing measurement paths.
Junction/MAT properties retain their existing semantics. Geometry is read from
the current drawing on every call; there is no persistent cache to become stale.

The kernels preserve input order within each sum, extended accumulation where
the platform provides it, squared-term rounding without contraction, axial
wrapping and upper-edge clamping. The arithmetic choices were checked against
R's [sum implementation](https://raw.githubusercontent.com/wch/r-source/R-4-5-branch/src/main/summary.c),
[arithmetic implementation](https://raw.githubusercontent.com/wch/r-source/R-4-5-branch/src/main/arithmetic.c),
and [row-sum implementation](https://raw.githubusercontent.com/wch/r-source/R-4-5-branch/src/main/array.c).
Native shape/index guards and finite positive segment-length checks reject
malformed manually edited drawings before bin indexing. Two pre-release safety
hazards found during source review are recorded separately from upstream
findings in [the bug ledger](../UPSTREAM_BUGS.md).

## Verification

- **1,224 passing assertions**, no failures, warnings or skips, including
  **210 focused native/statistics assertions**.
- **122 saved public-output witnesses** generated from an independently
  installed qualified 0.2.0 archive. Full outputs and metadata are compared
  exactly. Cases include physical calibration, empty drawings, mixed MAT
  missingness, repeated properties and all three estimator modes.
- Randomized independent grouped-sum references, cancellation, missing values,
  overflow, contour lengths across scales, axis/subnormal boundaries,
  conservation, edited geometry and invalid native inputs.
- All five final benchmark outputs and attributes are exactly identical to
  0.2.0: [receipt](evidence/fused-statistics/exact-outputs-qualified.csv).
- **94.6867% R/C++ line coverage**, including instrumented native source
  ([coverage](evidence/fused-statistics/coverage.json)).
- Full `R CMD check --as-cran --no-manual`: **0 errors, 0 warnings, 2 notes**.
  Notes are “New submission” and inability to verify the current time.
  Examples and rebuilt vignettes passed.
  [Check log](evidence/fused-statistics/00check.log),
  [test output](evidence/fused-statistics/testthat.Rout).
- Final installed-package comparison against retained executed reference
  fixtures: **151 Octave matches / 22 known differences**, and
  **75 Python matches / 41 known differences**. No unexplained comparison
  failures were introduced. [All rows](evidence/fused-statistics/parity-qualified.csv).

The 63 differences remain visible; they are not relabeled as parity.
The [parity report](PARITY.md) explains nearest-label ties and upstream
formula/binning defects. Default MAT skeletons differ on 4/6 masks because
equidistant boundary labels can differ; aligning those labels gives all six
matches. Actual MATLAB was not executed. This qualification is macOS arm64,
R 4.5.1; Linux, Windows and platforms with extended `long double` still need
their own executed evidence.

## Reproduction and provenance

Install the retained archives into separate temporary libraries, then run:

```sh
Rscript --vanilla tools/performance/fused_statistics.R BASELINE_LIB OUT_BASE baseline-0.2.0
Rscript --vanilla tools/performance/fused_statistics.R BASELINE_LIB OUT_PROTO grouping-prototype prototype
Rscript --vanilla tools/performance/fused_statistics.R FINAL_LIB OUT_FINAL fused-0.2.1-qualified
python3 tools/performance/histogram_python.py OUT_FINAL OUT_PYTHON
Rscript --vanilla tools/parity/compare.R FINAL_LIB OUT_PARITY.csv
```

Run from the repository root with the reference-data link available, using
`LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1`. The baseline archive SHA-256 is
`7c2d7e20e95d7b61991bcfe09e8372878d90a140b59c49d55049d1fb51954954`.
The original forest file is unchanged; the same one zero-length segment/contour
is dropped for both backends. Recreate saved witnesses with
[statistics_reference.R](../tools/performance/statistics_reference.R).

The bounded-check helper now accepts `--source-dir DIR` and `--work-name NAME`
to stage a frozen input separately from another agent's work. The final full
check and coverage both used the same frozen source and the existing 200 MiB
working-directory limit. Their observed workspace peaks were 14.3 MiB for the
check and 20.4 MiB for coverage. This cap describes on-disk workspace size,
not process RAM.

[Full provenance](evidence/fused-statistics/provenance.json) retains pinned
source revisions, archive hashes, runtime, preparation, scope and limitations.
No upstream report, package publication or remote push was performed.
