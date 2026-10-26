# Local performance qualification

Current follow-up: [0.1.1 performance](PERFORMANCE_0.1.1.md) and [expanded parity](PARITY.md). This page retains the 0.1.0 measurements.

The R port is faster on the measured, correctness-checked stages below.
These are synthetic local workloads, not an end-to-end natural-image study or
a timing claim about MATLAB itself.

## Environment and revisions

- Apple M3 Max, arm64, 36 GiB RAM; macOS Sonoma 14.3.
- R 4.5.1; installed C++17 kernels built with Apple clang 15.0.0, `-O2`.
- Octave 11.1.0, image 2.20.0, running MATLAB source
  `f12418a5e03c2f38e222fed1dffccf42b227cd29`.
- Python 3.14.7, NumPy 2.4.3, SciPy 1.17.1, source
  `14b3183b91146b1021eb7ef97db75ad703b06832`.
- R before optimization: `21a3b5d8ad076d39cb26e5c08f7dcafa6213aebe`;
  after: `b26cb52340447cb3faf9d15773b029ee2768504a`.
- Final package qualification adds tests at
  `141a4a44dcc6431faa53a6a4b8af31b06b306982`; its R/C++ code is identical
  to the benchmarked candidate.

Runs were sequential, with no overlapping builds or tests from this session.
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` were set to 1
for the final R, Python and Octave runs. The kernels measured are serial.
The earlier R baseline did not explicitly set those variables. This was a
shared desktop, without CPU isolation, so min/max variability is retained.

## Python comparison

Warm medians in milliseconds; speedup is vendored Python time divided by R time.

| Workload | Python ms | R ms | Speedup |
|---|---:|---:|---:|
| Geometry, 100 one-segment contours | 0.860 | 0.410 | 2.10× |
| Geometry, 1,000 one-segment contours | 8.554 | 1.336 | 6.40× |
| Geometry, one 10,000-segment contour | 5.683 | 1.094 | 5.20× |
| EDT, 32×39 | 0.692 | 0.0146 | 47.2× |
| EDT, 64×71 | 2.471 | 0.0498 | 49.6× |
| EDT, 128×135 | 9.182 | 0.193 | 47.5× |
| EDT, 512×519 | 146.639 | 3.297 | 44.5× |

Geometry computes length, directed orientation and curvature in both libraries.
Straight integer-coordinate inputs have analytically known unit segment lengths
and zero orientations/curvatures. Both outputs are checked before timing.
These inputs deliberately avoid the demonstrated Python subpixel truncation
defect: incorrect outputs are not credited as a fast result.

The rectangular EDT uses one feature pixel. Every distance is checked against
the analytical Euclidean distance (absolute error below 1e-12), and every nearest
label identifies that physical feature, accounting for row-major versus
column-major indexing. This benchmark does not qualify tie behavior.

Raw results: [Python](evidence/benchmark-python.csv),
[R](evidence/benchmark-r-after.csv).

## MATLAB source through Octave

Warm medians in milliseconds. The R MAT call defers branch tracing because the
upstream `computeMAT` returns distance, AOF and skeleton without tracing.
Both sides use 60 samples, radius 1, and an explicit 28-degree threshold.

| Workload | Octave ms | R ms | Speedup |
|---|---:|---:|---:|
| AOF, 32×32 | 1,064.1 | 0.344 | 3,095× |
| AOF, 64×64 | 4,760.5 | 1.453 | 3,276× |
| MAT core, 32×32 | 1,077.5 | 0.543 | 1,984× |
| MAT core, 64×64 | 5,001.2 | 2.094 | 2,389× |

The upstream AOF is an interpreted nested loop. The R version uses compiled
Rcpp loops; these ratios must not be extrapolated to MATLAB's JIT runtime.

AOF inputs use a point-distance field; MAT inputs use two vertical ink lines.
Before timing R, the harness compares its outputs with retained Octave arrays:
AOF absolute differences must be below 1e-10, signed-distance differences below
2e-6 (the Octave distance backend uses single precision), and skeletons must be
pixel-identical. All checks passed at both sizes.

The reference harness makes exactly two temporary syntax substitutions,
`jj++m_Neighbors8` → `jj + m_Neighbors8`, so Octave can parse the upstream
gradient code. The pinned vendor snapshot remains unchanged.
The default/explicit threshold discrepancy is recorded in the bug register;
both benchmark paths receive the explicit threshold.

Raw results: [Octave](evidence/benchmark-octave.csv),
[R](evidence/benchmark-r.csv). Corresponding output arrays are retained beside them.

## What changed in R

Rprof attributed 95.9% of the original contour workload to segment-table
construction, including repeated data frames and row binding. The optimized
wrapper flattens contours once, vectorizes segment calculations, and uses
contiguous groups instead of repeatedly scanning the whole table. Junction
criteria compute junctions once per selection operation. Logical image masks
also avoid repeated numeric conversion and membership scans.

| Workload | Before ms | After ms | Improvement |
|---|---:|---:|---:|
| 100 one-segment contours | 25.125 | 0.410 | 61.3× |
| 1,000 one-segment contours | 247.000 | 1.336 | 184.9× |
| One 10,000-segment contour | 1.797 | 1.094 | 1.64× |
| EDT, 512×519 | 6.594 | 3.297 | 2.00× |

A saved mixed-contour, subpixel witness from the checked pre-optimization package
is compared with `expect_identical()`, including complete property tables.
Thus the wrapper speedup does not rely on relaxed floating-point tolerances.
The compiled EDT, AOF, Pratt morphology, graph tracing, intersection sweep,
nearest-neighbor search and line-fitting kernels were already in Rcpp.
Profiling did not justify adding another native implementation for the wrapper.

Raw [before](evidence/benchmark-r-before.csv) and
[after](evidence/benchmark-r-after.csv) timings and
[original profile](evidence/geometry-before.Rprof.summary.txt) are retained.
The short post-optimization profile has too few samples for stable percentages.

## Scaling and measurement limits

R-only warm medians include AOF at 512×512 (102.5 ms), MAT core at 512×512
(139 ms), and junction detection on 10,000 spatially separated segments
(10.19 ms). The junction case is sparse; it does not establish performance on
dense all-intersecting arrangements.

R/Python timings use five warmed batches, selecting repetitions until a batch
takes at least 0.12 seconds or reaches 1,024 calls. Octave uses three warmed
individual calls. CSVs retain repetitions, minima and maxima. Their
`cold_seconds` field means the first timed call after setup and output
validation; it is not a fresh-process/startup measurement. R's first-call
clock is millisecond-resolution, so some small entries are zero.
Input generation, imports, model loading, disk I/O and validation are excluded
from the timed region. Peak RSS and parallel scaling were not measured.

The full Python image pipeline is incomplete at this pinned revision; missing
helpers, schema defects and side effects are recorded in
[UPSTREAM_BUGS](../UPSTREAM_BUGS.md). No timing win is claimed for a failing
Python pipeline. Rasterization and graph-tracing method differences also
preclude unqualified end-to-end equivalence claims.

## Reproduce

From this checkout, with the documented runtimes and optional R dependencies:

```sh
python3 tools/bounded_check.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 tools/benchmark_python.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 tools/benchmark_octave.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla tools/benchmark.R .agent-work/package-build/library docs/evidence/benchmark-r.csv
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla tools/benchmark_geometry.R .agent-work/package-build/library docs/evidence/benchmark-r-after.csv .agent-work/geometry-after.Rprof
```

Use a UTF-8 process locale (here `LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8`).
The bounded builder prestages only package files, avoiding the large reference
data tree. It monitors its owned workspace against the approved 200 MiB cap.
