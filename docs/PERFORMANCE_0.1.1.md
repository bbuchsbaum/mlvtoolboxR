# Performance follow-up: 0.1.1

Grouped histogram construction and MAT mapping are faster while preserving the
complete outputs of all five representative old/new workloads exactly.
The full installed test suite and the expanded [parity study](PARITY.md)
qualify the changes; performance does not override corrected behavior.

## Measured changes from 0.1.0

Warm medians on the same machine. The contour pair ran back-to-back; MAT used
alternating baseline/candidate runs with an unchanged raster control.

| Workload | R 0.1.0 | R 0.1.1 | Improvement |
|---|---:|---:|---:|
| Feature table, 1,000 contours × 20 segments, four properties | 1,609 ms | 30.875 ms | 52.1× |
| Long statistics table, same drawing/properties | 1,521 ms | 33.375 ms | 45.6× |
| Full MAT analysis, four elliptical rings on 256×256 canvas | 639–641 ms | 334 ms | 1.91–1.92× |

The feature baseline's five batches ranged from 1,558 to 1,618 ms; the candidate
ranged from 30.5 to 32.125 ms. Statistics ranged from 1,520 to 1,527 ms and
32.125 to 33.75 ms. Two final MAT runs each had median 334 ms. These are observed
local timings, not confidence intervals or guarantees for other datasets.

Raw [paired contour baseline](evidence/parity-perf/features-paired-before-timings.csv)
and [candidate](evidence/parity-perf/features-paired-final-timings.csv);
MAT [baseline 1](evidence/parity-perf/pair-before-1-timings.csv),
[candidate 1](evidence/parity-perf/final-timings.csv),
[baseline 2](evidence/parity-perf/pair-before-2-timings.csv),
[candidate 2](evidence/parity-perf/pair-final-2-timings.csv).

The unchanged raster control had medians 0.336–0.359 ms in those alternating
runs. Earlier unpaired runs were faster for both unchanged controls: junction
detection went from 109 to 153 ms and rasterization from 0.242 to 0.344 ms
between the initial baseline and final run. Their code and outputs are unchanged;
this is evidence of desktop timing variation. We retain those
[initial](evidence/parity-perf/before-timings.csv) and
[final](evidence/parity-perf/final-timings.csv) measurements.
The earlier unpaired feature ratio was about 35×, so 52× is specific to the
back-to-back pair, not a universal factor.

## Comparison with vendored libraries on the forest drawing

All implementations use the same explicitly cleaned 562-contour, 2,072-segment
drawing. Geometry and normalized eight-bin length/orientation/curvature
histograms pass the comparisons in [PARITY](PARITY.md#retained-natural-drawing).

| Stage | Octave / MATLAB source | Python | R 0.1.0 | R 0.1.1 |
|---|---:|---:|---:|---:|
| Geometry | 63.625 ms | 6.0605 ms | 1.0625 ms | 1.0469 ms |
| Per-contour histograms | 83.292 ms | 7.471 ms* | 619 ms | 6.6563 ms |

The new histogram time is 93.0× faster than R 0.1.0, 12.5× faster than the
MATLAB source interpreted by Octave, and 1.12× faster than adapted Python.
The small Python margin is a local observation on one drawing. Geometry was
not changed in this follow-up.

*Raw Python length histograms fail on their own `(n,1)` field under NumPy 2.4.3.
The timed adapted path reshapes that field to `(n,)` before timing.
The source numerical functions are unchanged; raw failure and adapted parity
are separately recorded. Histogram setup starts from computed geometry in all
three implementations. R currently recomputes the segment table internally.

[Reference timings](evidence/parity-perf/natural-timings-reference.csv),
[R baseline](evidence/parity-perf/natural-before-timings.csv),
[R candidate](evidence/parity-perf/natural-final-timings.csv).
These results do not estimate MATLAB JIT performance or qualify a complete
Python image pipeline. Earlier stage-specific Rcpp/EDT results remain in the
historical [0.1.0 benchmark report](BENCHMARKS.md).

## Why it improved, and what stayed exact

The initial feature profile attributed 90% of total time to
`contour_statistics`, including 47% in data-frame construction.
Statistics now accumulate grouped bin matrices and construct one final result;
`feature_table` directly constructs the wide output. Group sums retain R's
summation order/precision. Empty results, duplicate property requests,
normalization and per-contour junction participation retain explicit behavior.

The initial full MAT profile spent 67% in score mapping, including 37% in
nearest-neighbor searches. An expanded tangent bounding box discards only queries
that cannot lie within the requested radius; remaining queries still use exact
nearest-neighbor search and the strict radius boundary. Tangent generation
also computes invariant geometry once instead of repeating it for nine sample
offsets. Existing Rcpp kernels remain in use; this follow-up adds no new
compiled kernel because the measured overhead was in R orchestration.

[Initial feature profile](evidence/parity-perf/before-features_1000x20-profile.txt),
[initial MAT profile](evidence/parity-perf/before-mat_rings_256-profile.txt),
[final MAT profile](evidence/parity-perf/final-mat_rings_256-profile.txt).

The [exact-output receipt](evidence/parity-perf/exact-workflow-outputs.json)
compares complete inputs and outputs for features, statistics, a 19×19 junction
grid, full MAT and rasterization: all are `identical()`. Separate committed
regressions check 32 old statistics configurations, a feature-table witness,
six exact tangent witnesses, and mapping against exhaustive nearest-sample
search. Two intentional correctness fixes change outputs outside these
non-boundary benchmark domains: isolated ink signed distances and odd axial
histogram boundaries. They are documented as MLV-R001/R002.

## Environment, provenance and timing procedure

- Apple M3 Max, arm64, 36 GiB RAM; macOS Sonoma 14.3.
- R 4.5.1; existing C++17 kernels installed with Apple clang 15, `-O2`.
- Baseline archive 0.1.0, delivered at
  `0092e83c78dc9f357e781a355897d856747ef2cf`; runtime matches its checked source.
- Candidate archive 0.1.1, checked source
  `d83ed952fdc9f639b98e849b0d8697259894de46`.
- Octave 11.1.0/image 2.20.0, pinned MATLAB source
  `f12418a5e03c2f38e222fed1dffccf42b227cd29`.
- Python 3.14.7/NumPy 2.4.3, pinned source
  `14b3183b91146b1021eb7ef97db75ad703b06832`.

Benchmark processes were sequential, without overlapping tests/builds from this
session. OMP/OpenBLAS/MKL thread variables were set to one. There was no CPU
isolation on this shared desktop. Runtime/startup and input generation are excluded.
CSV `first_seconds` means the first timed call after loading/setup, not
fresh-process startup; the natural-drawing R harness also validates outputs
before timing. Millisecond-resolution first calls may round to zero.

R uses five warm batches, doubling repetitions until at least 0.15 seconds or
128 calls. Natural reference timings use one initial and three subsequent
individual calls. Output sizes are retained, but peak RSS, allocation totals,
parallel scaling and complete photograph-to-feature performance were not measured.
Do not multiply ratios from different workloads or releases.

## Reproduce

Build/install via `python3 tools/bounded_check.py`, then:

```sh
export LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla tools/benchmark_workflows.R .agent-work/package-build/library .agent-work/perf-new
```

To compare an installed baseline, run the same script against its library and
a different output directory. Optional arguments three/four select a workload
regular expression and record the installed package's exact source commit.
Then assert exact outputs:

```sh
Rscript --vanilla tools/parity/compare_workflows.R .agent-work/perf-old/outputs.rds .agent-work/perf-new/outputs.rds .agent-work/perf-exact.json
```

The optional natural-data commands are in [PARITY](PARITY.md#reproduce).
Retained [qualification](QUALIFICATION_0.1.1.md) identifies the checked archive,
source manifest and executed tests.
