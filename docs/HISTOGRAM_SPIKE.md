# Remaining histogram performance headroom

Integrated follow-up: [qualified fused statistics in 0.2.1](FUSED_STATISTICS.md).
The measurements below preserve the original 0.1.1 prototype experiment.

This is an isolated prototype, **not a new package release**. It answers whether
the earlier 1.12× Python histogram advantage is a practical limit. It is not:
two small Rcpp kernels bring the measured forest workload to about **5.4×
faster than adapted Python**, with exact preservation of tested R 0.1.1 outputs.

The live feature-schema implementation was concurrently editing
`R/statistics.R`. This experiment loaded the saved, checked 0.1.1 archive into
a private library and used separate function environments; it did not modify
the installed namespace or production sources. Integration with the evolving
feature schema and a full package check remain necessary before release.

## Fresh measurements

Median of three round medians, milliseconds. Each R round alternated method
order and used five warmed batches. All three histogram properties (orientation,
length, curvature), eight bins, normalized per contour.

| Workload | R 0.1.1 | Rcpp grouping only | Rcpp grouping + lengths | Gain over R 0.1.1 |
|---|---:|---:|---:|---:|
| Forest: 562 contours / 2,072 segments | 4.844 | 2.750 | 0.961 | 5.04× |
| 1,000 contours × 20 segments | 16.875 | 6.063 | 2.813 | 6.00× |
| One contour × 20,000 segments | 6.844 | 2.078 | 1.852 | 3.70× |

Fresh Python timing on the same forest drawing was **5.224 ms**, versus
**0.961 ms** for the combined prototype: **5.44×**. Python's three round medians
were 5.207–5.248 ms; prototype R medians were 0.953–0.984 ms.
The raw Python length-histogram function still fails under NumPy 2.4.3.
The timed path applies the previously documented `(n,1)` to `(n,)`
contour-length shape adapter before timing.

Fresh R 0.1.1 is only 1.08× faster than fresh adapted Python on this workload,
consistent with the earlier modest 1.12× result. This session's absolute times
differ from the earlier session, so the new ratios use fresh measurements.
The 1,000-contour case here uses three properties/eight bins, unlike the prior
four-property/16-bin benchmark; their absolute times should not be compared.

[R timings](evidence/histogram-spike/timings.csv),
[Python timings](evidence/histogram-spike/python-timings.csv),
[Python parity and raw failure](evidence/histogram-spike/python-validation.json).

## What the prototype changes

1. Replace `split(weights, key)` plus `vapply(sum)` with direct indexed
   accumulation in C++. Integer groups/bin indices already exist, so factor
   construction, sorting group labels and many small R objects are unnecessary.
2. Compute per-contour lengths in one compiled traversal, replacing hundreds
   of R calls to subtraction, squaring, `rowSums` and `sum`.
   Separate squared terms prevent multiply-add contraction from silently
   changing the established floating-point operations.

The R bin definitions, orientation wrapping, normalization, output construction
and scientific formulas remain unchanged. The accumulation retains each bin's
input order and uses `long double` (eight bytes on this tested arm64 platform).
Cross-platform equality must still be checked when integrating it.

[Native prototype](../tools/performance/histogram_candidates.cpp),
[R experiment](../tools/performance/histogram_spike.R),
[Python comparison](../tools/performance/histogram_python.py).

## Correctness evidence

The experiment passed **296 exact `identical()` comparisons** before timing:
saved 0.1.0/0.1.1 statistics and feature witnesses including MAT/junction
properties; 150 adversarial grouped-sum cases with missing weights, cancellations,
empty bins and varying magnitudes; all five integer/subpixel contour fixtures
across 1/5/8/13/31 bins, raw/normalized values and drawing/contour summaries;
and complete feature outputs plus contour lengths for all three timed drawings.

Fresh Python normalized histogram output differs by at most **5.46e-12** from
the R output exported to CSV, within the existing 1e-8 comparison tolerance.
The prototype is exactly equal to the baseline R result before CSV export.
The forest file is unchanged; the same one zero-length segment/contour is
explicitly dropped.

[Validation](evidence/histogram-spike/validation.json),
[full experiment log](evidence/histogram-spike/experiment.log),
[provenance](evidence/histogram-spike/provenance.json).
This is targeted prototype validation. It is not a replacement for integrated
package tests, instrumented C++ coverage or R CMD check.

## What is still on the table

The combined prototype's remaining profile includes segment-table construction
and repeated property/native-call setup. A fused geometry/bin accumulation path
could avoid building the full segment data frame when only histograms are
requested. That is the next plausible single-drawing optimization; no further
speedup is measured or promised here.

Explicit reuse of validated geometry may also help repeated feature requests,
but caching must respect mutation and the feature-schema work. Parallelism
belongs at the independent-drawing level for sufficiently large batches;
a roughly one-millisecond single drawing is too small to justify worker startup.

Other already measured stages were faster than Python before this prototype:
forest geometry about 5.8× and synthetic EDT about 44–50× in their respective
earlier runs. The 1.12× number never described the whole toolbox. See the
[0.1.1 report](PERFORMANCE_0.1.1.md) and [earlier stage report](BENCHMARKS.md).

## Reproduce and measurement limits

Apple M3 Max, 36 GiB RAM, macOS Sonoma 14.3, R 4.5.1, Apple clang 15,
C++17 `-O2`; Python 3.14.7 / NumPy 2.4.3.
Baseline archive SHA256:
`1860a64ac5ac92284972f2e9ba76b8f78e97064e3634dd1fedc79d42e11b668d`,
checked source `d83ed952fdc9f639b98e849b0d8697259894de46`.
Python reference: `14b3183b91146b1021eb7ef97db75ad703b06832`.

Run from the checkout after installing that archive into a private library:

```sh
export LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla tools/performance/histogram_spike.R PRIVATE_LIBRARY OUTPUT_DIR WORK_DIR
PYTHONDONTWRITEBYTECODE=1 python3 tools/performance/histogram_python.py WORK_DIR OUTPUT_DIR
```

The forest reference file is required at its canonical local path. The captured
invocations set a private compiler cache, use `-O2`, and exclude import,
compilation, input generation and file I/O from warmed timings. R and Python
were measured sequentially; Python rounds followed the R rounds, so they are
not cross-process alternating pairs. This is a shared desktop without CPU
isolation. The absolute difference is large relative to within-run variability,
but no cross-machine, peak-memory or parallel-scaling claim is made.

The prototype and receipts are retained. Its private installation, compiler
cache and temporary numeric data are removed after verification.
