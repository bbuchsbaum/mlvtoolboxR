# Implementation status

The local R package is version 0.1.0. This file supersedes the original
assessment's statements that no implementation exists.

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

- Independent analytical/metamorphic geometry tests and exhaustive small EDT/KNN
  witnesses, plus variable-field scalar AOF comparisons.
- Exact Pratt agreement on three executed Octave fixtures (box, ring, noise).
- On the 25-by-35 reference MAT fixture: skeletons match at every pixel
  (44 pixels); maximum flux discrepancy about 1.1e-14.
- Angle construction and one-based grid sampling order are tested explicitly:
  mathematically equivalent reassociation can select different cells near ties.
- R default and explicit 28-degree threshold share one rule. Degenerate masks,
  short branches and valid-zero mappings have defined behavior.

## Qualification

Source candidate 141a4a44dcc6431faa53a6a4b8af31b06b306982 passes the installed
250-assertion regression suite and R CMD check --as-cran --no-manual with
0 errors, 0 warnings and one informational New submission NOTE. Coverage is
92.6117% across R and instrumented C++, up from 82.8125%.

Matched local benchmarks show 2.1-6.4x faster contour geometry and 44-50x faster
EDT than the runnable Python paths. The MAT core is about 1,984-2,389x faster
than the MATLAB source interpreted by Octave on the tested 32/64-pixel cases.
These ratios do not describe MATLAB JIT performance or a complete image pipeline.
Vectorized wrappers make the 1,000-contour R workload 184.9x faster than the
earlier checked R candidate, preserving a mixed-contour result exactly.

See [qualification](docs/QUALIFICATION.md), [benchmarks](docs/BENCHMARKS.md)
and the machine-readable receipts in docs/evidence. The local installable archive
is mlvtoolbox_0.1.0.tar.gz; it is ignored by Git and can be rebuilt with
tools/bounded_check.py.

Prepared CI covers macOS, Windows and Linux. Hosted CI and actual Windows/Linux
runs are unexecuted. No CRAN/R-universe submission, site deployment or upstream
issue submission has occurred.

MATLAB itself and the full Python image pipeline have not been run. Natural
image scientific acceptance, general SVG paint effects, MAT v7.3/HDF5 adapters
and bundled SAM extraction are outside this local qualification. Rasterization,
graph tracing, EDT ties, taper smoothing and mean orientation have explicit
method choices documented in UPSTREAM_BUGS.md.

## Storage and provenance

Cleanup verified 159 source files and 39 retained data files, removed about
2.63 GiB of logical content, and preserved one local reference-data collection.
See CLEANUP_RECEIPT.json, vendor/SOURCES.json and tools/restore_upstream.py.

The R implementation and attributed derived kernels are GPL-3.0-or-later.
Original reference notices and mixed rights remain recorded; this local port
does not assert dataset redistribution rights or blanket upstream relicensing.
