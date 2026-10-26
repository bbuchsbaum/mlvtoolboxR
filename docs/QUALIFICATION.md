# Local qualification receipt

Current follow-up: [0.1.1 qualification](QUALIFICATION_0.1.1.md). This page retains the 0.1.0 receipt.

The local R package `mlvtoolbox` 0.1.0 was built, installed and checked at
source commit `141a4a44dcc6431faa53a6a4b8af31b06b306982`.
The final source archive is [mlvtoolbox_0.1.0.tar.gz](../mlvtoolbox_0.1.0.tar.gz)
(107,812 bytes), SHA256
`ed7d49c60174011d80b299d1dda7f73be5cb8f049e26680034c9002f08cc30c7`.
The archive is a local ignored build artifact; regenerate it with the bounded
builder after cloning. Later delivery commits contain documentation and evidence.

## Executed gates

| Gate | Result | Evidence |
|---|---|---|
| Source build and optimized install | Exit 0; C++17, clang 15, -O2 | [build](evidence/build-07.log), [install](evidence/install-06.log) |
| `R CMD check --as-cran --no-manual` | 0 errors, 0 warnings, 1 NOTE | [full log](evidence/R-CMD-check-final.log) |
| Regression tests in installed package | 250 passing assertions; 0 failures/warnings/skips | [test output](evidence/testthat-final.Rout) |
| Examples and vignette build/rebuild | Passed within R CMD check | [check log](evidence/R-CMD-check-final.log) |
| R and C++ line coverage, tests + examples | 92.6117%, from 82.8125% initially | [summary](evidence/coverage.json), [lines](evidence/coverage-lines.csv), [initial run](evidence/coverage-initial.log) |
| Numerical/reference tests | Passed analytical, exhaustive, metamorphic and executed Octave witnesses | tests/testthat and inst/fixtures |
| Final package source manifest | All 91 staged files match the live source bytes | [manifest](evidence/package-source-manifest.json) |
| Vendor/data integrity | 159 source files, 39 data files; zero mismatches | [verification](evidence/reference-verification.json) |
| Matched stage performance | Faster on all measured common correct stages | [benchmark report](BENCHMARKS.md) |
| Plot inspection | Synthetic length, curvature, junction and local MAT panels inspected | [preview](evidence/workflow-preview.png) |

The remaining NOTE is CRAN's informational “New submission”; no submission was
made. The first retained successful check also had 0 errors, 0 warnings and
this one NOTE. A local testthat binary reports that it was built under R 4.5.2
while this runtime is R 4.5.1; the assertions themselves emit no warnings.
No warning suppression or test skips were added to achieve the result.

Tests cover public drawing/image workflows and error boundaries, conservative
selection budgets, circle crossings with both endpoints outside, negative
coordinates, transitive junction merging, axial angle conventions, input
immutability, valid-zero score mapping and malformed adapters. The final added
workflow runs a PNG RGB provider through contour tracing, MAT scoring, feature
mapping, selection and plotting. The wrapper optimization has an exact saved
pre-change result witness.

The 25×35 Octave MAT fixture has identical skeletons (44 pixels), nearest labels
agree on that fixture, and maximum AOF discrepancy is approximately 1.1e-14.
Three Pratt morphology fixtures also agree exactly. Further benchmark fixtures
at 32×32 and 64×64 pass explicit numerical and pixel comparisons.
These are bounded synthetic witnesses, not universal equivalence claims.

## Resources and retained adverse evidence

The user explicitly approved running below the parent directory's 20 GB free
space threshold. The builder monitors only its owned workspace, with a 200 MiB
limit and process-group-specific cancellation. The final check peak was
9,060,743 logical bytes; the coverage peak was 11,787,058 bytes.
See [check metadata](evidence/check-06.json) and
[coverage metadata](evidence/coverage-03.json). This is a sampled monitor,
not an operating-system disk quota.

The initial direct source-build attempt copied retained data before applying
R's exclusions and was stopped by the monitor. Its
[failure receipt](evidence/build-cap-failure.json) and
[log](evidence/build-cap-failure.log) are retained.
The fixed builder prestages the tracked package allowlist, then invokes R.
A fresh staging tree prevents stale files from surviving between runs.
The final tarball excludes vendor snapshots, reference data, tracker state and
development evidence, while including the prebuilt vignette index.

The first Octave benchmark attempt omitted loading the image package and failed.
Its [log](evidence/benchmark-octave-first-attempt.log) is retained; the corrected
run completed all four cases. Unsupported or incorrect Python paths were
documented rather than assigned benchmark scores.

The cleanup removed 2,824,190,787 logical bytes (about 2.63 GiB) of redundant
reference/history/generated content. One local canonical collection remains at
`reference-data/mlv`, outside the R package and Git. Its 39 data files total
698,453,722 bytes; the cleanup's directory total additionally includes its
manifest. Source SHAs, hashes, volume information and restoration instructions
remain in [SOURCES](../vendor/SOURCES.json),
[cleanup receipt](../CLEANUP_RECEIPT.json) and
[restore helper](../tools/restore_upstream.py).
Run `python3 tools/verify_references.py` to recheck custody.

After retaining these receipts and the final archive, this session's owned
temporary check trees, coverage objects, debug objects and temporary install
library are removed. No reference sources or canonical data are removed.

## Qualification boundary

- Actual MATLAB execution, Windows/Linux execution and hosted CI are unexecuted.
  The three-OS CI workflow is prepared, without a claimed hosted verdict.
- Natural-image scientific acceptance, peak RSS, parallel scaling and complete
  image-pipeline performance remain unmeasured.
- MATLAB rasterization/morphology and edge-linker equivalence are not asserted.
  Pratt morphology, pixel graph tracing, tie rules, axial mean orientation and
  taper smoothing are explicit methods.
- General SVG paint effects and MAT v7.3/HDF5 are outside the implemented adapters.
  RGB extraction requires an explicit provider; SAM weights and legacy serialized
  MATLAB models are not bundled.
- Manual/PDF generation was omitted with `--no-manual`.
- No upstream issue, repository push, package submission or site publication
  occurred. Findings distinguish executed bugs from source findings and
  scientific contract questions in [UPSTREAM_BUGS](../UPSTREAM_BUGS.md).

Machine-readable details: [qualification.json](evidence/qualification.json).
Runtime metadata: [sessionInfo](evidence/sessionInfo.txt).
