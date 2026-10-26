# Release A: feature contracts and optional estimators

This local candidate is qualified independently of the stimulus designer. Its
378 package files match the retained [source manifest](evidence/release-a/package-source-manifest.json).
The source is an uncommitted candidate over `85363b4f45efe14303e37d590fa69dfb8c058937`;
that Git revision alone does not identify these changes. A concurrent commit
added only histogram-spike tools, evidence and tracker records, without changing
the package source. The content manifest identifies exactly what was checked.

[Source archive](../mlvtoolbox_0.2.0.tar.gz): 295,163 bytes, SHA256
`7c2d7e20e95d7b61991bcfe09e8372878d90a140b59c49d55049d1fb51954954`.
It is ignored by Git. [Machine-readable receipt](evidence/release-a/qualification.json).

| Gate | Result |
|---|---|
| Build and optimized install | Passed with local R 4.5.1 on macOS arm64 |
| Installed tests | 901 assertions; zero failures, warnings or skips |
| R CMD check --as-cran --no-manual | Zero errors/warnings; two triaged notes |
| R and C++ coverage, tests plus examples | 93.6716%, above the 90% target |
| Existing numeric parity fixtures | Passed unchanged; metadata checks extended explicitly |
| Analytical and invariance gate | Passed single-count turning, cyclic closure, tent integrals, subdivision, reversal, rotation and unit scaling |
| Tracing/simplification sensitivity | All 16 frozen circle cases passed; observed errors 0.1%–14.3% under the 25% limit |
| Workflow integration | Checked physical calibration, affine histogram integration, prediction refusal, serialization, batch/export metadata, missingness and stale MAT invalidation |

The check notes are “New submission” and inability to verify current time.
The installed testthat binary reports being built under R 4.5.2; assertions ran
under R 4.5.1 without test warnings. The vignette rebuilt and examples passed.
Full [check log](evidence/release-a/R-CMD-check.log),
[test output](evidence/release-a/testthat.Rout),
[coverage](evidence/release-a/coverage.json), and
[runtime](evidence/release-a/sessionInfo.txt) are retained.

The existing approved low-disk exception was bounded by a 200 MiB owned build
workspace. Sampled peak usage was 12,342,071 bytes during check and 15,394,056
bytes during coverage. This measures workspace files, not peak RSS or an OS
quota. Build, install, check and coverage took approximately 8.7, 5.0, 126.2 and
14.9 seconds. The coverage monitor initially raced R's temporary-file cleanup;
the failure is retained, the monitor now handles disappearing files, and the
rerun completed with an executed exit-zero receipt.

`mlv2023` remains the default for one release. New branch scores and signed-turn
tent curvature are alternative geometric definitions. Smoothing can cancel
opposite turns; clipping kernels introduces endpoint bias. Raster tracing,
thinning, MAT extraction and mapped-pixel histograms remain resolution-sensitive.
The adverse zigzag and raster cases are retained, without a convergence claim.
See the [measurement contract](FEATURE_SPEC.md) and
[frozen acceptance specification](ESTIMATOR_ACCEPTANCE.md).

MATLAB, Linux, Windows, hosted CI, PDF manual generation, peak RSS and scientific
predictive validation were not executed. This is local candidate qualification,
not cross-platform or perceptual validation. No push, tag, deployment or package
submission occurred.
