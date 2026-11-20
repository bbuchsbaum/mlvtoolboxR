# Local qualification: 0.1.1

Checked package source: `d83ed952fdc9f639b98e849b0d8697259894de46`.
Later delivery changes add development tools, documentation and receipts;
the 359 package-source files still match the checked manifest.

Local source archive: [mlvtoolbox_0.1.1.tar.gz](../mlvtoolbox_0.1.1.tar.gz),
270,147 bytes, SHA256
`1860a64ac5ac92284972f2e9ba76b8f78e97064e3634dd1fedc79d42e11b668d`.
It is ignored by Git and can be regenerated with `tools/bounded_check.py`.
The previous 0.1.0 archive and its historical receipts are preserved.

| Gate | Executed result | Evidence |
|---|---|---|
| Build and optimized install | Exit 0, C++17 clang 15 / -O2 | [build](evidence/parity-perf/build-09.log), [install](evidence/parity-perf/install-08.log) |
| R CMD check --as-cran --no-manual | 0 errors, 0 warnings, 1 NOTE: New submission | [full check](evidence/parity-perf/R-CMD-check-final.log) |
| Installed test suite | 583 assertions, 0 failures/warnings/skips; previously 250 | [test output](evidence/parity-perf/testthat-final.Rout) |
| R and C++ coverage, tests + examples | 92.8989%; previously 92.6117% | [summary](evidence/parity-perf/coverage.json), [lines](evidence/parity-perf/coverage-lines.csv) |
| Synthetic reference comparisons | 18 cases, 289 comparisons; 226 matches, 63 classified differences, none unexplained | [full parity report](PARITY.md) |
| Retained forest drawing | Four geometry/histogram comparisons pass after one explicit zero-length removal | [receipt](evidence/parity-perf/natural-validation-comparison.csv) |
| Performance output preservation | All five old/new workload outputs identical | [exact receipt](evidence/parity-perf/exact-workflow-outputs.json) |
| Source and fixture integrity | 359 package files and 264 fixture CSVs, zero mismatches | [verification](evidence/parity-perf/source-verification.json), [source manifest](evidence/parity-perf/package-source-manifest.json) |
| Vendor/data custody | 159 source files, 39 data files, zero mismatches | [verification](evidence/parity-perf/reference-verification.json) |
| Scientific plot inspection | Four differing MAT cases inspected; controlled label substitution resolves all six skeleton comparisons | [figure](evidence/parity-perf/mat-parity.png), [assertions](evidence/parity-perf/visual-parity.log) |

The new tests combine executed Octave/Python fixtures with analytical axial
boundaries, length conservation, exact old-output witnesses, exhaustive spatial
mapping, exact nearest-label distance checks and independent line-fit projections.
The full check also passed examples and vignette rebuilds. No test skip or warning
suppression was introduced. The installed testthat binary reports being built
under R 4.5.2 while this runtime is 4.5.1; test assertions themselves emit no warnings.

Two R-port revisions were made: isolated-pixel signed distances and odd-bin
axial histogram boundaries. Six new upstream findings are recorded in
[DIFFERENCES](../DIFFERENCES.md), with raw failures and numerical differences
retained. No upstream issue was submitted.

## Resources and evidence

The already approved low-disk exception remains bounded by the owned-workspace
monitor at 200 MiB. Check peak usage was 11,822,896 logical bytes; instrumented
coverage peak was 16,994,824 bytes. This is a sampled monitor, not an OS disk
quota. [Check metadata](evidence/parity-perf/check-08.json),
[coverage metadata](evidence/parity-perf/coverage-04.json).

After saving the archives and receipts, session-owned temporary installations,
check trees, raw natural-data working copies and compiled objects are removed.
The verified vendor snapshots and canonical reference data remain untouched.
Regeneration scripts and committed synthetic fixtures retain the reproducible inputs.
[Cleanup receipt](evidence/parity-perf/temporary-cleanup.json): 14,612,232 logical
bytes removed from this session's owned temporary paths.

All reports and receipts for this follow-up are under `docs/evidence/parity-perf`;
historical 0.1.0 evidence is not overwritten. Initial harness failures and
pre-fix mismatches are retained beside successful runs. Logs identify commands,
exit codes and durations. [Machine-readable qualification](evidence/parity-perf/qualification.json);
[runtime metadata](evidence/parity-perf/sessionInfo.txt).

## Boundaries

Default R/Octave MAT skeletons differ on four of six new masks at equidistant
nearest-boundary ties. Shared-label flux and skeletons agree; universal end-to-end
MAT parity is not claimed. The forest witness establishes geometry/histograms,
not scientific natural-image acceptance. Python shape-adapted paths are labeled.

Actual MATLAB, Windows/Linux, hosted CI, full Python image processing, peak RSS
and parallel scaling remain unexecuted. PDF manual generation was omitted with
`--no-manual`. No push, package submission, site deployment or upstream filing occurred.

See [performance](PERFORMANCE_0.1.1.md) for paired measurements and
[parity](PARITY.md) for detailed compatibility boundaries.
