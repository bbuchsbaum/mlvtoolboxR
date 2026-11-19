# CRAN readiness: 0.2.1

Local audit on 2026-10-26: no package errors or warnings remain in the full
CRAN-like check. Submission is not yet qualified on current R-devel or on
Windows/Linux. No submission or publication was performed.

## Changes

- Added the two original Octave morphology lookup-table authors to `Authors@R`
  as contributors, consistent with the preserved notices in `src/pratt_lut.h`
  and `inst/COPYRIGHTS`.
- Excluded `.impeccable`, `.ruff_cache`, `PRODUCT.md`, and `cran-comments.md`
  from source builds. The direct build caught `.ruff_cache`, which the earlier
  curated staging process did not include.
- Added draft submission comments without claiming unexecuted platform checks.

## Evidence

Built directly from the repository with `R CMD build .`, then ran:

```
R CMD check --as-cran --output=.agent-work/cran-readiness-final mlvtoolbox_0.2.1.tar.gz
```

The final build/check used `LC_ALL=en_US.UTF-8`, `R_PROFILE_USER=/dev/null`,
and `R_ENVIRON_USER=/dev/null`. Platform: macOS Sonoma arm64, R 4.5.1.
Unlike the previous qualification, this check did not use `--no-manual`.

- Final: **0 errors, 0 warnings, 3 notes**; exit status zero.
- All **1,332 assertions passed**, with no test warnings or skips.
- Examples, both vignette rebuilds, compiled-code checks and PDF manual passed.
- URL checker passed (one external URL discovered).
- Session-state scan found no unguarded `par`, `options`, or `setwd` changes in
  examples/vignettes. The plotting function restores its graphics parameters.
- Archive inspection confirmed no vendor/reference datasets, local caches,
  planning documents, tracker state, research tools or evidence directories.
- The first direct check had 0 errors, 0 warnings and 3 notes: new submission,
  `.ruff_cache`, and unavailable HTML validation. The cache note is fixed;
  the final run additionally reported inability to verify current time.

The remaining notes are new submission, inability to verify current time, and
skipped HTML validation because `/usr/bin/tidy` is too old. HTML validation
remains unverified; no checks were suppressed to remove these notes.

Final archive: `mlvtoolbox_0.2.1.tar.gz`, 352,835 bytes, SHA256
`d00d44ea0f450c7d6a23fc507e6805e9e4f9ff4a144fb6e607db1699b6354ef5`.
This is based on `8a1452b3af4d8953c1edc6d7f53a49f771f2ab54` plus the
DESCRIPTION and build-ignore changes captured in the evidence patch.

[Evidence](evidence/cran-readiness/) includes initial/final logs, tests,
session information, archive contents/hash and the source patch.

## Before submission

Run the final candidate on current R-devel and Windows/Linux, and validate
HTML with a current HTML Tidy. The existing CI matrix covers release R on
three platforms but was not triggered here. Review the draft
`cran-comments.md` after those results. The scientific source attribution
and existing license notices were inspected; this pass does not establish
new redistribution rights or replace the existing upstream rights inventory.

These boundaries follow the [CRAN submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html)
and [repository policy](https://cran.r-project.org/web/packages/policies.html).
