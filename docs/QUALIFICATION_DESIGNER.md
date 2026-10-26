# Release B: bounded whole-contour subset design

The precommitted synthetic comparative gate passed without changing the frozen
fixtures, objective or tolerances. Both required knapsack cases improved raw
integrated-curvature separation from ranked 100 to optimal 270 geometric degrees
under the same length, count and orientation constraints. Feasible controls tied.
These deliberately constructed cases establish combinatorial behavior, not
natural-image usability or a perceptual advantage.

| Frozen case | Exact designer | Ranked comparator | Outcome |
|---|---:|---:|---|
| One-drawing knapsack | 270 | 100 | Required improvement passed |
| Set-level knapsack, unequal canvases | 270 | 100 | Required improvement passed |
| Tied straight contours | 0 | 0 | Tie retained |
| Two one-contour drawings, pooled | 90 | 90 | Tie retained |
| Same pair, per-drawing matching | Proven infeasible | Infeasible | No forced split |
| Smaller retained-length budget | 178 | Infeasible | Reported; not counted as a win |
| One-node budget | Exhausted, no solution | 100 | No infeasibility/optimality claim |

Every claimed feasible output was independently remeasured and checked against
all constraints. Completed searches match a separate exhaustive ternary oracle.
No heuristic fallback or silent relaxation was used. Complete residuals, source
IDs, bound/gap status, wall times and memory receipts are retained for every case
in [benchmark evidence](evidence/designer-v1/results.csv); the serialized results
and readable `*-diagnostics.txt` files sit beside it.

Search took 0–2 ms and 4–259 nodes on these tiny cases; end-to-end API calls took
17–55 ms. The whole measured processes, including R startup, the comparator and
oracle, peaked at roughly 90–99 MiB RSS. A supplemental 2 GiB process cap was
declared before comparative execution; this did not alter the frozen cases or
tolerances. The first memory-monitor attempt was blocked by the sandbox; its
failure is retained, and the authorized rerun measured every case successfully.
These timings are small synthetic observations, not large-experiment guarantees.

## Exact candidate and package checks

The checked 388-file candidate includes the concurrent fused-statistics work
and uses package version 0.2.1. It is an uncommitted snapshot over
`85363b4f45efe14303e37d590fa69dfb8c058937`, identified by the
[source manifest](evidence/designer-v1/package-source-manifest.json), not by that
base commit alone. Source and frozen-fixture verification found zero mismatches.

[Retained archive](../mlvtoolbox_0.2.1-designer.tar.gz): 344,296 bytes, SHA256
`18aebd0f8f7434d51f9417729408f54893db567aa526c262a4dd6da386cb6451`.
It is ignored by Git. [Qualification receipt](evidence/designer-v1/qualification.json).

* Installed tests: 1,214 assertions, zero failures, warnings or skips.
* R CMD check: zero errors/warnings. Notes: New submission; unable to verify
  current time. Examples and both rebuilt vignettes passed.
* R and C++ coverage: 94.6411%, including new designer and fused-statistics code.
* Published numeric parity, feature guards, estimator gates, physical-unit and
  stale-cache tests continue to pass.
* The separate owned `designer-build` workspace stayed below the existing
  200 MiB cap: sampled check peak 14,782,391 bytes; coverage peak 21,600,997 bytes.

Build/install/check took about 9.3/6.1/136.3 seconds; coverage took 17.2 seconds.
The runtime was local macOS arm64, R 4.5.1. The testthat binary's R 4.5.2 build
notice is distinct from test warnings. Full logs, source/fixture hashes and
runtime metadata are retained under `docs/evidence/designer-v1`.

The public workflow preserves original contour identity, specification,
per-input calibration, provenance and geometry through R serialization.
It reconstructs selected drawings before verifying residuals. Optional junction
and MAT diagnostics are recomputed afterward and are never optimized constraints.
See the [workflow guide](SUBSET_DESIGNER.md) and the runnable stimulus-design
vignette for feasible, infeasible and budget-limited examples.

Linux, Windows, MATLAB, hosted CI, PDF manual generation, natural-image usability
and perceptual validity remain unexecuted. Graph descriptors and the Voronoi
prototype are outside this candidate. No push, tag, deployment or submission
occurred. Later source changes require their own qualification.
