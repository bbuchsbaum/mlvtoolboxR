# Contour graph qualification

The frozen matched-marginal gate passed. All four paired drawings have exactly
matching orientation, length and curvature histogram masses, while the graph
descriptors distinguish the predeclared connectivity, continuation, parallel-pair
and enclosure differences. Hand-enumerated cross, T, bowtie, near-miss, empty and
grid topology agrees. Overlap rejection and domain limits are explicit.

The final focused graph suite has 108 assertions. It covers transformations,
calibrated scaling, collinear subdivision, transitive point welding, perturbation
boundaries, source identities, serialization and stale geometry/parameter/
calibration rejection. The original frozen suite completed in 1.559 seconds
under its 10-second local gate. A 10-by-10 line grid took 0.029 seconds and
returned 100 vertices, 180 edges, one component and cycle rank 81.

The checked package includes Release A, the subset designer, final concurrent
fused-statistics changes and the graph implementation. Its 394 files match the
[source manifest](evidence/graph-v1/package-source-manifest.json), and the frozen
graph files still match their pre-implementation hashes. It was checked as an
uncommitted 0.2.1 candidate over `85363b4f45efe14303e37d590fa69dfb8c058937`;
the manifest/archive, not that base commit alone, identifies these bytes.

[Source archive](../mlvtoolbox_0.2.1-graphs.tar.gz): 352,771 bytes, SHA256
`de072ae04587283c3d2f5f49ac6d0e9101341282cac4aa1ac92524a8584d0283`.
[Machine-readable receipt](evidence/graph-v1/qualification.json).

| Gate | Executed result |
|---|---|
| Installed package tests | 1,332 assertions, zero failures/warnings/skips |
| R CMD check --as-cran --no-manual | Zero errors/warnings; New submission and current-time verification notes |
| R and C++ coverage | 95.0669% |
| Examples and vignette rebuilds | Passed |
| Published parity and earlier feature/designer gates | Passed |
| Source and graph-freeze verification | Zero mismatches |

Build/install/check took about 9.9/6.3/131.0 seconds, coverage 18.5 seconds.
Sampled workspace peaks were 14,921,402 bytes for check and 21,259,936 for
coverage, below the existing 200 MiB owned-workspace limit. Full logs and adverse
earlier check receipts are retained in `docs/evidence/graph-v1`. An initial
successful check preceded the stale-calibration correction; the final result
above comes from the rerun after that fix. Tests ran on macOS arm64, R 4.5.1.

Linux, Windows, MATLAB, hosted CI and PDF manual generation remain unexecuted.
These are geometric descriptors with explicit domain rules, not demonstrated
border ownership, global symmetry or natural-image predictive validity. No
push, tag, deployment or package submission occurred. The later Voronoi research
prototype does not modify package methods or dependencies.
