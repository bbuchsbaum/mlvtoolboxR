# Medial-axis (MAT) qualification against stored MATLAB outputs (track B)

Plan: docs/plans/mlv-reference-qualification.md, phase P3. Date 2026-11-19.
Principle: agreement with MLV is evidence, not the goal; each divergence names its oracle and who was wrong. IDs follow track A's (M018–M026, R005–R009): MLV-M027–M031 and MLV-R010–R013.

## Summary
- MAT stage: with aof="mlv2023" (morphology default "hmt8"), R reproduces MATLAB's stored Forest_DBW distance map (bit-exact, single), AOF (2.1e-14) and skeleton (0 px).
- New branch_rule axis. Default branch_rule="invariant": tracing, scoring and mapping are exactly invariant to branch reversal, branch renumbering, contour renumbering and the 8 lattice symmetries (tests, 1e-12). branch_rule="mlv2023": the published chain (Kovesi edgelink + cleanedgelist ported line by line, forward MLV scores, MLV fitLineSegments, last-writer rating image, MATLAB tangent enumeration and tie order, zeros as missing).
- Published-compatible mode (aof="mlv2023", morphology="hmt8", branch_rule="mlv2023") on ALL 2,757 drawings with stored MAT samples: per-pixel within 1e-4 = 0.991 parallelism, 0.994 separation, 0.993 mirror; per-contour mean r = 0.9993 / 0.9999 / 0.9986; per-drawing mean r = 1.000. Residual = MLV single precision (section 11).
- Defaults (centered, hmt8, invariant) deliberately differ from published values (per-contour r 0.891 / 0.978 / 0.842).

## Inputs

| Input | SHA-256 | Role |
|---|---|---|
| reference-data/mlv/examples/Forest_DBW.mat | da1f7ed1e512913b0210f8aefcb83f8ef673b82fdcb7e1c5d62da2bd778fb5 | MATLAB skeleton, distance_map (single), AOF (768x1024); per-pixel *_allX/_allY/_allScores |
| reference-data/mlv/examples/MAT.mat | 96bd08c78f68d0c1ac53e99c23cd98c40d2e936ec9179a4b0a1056226e8c47a2 | pairs with imgLD.mat; other code revision |
| reference-data/mlv/examples/imgLD.mat | 3209102e42c0d6372900f97cc7b0f821f501089851b6cad4daa5e50a9a7d5e74 | 600x800 binary drawing |
| .agent-work/mlvref/<dataset>/*.csv.gz | P0 export manifest | 29 datasets with MAT samples (2,757 drawings) |

Natural data never leave reference-data/ and .agent-work/; this directory holds aggregate CSVs and synthetic PNGs only.

Export note: vecLD.parallelism/separation/mirror{c} are per-PIXEL vectors (MATpropertiesToContours rasterizes each contour), never per-segment, so the score columns of segments.csv.gz are NaN by construction. All comparisons use the _all samples; per-contour MATLAB means are recomputed from them over each contour's raster (exact, section 1).

## Harness (repository root, MLVREF_LIB=<library>)
- tools/mlvref/mat_common.R: stage helpers, controlled-substitution MAT objects, MLV mapping emulation (map_at, tangent_matlab).
- mat_ties.cpp: exact equidistant tie sets; AOF in MATLAB arithmetic order; tie-mean and centred variants.
- mat_thin.cpp: Guo-Hall, Zhang-Suen, ordered hit-or-miss thinning, thinness oracle.
- mat_octave_chain.m: pinned traceSkeleton -> computeMATproperty -> mapMATtoContour under Octave on supplied MAT fields. mat_octave_aof.m: pinned computeGradientVectorField + computeAOF.
- mat_forest_ties.R, mat_equivariance.R, mat_raster.R, mat_forest_stages.R, mat_branch_invariance.R, mat_forest_drift.R, mat_example_pair.R, mat_thin_validate.R, mat_trace_merge.R (prototype tracer), mat_datasets.R + mat_datasets_summary.R (end to end), mat_figures.R, mat_refresh_witnesses.R.

Octave 11.1 + image 2.20.0 adapters (temporary copies, never in vendor/): jj++m_Neighbors8 -> jj + m_Neighbors8; fixed-window smoothdata(x,'movmean',k) (the only forms reached by parallelism/separation/mirror); vectorised exact knnsearch (exact ties -> lowest row); edgelink.m and findendsjunctions.m with bwmorph(...,'skel',Inf) as identity, because the input is MATLAB's own skeleton and Octave's Pratt skel adds 12 bridge pixels to it.

## Tolerances
MATLAB stores distance_map as single, and MLV computes radii, branch scores and tangent points in single; the pinned chain reproduces stored scores to <=1.3e-6. Scores: |error| <= 1e-4. Distances: bit-exact after rounding R to single. AOF: 1e-12 (R and MATLAB differ only in normalise-then-dot vs dot-then-divide; observed <= 2.2e-14). No tolerance was loosened to make a comparison pass.

## 1. Raster recovery and stage isolation (Forest_DBW)
- Raster recovery: ink <=> D <= 0 holds. Recomputing the MLV signed distance (computeGradientVectorField construction) from that mask reproduces all 786,432 stored values bit-exactly in single. The 54 pixels at D = -1 are ink pixels on the image border (never "outer"), not isolated pixels. The recovered set is MATLAB's thinned ink, and it is a fixed point of both thinnings.
- Nearest-boundary ties: R's lowest-linear-index labels with MLV sampling reproduce the stored AOF on every pixel (max 1.6e-14); 178,273 pixels sample a tied cell (tie sets: 36,043 of size 2, 1,000 of 3, 76 of 4, 4 of 5). Highest-index labels would change 163,437 pixels. MATLAB bwdist breaks ties by lowest index.
- Rasterization (proxy for MATLAB's raster = thin ink plus every stored sample position, 23,751 px; all endpoints are integers): the 0.2.1 Bresenham raster had 23,754 px, 427 proxy pixels missing, 430 outside, and reversing every segment changed 1,984 px. All 1,040 half-pixel tie cases in MATLAB follow one rule: positive slope takes the larger minor coordinate, negative slope the smaller. The new package rasterizer (major-axis DDA with that rule): 23,756 px, 0 missing, 5 outside (pixels with zero score), reversal changes 0 px (forest_raster_rules.csv).
- Thinning: Octave/R Pratt skel is not thin (on the flux mask, 18,035 redundant simple non-end pixels vs MATLAB 327; 2,376 2x2 blocks vs 2; 5,982 MATLAB pixels missing, 15,672 extra). Sequential 8-direction hit-or-miss thinning (edge + corner elements) was searched over all 40,320 direction orders: exactly 3 reproduce MATLAB's 24,904-pixel skeleton with 0 differences, all with prefix S, SW, E, SE, W. The adopted order (S,SW,E,SE,W,NW,NE,N) also reproduces the independent ink-thinning stage (R raster -> MATLAB thin ink) exactly (mat_thin_validate.R). Guo-Hall (= Octave 'thin') differs by 11,854 px; Zhang-Suen by 8,507.
- Net: with the new rasterizer and morphology="hmt8", aof="mlv2023", compute_mat on the Forest vecLD reproduces MATLAB's distance map (0 single mismatches), AOF (max 2.1e-14) and skeleton (0 px) exactly.
- MAT.mat + imgLD.mat: they pair (imgLD thin ink vs MAT.mat D<=0 differs by 19 px under hmt8) but MAT.mat was produced by another code revision: distance_map is double, 4 outer-boundary pixels differ, and its "AOF" field is a binary mask containing the skeleton. Not used for pinned-code claims.

## 2. Version drift (P1)
- Tracing/scoring/mapping: stored = pinned. The pinned chain under Octave (traceSkeleton -> computeMATproperty -> mapMATtoContour) on MATLAB's stored fields with single radii reproduces every stored per-pixel score to <=1.3e-6 (mirror bit-exact). With double radii, parallelism errors reach 0.60 via flipped exact nearest-tangent ties.
- D and AOF: consistent with the pinned algorithm (R's exact implementation reproduces them bit-for-bit). Octave is not a valid stand-in: its bwdist is inexact (max |dD| = 0.042) and its labels differ, so Octave's AOF differs at 67,161 px. Octave portability finding, not MLV.
- Morphology: MATLAB bwmorph('skel') is undocumented; identified empirically above (exact on two Forest stages).

## 3. Ties and the AOF definition
Effect of the tie choice on Forest, hmt8 morphology (forest_tie_policy.csv, forest_tie_downstream.csv):

| Policy | flux-mask px changed | skeleton px changed (of 24,904) | mapped parallelism changed |
|---|---|---|---|
| lowest (MATLAB) | 0 | 0 | 0 |
| highest | 10,224 | 8,632 | 68% of px |
| tie-mean, MLV sampling | 5,204 | 2,626 | 34% |
| centred tie-mean (new default) | 6,122 | 3,285 | 34% |

MLV-M027: computeAOF centres the sampling circle on the pixel corner (i+.5,j+.5) but reads labels from floor() cells, so the flux is misregistered by half a pixel and not lattice-equivariant. Oracle 1 (analytical): ink rows 10 and 22, medial row 16: mlv2023 flux exceeds threshold on rows 16 AND 17 (34.4, 32.4); centred tie-mean exceeds it only on row 16, with symmetric neighbours (17.58 at rows 15 and 17) (aof_registration_synthetic.png). Oracle 2 (D4 equivariance, 6 masks x 4 transforms on the correct lattice, equivariance.csv): every MLV-sampled policy fails (lowest, highest, tie-mean, symmetric tie-mean; up to 52 flux units); centred + lowest also fails (466 px); only centred + nearest-pixel cell + mean unit vector over all equidistant boundary pixels passes (max 5.0e-14, 0 mask changes). The package tests repeat this for all 8 symmetries on 5 masks.

Decision implemented (user-approved): default aof="centered"; aof="mlv2023" kept as the explicit published-compatible option.

## 4. End-to-end agreement with stored MATLAB
Setup: stratified 20 drawings per dataset (seed 20261119), 580 drawings, ~6.8 M stored pixels, 146,645 contours; threshold 0.3 (the omitted-threshold value used by computeAllMATfromVecLD, MLV-M004); geometry from segments.csv, zero-length segments dropped. Per-contour MATLAB mean = mean nonzero mapped score over the contour's raster, reproduced exactly by the package rasterizer. Receipts: e2e_summary.csv (per family x config x property), e2e_contour_summary.csv. Runtime ~2-9 s per drawing per config; a full 2,757-drawing run is feasible (~1 h on 6 cores per config).

Per pixel (fraction within 1e-4 / within 1e-2 / mean abs / per-drawing-mean r):

| config | parallelism | separation | mirror |
|---|---|---|---|
| hmt8 + mlv2023 (published-compatible) | 0.343 / 0.505 / 0.056 / 0.996 | 0.428 / 0.702 / 0.015 / 1.000 | 0.299 / 0.443 / 0.048 / 0.985 |
| hmt8 + centered (new default) | 0.215 / 0.379 / 0.071 / 0.995 | 0.360 / 0.663 / 0.019 / 1.000 | 0.171 / 0.372 / 0.054 / 0.986 |
| pratt + mlv2023 (old default) | 0.055 / 0.112 / 0.221 / 0.890 | 0.093 / 0.367 / 0.068 / 0.991 | 0.025 / 0.073 / 0.141 / 0.841 |

Coverage of stored pixels: 0.994 (hmt8), 0.995 (centered), 0.930 (pratt).

Per contour (r / median abs / within 0.05 / mean bias R-MATLAB):

| config | parallelism | separation | mirror |
|---|---|---|---|
| mlv2023 | 0.898 / 0.021 / 0.703 / -0.025 | 0.985 / 0.006 / 0.922 / -0.001 | 0.829 / 0.027 / 0.674 / -0.008 |
| centered | 0.871 / 0.029 / 0.648 / -0.019 | 0.979 / 0.008 / 0.900 / -0.002 | 0.809 / 0.031 / 0.644 / -0.004 |
| pratt | 0.600 / 0.124 / 0.242 / -0.123 | 0.886 / 0.033 / 0.620 / -0.025 | 0.538 / 0.104 / 0.263 / -0.099 |

Per-family ranges of per-contour r (mlv2023 / centered): parallelism 0.856-0.938 / 0.824-0.925; separation 0.975-0.995 / 0.965-0.993; mirror 0.773-0.887 / 0.751-0.875. Contours defined in MATLAB but not R: 305 (mlv2023), 229 (centered), 3,005 (pratt). pratt failed on 7 of 580 drawings (MLV-R012); hmt8 configs had no failures.

Stored-data notes: 11 separation samples (one IAPS drawing) are +Inf (radius 0: (1-1/0)^10), excluded from mean-abs; 41,284 parallelism samples exceed 1 by <=1.2e-6 (MLV single-precision arithmetic).

### Per dataset family (pixels within 1e-4 / per-contour r)

aof="mlv2023" (published-compatible), hmt8:

| family | n | parallelism | separation | mirror |
|---|---|---|---|---|
| beaches | 20 | 0.349/0.915 | 0.438/0.995 | 0.291/0.847 |
| cities | 20 | 0.344/0.883 | 0.451/0.976 | 0.280/0.841 |
| forests | 20 | 0.357/0.926 | 0.461/0.989 | 0.322/0.883 |
| high_ceiling | 20 | 0.395/0.906 | 0.487/0.983 | 0.318/0.820 |
| highways | 20 | 0.351/0.856 | 0.434/0.983 | 0.291/0.839 |
| IAPS | 240 | 0.341/0.914 | 0.416/0.988 | 0.305/0.830 |
| low_ceiling | 20 | 0.411/0.920 | 0.504/0.985 | 0.301/0.831 |
| mountains | 20 | 0.338/0.938 | 0.400/0.994 | 0.299/0.887 |
| OASIS | 180 | 0.310/0.865 | 0.398/0.975 | 0.283/0.773 |
| offices | 20 | 0.336/0.894 | 0.432/0.984 | 0.300/0.848 |
| ALL | 580 | 0.342/0.898 | 0.428/0.985 | 0.299/0.829 |

aof="centered" (default), hmt8:

| family | n | parallelism | separation | mirror |
|---|---|---|---|---|
| beaches | 20 | 0.244/0.882 | 0.382/0.993 | 0.167/0.833 |
| cities | 20 | 0.232/0.845 | 0.380/0.973 | 0.147/0.808 |
| forests | 20 | 0.227/0.890 | 0.393/0.985 | 0.179/0.858 |
| high_ceiling | 20 | 0.272/0.894 | 0.415/0.976 | 0.182/0.795 |
| highways | 20 | 0.230/0.824 | 0.360/0.979 | 0.157/0.818 |
| IAPS | 240 | 0.210/0.880 | 0.351/0.983 | 0.178/0.813 |
| low_ceiling | 20 | 0.296/0.902 | 0.436/0.976 | 0.175/0.815 |
| mountains | 20 | 0.210/0.925 | 0.343/0.993 | 0.164/0.875 |
| OASIS | 180 | 0.175/0.845 | 0.329/0.965 | 0.167/0.751 |
| offices | 20 | 0.208/0.844 | 0.353/0.980 | 0.163/0.831 |
| ALL | 580 | 0.215/0.871 | 0.360/0.979 | 0.171/0.809 |

Forest stage attribution by controlled substitution (within 1e-4 against stored MATLAB samples; forest_stage_attribution.csv):

| Stage | Configuration | Parallelism | Separation | Mirror |
|---|---|---|---|---|
| A | MATLAB skeleton + edgelink branches + MLV mapping order (pinned-equivalent) | 0.999 | 0.999 | 0.927 |
| B | MATLAB skeleton + edgelink branches + R mapping | 0.848 | 0.935 | 0.907 |
| C | MATLAB skeleton + R tracing + MLV mapping order | 0.352 | 0.460 | 0.292 |
| D | MATLAB skeleton + R tracing + R mapping | 0.304 | 0.458 | 0.293 |
| E | R hmt8 + mlv2023 (MAT identical to MATLAB) + R tracing/mapping | 0.304 | 0.458 | 0.293 |
| F | R pratt + mlv2023 | 0.034 | 0.123 | 0.021 |
| G | R defaults (hmt8 + centered) | 0.207 | 0.398 | 0.182 |

A's residual is single precision; mirror 0.927 reflects R's intentional MLV-M016 fix. E = D: in published-compatible mode the MAT stage contributes nothing. The default flux (E -> G) changes 3,285 skeleton pixels on Forest, by design.

So in published-compatible mode the whole remaining gap is three conventions downstream of the MAT: branch tracing, branch traversal direction, mapping order. Branch direction alone: reversing every edgelink branch (same geometry) drops parallelism agreement from 0.999 to 0.044 (separation 0.099, mirror 0.111); removing edgelink's duplicated merge points drops it to 0.833 (forest_branch_invariance.csv). A spur-pruned, junction-merged tracer prototype (mat_trace_merge.R) matches 462 of 519 edgelink branch geometries but still scores ~0.35, for the direction/duplicate reasons.

## 5. Adjudications

**MLV-R010 (R wrong, FIXED): width-1 rasterization depended on segment direction.** Oracle: reversal invariance; (1,1)-(5,2) gave pixel 22 forward and 21 reversed; Forest: 1,984 px changed under reversal. Fix (src/kernels.cpp, cpp_raster_indices): major-axis DDA between rounded endpoints, exact half-pixel ties to the larger minor coordinate for positive slope and smaller for negative; reversal-invariant and 0-missing against MATLAB's raster proxy. Tests: analytic pixel lists both slopes and steep cases, 100 random segments reversed, pixel count = max(|dx|,|dy|)+1.

**MLV-R011 (R wrong, FIXED): .tangent_points produced NaN subscripts and errored on coincident consecutive branch points** (edgelink branches carry 272 on Forest). Fix (R/mat.R): a zero-length pair contributes no tangent samples (MLV instead uses the centre point). Test included.

**MLV-R012 (R wrong, OPEN): Pratt skeletons can include ink pixels at D = -1;** score_branch rejects negative radii, aborting compute_mat_properties on 7 of 580 natural drawings (pratt path). Not observed with hmt8 in 580 drawings. Needs a contract decision.

**MLV-R013 (R fidelity defect, FIXED by default switch): Pratt thinning is neither thin nor MATLAB's skel.** Oracles: thinness (redundant simple non-end pixels, 2x2 blocks), homotopy, and exact agreement with MATLAB on two Forest stages. hmt8 is now the default; Pratt kept as morphology="pratt"/method="pratt" and still covered by the Octave parity tests. Note: both MATLAB skel and hmt8 are sequential and therefore not lattice-equivariant on thick input; equivariance is asserted for the flux.

**MLV-M027 (MLV wrong): computeAOF half-pixel misregistration;** not D4-equivariant; ridge on rows {16,17} for exact medial row 16. Oracles: analytical midline + D4 equivariance (section 3). Default is now aof="centered"; mlv2023 kept for published values.

**MLV-M028 (MLV design defect): mlv2023 branch scores depend on traversal direction** (forward differences padded at the end, movmean(2) radius window, tangent pair (i,i+1) scored by point i). Oracle: reversal invariance; on a synthetic branch reversal changes parallelism by 0.029, separation 0.009, mirror 0.21; on Forest it drops stored-value agreement to 0.044. R's arclength estimator is reversal-invariant (<=4e-15). Test included.

**MLV-M029 (MLV defect): mapMATtoContour reads scores from a last-writer-wins skeleton image,** so a branch's tangent point can receive another branch's score at shared pixels; and exact nearest-tangent ties are resolved by enumeration order ([FP1;FP2], i -> r1 -> r2). Emulating both moves Forest parallelism agreement from 0.848 to 0.999. R keeps each branch's own score; ties are underdetermined (recommend max over exact ties).

**MLV-M030 (MLV defect, source-confirmed): Kovesi cleanedgelist mergenodes concatenates edges including the shared node,** duplicating junction points (272 on Forest); duplicates inflate mirror's point counts and shift windows (removing them: parallelism agreement 0.833, mirror 0.710).

**MLV-M031 (MLV defect): separation = (1-1/R)^10 is unbounded:** radius 0 gives +Inf (11 stored samples) and radii below 0.5 give scores above 1. R returns NA below radius 1.

### Underdetermined, documented
- R tracing (split at every junction, drop edges <10 points) vs MLV (prune spurs, merge through degree-2 nodes; topology more principled, but duplicates and order dependence are not).
- Flux equal to threshold (MATLAB keeps it, R uses >; 0 cases on Forest).
- Single precision in MLV (scores up to 1.2e-6 above 1).
- Rasterizer ties for non-integer endpoints: endpoints are rounded first; MATLAB's tie rule was observed only on integer endpoints.

### Octave portability (not MLV)
Image 2.20.0 bwdist is inexact (0.042) with different labels; Octave Pratt skel != MATLAB skel.

## 6. Default switch, versioning and witnesses

*Superseded in part by section 15: the default morphology is now `"minimal"`; `"hmt8"` is the published-compatible thinning.*
- compute_mat(aof=c("centered","mlv2023"), morphology=c("hmt8","pratt")); skeletonize(method=c("hmt8","pratt")); trace_contours thins with hmt8 or the spec's morphology and records the morphology definition in its provenance.
- feature_spec: mat$aof (default "centered") and mat$morphology (default "hmt8"), validated; definitions aof = "centered-tie-mean-1" | "mlv2023-computeAOF-lowest-index", morphology = "sequential_hmt8_matlab_skel_order-1" | "octave_image_pratt_2.20.0". MAT parameters record aof, morphology, aof_definition, morphology_definition; .mat_semantics includes the definitions.
- Cache validation (existing mechanism): 0.2.x specs (Pratt definition, no aof) fail validate_feature_spec as noncanonical; compute_mat_properties(spec=) rejects a MAT whose methods or definitions differ from the spec; drawings mapped under other MAT semantics fail the existing semantics check.
- Witnesses. No numerical witness changed: the statistics-before-optimization feature_table(..., "parallelism") witness is identical under the new defaults. statistics-0.2.0.rds cases 77-106 store serialized 0.2.x feature_specs, which the version bump rejects by design; tools/mlvref/mat_refresh_witnesses.R rebuilds each spec from its own settings and asserts the recomputed output is identical after swapping the embedded (calibration-free) spec; the only differences are mat$aof/mat$morphology and the two definition strings; all other cases byte-identical; errata line appended. Track A's junction rows are preserved (track A re-applied this script after their R007/R008 refresh).
- Parity coverage of the old paths: test-parity-performance.R (6 masks: distance, AOF, Pratt skeleton) now selects aof="mlv2023", morphology="pratt"; test-reference-io.R selects method="pratt"; both pass.

## 7. Files changed (track B)
R/mat.R, R/images.R, R/feature-spec.R, R/measurement.R; src/kernels.cpp (cpp_raster_indices rasterizer; new cpp_skeletonize_hmt, cpp_aof_centered); R/RcppExports.R, src/RcppExports.cpp (regenerated); man/compute_mat.Rd, man/skeletonize.Rd, man/feature_spec.Rd, man/trace_contours.Rd; tests/testthat/test-mat-reference.R (new: 7 blocks, 276 expectations: rasterizer reversal, hmt8 thinness/homotopy, analytical midline, D4 equivariance of the default flux on 5 masks x 8 symmetries with a thin-input precondition, method recording/spec control/version rejection, coincident points, branch-direction dependence); tests/testthat/test-parity-performance.R, test-reference-io.R (explicit legacy methods); inst/fixtures/statistics-0.2.0.rds (cases 77-106, spec definitions only); tools/mlvref/mat_* (harness); docs/evidence/mlvref/mat/: e2e_summary.csv, e2e_contour_summary.csv, equivariance.csv, forest_tie_policy.csv, forest_tie_downstream.csv, forest_stage_attribution.csv, forest_branch_invariance.csv, forest_raster_rules.csv, thinning_synthetic.png, aof_registration_synthetic.png.

### Files changed in the branch-rule round
- R/mat.R: compute_mat(branch_rule = c("invariant","mlv2023")) recorded with definition, tracing id and skeleton_ink_excluded; trace_skeleton(mat, min_points, rule) with ink exclusion and both tracers; compute_mat_properties rule-aware scoring (direction-symmetric vs forward/MLV mirror) and score images (max vs last writer); map_mat_properties rule-aware mapping; spec check on branch_rule definition; docs.
- R/mat-compat.R (new): .mlv_cleanedgelist (R reference port), .trace_edgelink, .tangent_points_mlv, .mirror_mlv.
- src/kernels.cpp (appended): cpp_edgelink_trace, cpp_cleanedgelist, cpp_fit_lines_mlv, cpp_nearest_max; KD-tree range query added inside the existing KDTree class. RcppExports regenerated.
- R/feature-spec.R: mat$branch_rule, definitions branch_rule, .branch_rule_definition(). R/measurement.R: .mat_semantics includes the branch_rule definition.
- man/compute_mat.Rd, feature_spec.Rd, trace_skeleton.Rd (roxygen).
- inst/fixtures/edgelink/{grid,fan,kite}.mask.csv and .octave.branches.csv (synthetic; generated by tools/mlvref/mat_edgelink_fixture.R + .m).
- tests/testthat/test-mat-reference.R: +5 blocks (ink exclusion; edgelink exactness; MLV mirror; reversal/renumbering invariance; D4 equivariance) and branch_rule spec/versioning assertions. tests/testthat/test-parity-performance.R: mapping oracle updated (section 6 amendment).
- inst/fixtures/statistics-0.2.0.rds: cases 77-106 spec definitions only.
- Harness: tools/mlvref/mat_forest_compat.R, mat_compat_residual.R, mat_edgelink_fixture.R/.m; mat_datasets.R configs compat/default (branch_rule passed). Evidence: forest_compat.csv; e2e_summary.csv and e2e_contour_summary.csv now hold the full 2,757-drawing run.


## 8. Tests
Full suite on the current working tree (all tracks), R CMD INSTALL --no-test-load to a scratch library, testthat::test_dir(load_package = "installed"): 105 test blocks, 2,583 expectations, 0 failed, 0 errors, 0 skipped (test-mat-reference.R: 13 blocks, 680 expectations). R CMD check not run locally; CI required.

## 9. Open questions
1. Whole-pipeline lattice equivariance is limited by sequential thinning (hmt8 = MATLAB skel). A symmetric parallel thinning would restore it at the cost of MATLAB-exact skeletons; hmt8 is kept as the default (thin, homotopic, MATLAB-exact).
2. The invariant mapping keeps R's tangent geometry and the maximum across branches; MLV's d==0 centre points are not used because they carry no tangent direction (confirmed).

Resolved: isolated closed loops (section 14), provisional IDs (renumbered), negative-distance skeleton pixels (section 10).

## 10. Branch rules and skeleton pixels on ink

*All compat runs in this report use `threshold_angle = asin(0.3*pi/2)*180/pi` (threshold 0.3; `tools/mlvref/mat_datasets.R`). Section 15 supersedes the default morphology named here.*

**Skeleton pixels on ink.** Skeleton pixels with signed distance <= 0 are ink, not medial points. trace_skeleton() removes them before tracing on every rule; the count is attribute "skeleton_ink_excluded" of the branch list and mat$parameters$skeleton_ink_excluded. This removes MLV-R012: the 7 pratt-path crashes no longer occur (no drawing failed in any run). MATLAB's Forest skeleton contains 0 such pixels, so compat reproduction is unaffected.

**branch_rule="invariant"** (default; definition "direction-symmetric-tiemax-graph-1"):
- Tracing: pixel graph (diagonal steps only where no orthogonal connection exists) decomposed into maximal paths between end/junction pixels; paths with >= min_branch_points pixels kept; junction pixels shared by meeting branches, never repeated within one.
- Scoring: each score = mean of the forward score and the score of the reversed branch mapped back (any estimator; arclength is already symmetric). Score image: maximum over branches.
- Mapping: a tangent sample of segment (i,i+1) carries the mean of the branch's own scores at i and i+1; exact nearest-sample ties (within 1e-9 px) take the maximum; inclusion d < radius - 1e-9; maximum across branches; valid zeros kept.
- Known limitation: an isolated closed loop (no end/junction pixel) is cut at a lattice-order start pixel; documented, not covered by the invariance claim.

**branch_rule="mlv2023"** (definition "mlv2023-edgelink-lastwriter-1"):
- Tracing: Kovesi edgelink (clean, skel via hmt8, junction LUT crossings >= 6, row-major raster scan, trackedge/nextpoint including the first edge labelled -0, remember-last fallback, loop closure repeating the first point) and cleanedgelist(minlength = min_branch_points, Euclidean length) with its merge/spur logic and duplicated merge points. C++ (cpp_edgelink_trace, cpp_cleanedgelist) plus an R reference port (.mlv_cleanedgelist).
- Scoring: forward mlv2023 estimator; mirror via a port of MLV fitLineSegments/getDistanceFromLineSegment on [X,Y] = (row, column), including the forced terminal span (MLV-M016). Score image: last writer wins (branch order, then point order).
- Mapping: tangent samples in MATLAB order ([FP1; FP2], segment -> r1 -> r2), coincident points contribute the skeleton point itself; knnsearch-equivalent nearest with lowest-index ties; strict d < radius; scores read from the last-writer image at the segment's first point; max across branches; zeros treated as missing (MLV-M011 reproduction).

**Spec/versioning:** feature_spec mat block gains branch_rule (default "invariant", validated); definitions gain branch_rule; MAT parameters record branch_rule, branch_rule_definition, tracing ("pixel_graph_no_corner_cut" | "kovesi_edgelink_cleanedgelist") and skeleton_ink_excluded; .mat_semantics includes the branch_rule definition; compute_mat_properties(spec=) rejects a MAT with a different branch_rule definition; specs lacking the branch_rule definition fail validate_feature_spec ("noncanonical"). Outside the mat block, the only feature-spec.R edit is one line adding branch_rule to the definitions vector next to aof/morphology (plus helper .branch_rule_definition at file end).

## 11. Published-compatible chain: validation
Oracles, all exact:
- Tracing on MATLAB's Forest skeleton: the R port returns the same 519 branches as the pinned Octave GetConSeg/edgelink, identical in order, direction, points and the 272 duplicated merge points.
- Synthetic fixtures inst/fixtures/edgelink/{grid,fan,kite}: masks are hmt8 skeletons of synthetic drawings (fixed points of hmt8); expected branches from pinned Octave GetConSeg (tools/mlvref/mat_edgelink_fixture.R/.m). R matches all three exactly (8 branches each; fan has 3 and kite 1 duplicated merge points); C++ cleanedgelist equals the R reference port.
- Mirror: .mirror_mlv reproduces Octave computeMATpropertyPerBranch mirror on the 5 parity branches within 1e-12, including terminal_bend (MLV-M016).
- Precision attribution (tools/mlvref/mat_compat_residual.R): on MATLAB's stored Forest MAT fields, the R compat chain equals the pinned Octave chain run in double at 23,567 of 23,568 stored pixels (1 pixel differs for parallelism/separation, 0 for mirror; consistent with a nearest-tangent tie). Against the Octave chain run in single (as MATLAB), 20,122 parallelism pixels differ by >1e-9 but only 24 by >1e-4. The compat residual against stored MATLAB values is therefore MLV's single-precision arithmetic (radii, tangent points), which flips a small fraction of exact nearest-tangent ties.

Forest_DBW end to end from the vector drawing (forest_compat.csv), pixels within 1e-4 (par / sep / mir):

| config | par / sep / mir | branches |
|---|---|---|
| compat (mlv2023, hmt8, mlv2023) | 0.9986 / 0.9985 / 0.9997 | 519 branches, 1.4 s |
| defaults (centered, hmt8, invariant) | 0.054 / 0.108 / 0.072 | 589 branches |
| aof=mlv2023, hmt8, invariant | 0.058 / 0.113 / 0.091 | 605 branches |

(The earlier stage-A figure 0.999 used MATLAB's fields; this row starts from the vecLD.)

## 12. End-to-end on all 2,757 drawings
This supersedes the stratified tables in section 4 as the headline; section 4's 580-drawing tables remain as the per-stage history. All drawings with stored MAT samples in the 29 exported datasets; threshold 0.3; 0 failed drawings; runtime 0.74 s (compat) and 1.39 s (default) per drawing. Receipts: e2e_summary.csv, e2e_contour_summary.csv (now the full run).

Per pixel: within 1e-4 / within 1e-2 / mean abs / per-drawing-mean r

| config | parallelism | separation | mirror |
|---|---|---|---|
| compat | 0.9909 / 0.993 / 5.2e-4 / 1.000 | 0.9944 / 0.998 / 8.6e-5 / 1.000 | 0.9931 / 0.995 / 3.6e-4 / 1.000 |
| default | 0.0836 / 0.378 / 0.064 / 0.995 | 0.1071 / 0.665 / 0.018 / 0.999 | 0.0849 / 0.324 / 0.050 / 0.991 |

Coverage of stored pixels: compat 1.000, default 0.995.

Per contour (671,629 contours): r / median abs / within 0.01 / within 0.05 / mean bias (R - MATLAB)

| config, property | r / median abs / within 0.01 / within 0.05 / bias |
|---|---|
| compat parallelism | 0.9993 / 1.0e-7 / 0.987 / 0.998 / -2.1e-5 |
| compat separation | 0.9999 / 3.2e-8 / 0.997 / 1.000 / -4.4e-6 |
| compat mirror | 0.9986 / 4.1e-15 / 0.990 / 0.997 / 6.8e-6 |
| default parallelism | 0.8906 / 0.026 / 0.264 / 0.670 / -0.023 |
| default separation | 0.9779 / 0.0075 / 0.577 / 0.905 / -0.0025 |
| default mirror | 0.8416 / 0.029 / 0.233 / 0.670 / -0.011 |

Contours defined in MATLAB but not R: compat 4, default 1,242; defined in R only: compat 5-7, default 300-364.

### 12.1 Per dataset family, all 2,757 drawings: pixels within 1e-4 / per-contour r

compat (aof=mlv2023, hmt8, branch_rule=mlv2023):

| family | n | parallelism | separation | mirror |
|---|---|---|---|---|
| beaches | 80 | 0.9916/0.9997 | 0.9942/1.0000 | 0.9929/0.9992 |
| cities | 79 | 0.9914/0.9994 | 0.9953/0.9999 | 0.9922/0.9991 |
| forests | 80 | 0.9936/0.9997 | 0.9969/1.0000 | 0.9960/0.9996 |
| high_ceiling | 100 | 0.9927/0.9994 | 0.9957/0.9998 | 0.9946/0.9988 |
| highways | 80 | 0.9911/0.9990 | 0.9945/1.0000 | 0.9925/0.9992 |
| IAPS | 1182 | 0.9905/0.9993 | 0.9937/0.9999 | 0.9938/0.9989 |
| low_ceiling | 100 | 0.9927/0.9995 | 0.9954/1.0000 | 0.9939/0.9985 |
| mountains | 76 | 0.9919/0.9998 | 0.9939/1.0000 | 0.9956/0.9996 |
| OASIS | 900 | 0.9898/0.9990 | 0.9942/0.9999 | 0.9915/0.9977 |
| offices | 80 | 0.9903/0.9992 | 0.9947/0.9999 | 0.9896/0.9986 |
| ALL | 2757 | 0.9909/0.9993 | 0.9944/0.9999 | 0.9931/0.9986 |

default (aof=centered, hmt8, branch_rule=invariant):

| family | n | parallelism | separation | mirror |
|---|---|---|---|---|
| beaches | 80 | 0.1241/0.8853 | 0.1522/0.9906 | 0.0740/0.8755 |
| cities | 79 | 0.1483/0.8472 | 0.1734/0.9751 | 0.0986/0.8308 |
| forests | 80 | 0.1182/0.9007 | 0.1447/0.9865 | 0.0799/0.8926 |
| high_ceiling | 100 | 0.1311/0.9049 | 0.1533/0.9768 | 0.0909/0.8106 |
| highways | 80 | 0.1237/0.8524 | 0.1374/0.9828 | 0.0953/0.8714 |
| IAPS | 1182 | 0.0634/0.9020 | 0.0882/0.9799 | 0.0816/0.8547 |
| low_ceiling | 100 | 0.1547/0.9147 | 0.1739/0.9782 | 0.0879/0.8291 |
| mountains | 76 | 0.0393/0.9440 | 0.0666/0.9908 | 0.0527/0.9050 |
| OASIS | 900 | 0.0337/0.8691 | 0.0574/0.9665 | 0.0783/0.7977 |
| offices | 80 | 0.1353/0.8617 | 0.1528/0.9803 | 0.1141/0.8573 |
| ALL | 2757 | 0.0836/0.8906 | 0.1071/0.9779 | 0.0849/0.8416 |

Reading: the default departs from published values by design, through four corrections (centred flux, MLV-M027; direction-symmetric scores, M028; per-branch scores and order-free ties, M029; no duplicated merge points, M030) plus R's corrected mirror fit (MLV-M016). Per-drawing means remain close (r 0.991-0.999), per-contour means moderately (0.84-0.98), per-pixel values rarely within 1e-4.

## 13. Invariance oracles for the default (all in tests/testthat/test-mat-reference.R)
- Branch reversal: reversing every branch leaves all per-pixel mapped scores (parallelism, separation, mirror, taper) unchanged within 1e-12; under branch_rule="mlv2023" the same reversal changes parallelism by >1e-3 (MLV-M028).
- Branch renumbering (random permutation) and contour renumbering (permuted contours of the drawing): per-pixel images unchanged within 1e-12; per-contour vectors permute accordingly.
- D4: for 5 synthetic drawings (parallel lines, square, plus, X, L), transforming the MAT fields by each of the 8 lattice symmetries leaves the mapped scores equivariant (max difference 0; tolerance 1e-12); the published rule shows differences up to 0.93. Whole default pipeline on a diamond drawing whose thinned flux is itself equivariant: skeleton and scores equivariant under all 8. hmt8 thinning, like MATLAB skel, is sequential and order-dependent, so whole-pipeline equivariance is not claimed for arbitrary drawings.
- Ink exclusion: skeleton pixels forced onto ink are removed under both rules and counted.

## Amendments to sections 5-8
- Section 5: MLV-R012 FIXED (ink skeleton pixels excluded, section 10). MLV-M028/B3/B4 are now both corrected by default and reproduced by branch_rule="mlv2023". The underdetermined tracing item is resolved as in section 10 (default graph decomposition; published edgelink available).
- Section 6 (witnesses): statistics-0.2.0.rds cases 77-106 rebuilt again for the branch_rule definition by the same script (every value asserted identical; cases 1-76 and 107-122 byte-identical; errata line replaced by one covering aof/morphology/branch_rule). statistics-before-optimization.rds is unchanged and its feature_table(...,"parallelism") witness still passes under the new default. test-parity-performance.R "spatially filtered MAT mapping agrees with exhaustive nearest samples": the exhaustive oracle now implements the default invariant definition (segment-mean sample scores, max over ties within 1e-9, d < radius - 1e-9); justification: the old oracle encoded the order-dependent first-point/lowest-index rule that this change replaces by design.

## 14. Final round: closed loops and the retained mapping oracle
- **Isolated closed loops** (skeleton cycles with no end or junction pixel; R/mat.R `.is_isolated_loop`, `.score_loop`). This round introduced a start/direction average; section 15 replaced it with the final definition: an isolated closed loop has no start. Under branch_rule = "invariant", local estimators (parallelism, separation, taper, and all arclength scores) are evaluated with circular windows on the periodic point sequence, averaged over both traversal directions; this is a definition, not an average over cut points. The mlv2023 mirror score, a global greedy line-segment fit without a local window, is the mean over every start pixel and both directions. Both are independent of the traced start pixel and direction. Junction-anchored loops are unaffected. Test: an elliptical ring, 3 start pixels × 2 directions, all four properties, per-pixel mapped scores identical within 1e-12; the published rule keeps MLV's start dependence (>1e-6). Full 2,757-drawing default rerun: per-pixel within 1e-4 0.0835 / 0.1071 / 0.0849 (was 0.0836 / 0.1071 / 0.0849); per-contour r unchanged at 0.8906 / 0.9779 / 0.8416. Compat results are unaffected.
- **Original mapping oracle retained:** test-parity-performance.R "spatially filtered MAT mapping agrees with exhaustive nearest samples" asserts both rules on the same synthetic branches and radii 0.1, 1, 4, 100: (a) the default invariant rule (segment-mean scores, max over ties within 1e-9); (b) branch_rule="mlv2023" with the original exhaustive first-point / lowest-index oracle, MLV tangent enumeration, last-writer rating image, strict radius, maximum across branches and zeros as missing (MLV-M011). Exact equality and identical mapping masks in both.
- **IDs in code and tests:** MLV-R010 (rasterizer), R011 (coincident tangent points), R012 (ink skeleton exclusion), R013 (Pratt replaced by hmt8), M027 (AOF misregistration), M028 (direction dependence), M029 (last-writer mapping and order ties), M030 (duplicated merge points), M031 (unbounded separation).

## 15. Independent review and fixes
A fresh-context reviewer verified five defects and four concerns in the changes above; all are fixed and tested. Full suite: 116 test blocks, 3,355 expectations, 0 failed (test-mat-reference.R: 18 blocks, 955 expectations).

1. **hmt8 is not always one pixel wide** (a solid 4×4 square leaves a 2×2 block; 16 of 40 random polygon drawings had blocks in the default skeleton). New default `morphology = "minimal"` (definition `hmt8-then-sequential-simple-pixel-removal-1`): hmt8, then sequential column-major removal, to a fixed point, of every 8-simple pixel (Yokoi number 1) with at least two neighbours. Oracles on the 4×4 square plus 40 random masks: 0 simple pixels with two or more neighbours, subset of the input, identical 8-component and hole counts, idempotent. A 2×2 block can remain only where every one of its pixels is non-simple (removing any would change the 8-connected topology); the second review found such irreducible blocks in 31 of 400 noise masks, and none in natural drawings. `morphology = "hmt8"` remains the MATLAB-exact thinning. The published-compatible recipe is `aof = "mlv2023"`, `morphology = "hmt8"`, `branch_rule = "mlv2023"`, and it still reproduces the Forest MAT exactly. Forest default skeleton: hmt8 24,885 px (307 redundant, 2 blocks) → minimal 24,539 px (0, 0).
2. **Lollipop branches** (start and end on one junction pixel) wrote two different scores to that pixel, last writer winning, so reversal changed `score_images` by 1.8e-4. The invariant rule now takes the maximum within a branch as well as across branches.
3. **Isolated-loop detection** counted degree on the unfiltered skeleton; it now uses the medial mask (skeleton and D > 0) that tracing uses.
4. **The centred AOF is lattice-equivariant only when `samples %% 4 == 0`.** `aof = "centered"` now requires it, with an informative error in `compute_mat()` and `feature_spec()`; `aof = "mlv2023"` accepts any count.
5. **`cpp_aof_centered` border rule was asymmetric** (one row skipped at the top, two at the bottom). All sample points must now lie inside the image; flip and transpose tests include border rows.

Concerns resolved: `compute_mat()` documents exactly what is invariant (the flux field for `samples %% 4 == 0`; the branch rule under reversal, renumbering, loop start and lattice transforms of the MAT fields) and that whole-pipeline equivariance is not guaranteed (sequential thinning, slope-sign rasterizer ties). Isolated-loop scoring is now O(m) via circular windows, and the global mlv2023 mirror fit is averaged over starts and directions in C++ (`cpp_mirror_loop`), equal to the brute-force average within 1e-12. `cpp_cleanedgelist` uses ordered adjacency sets and a coordinate hash (memory O(nodes + links); a 20,000-edge input runs in 0.02 s) with output identical to the R reference port. `feature_spec()` records `tracing` according to `branch_rule`.

Full 2,757-drawing rerun (receipts e2e_summary.csv, e2e_contour_summary.csv, forest_compat.csv):

| config | per-pixel within 1e-4 (par / sep / mir) | per-contour r (par / sep / mir) |
|---|---|---|
| compat (mlv2023, hmt8, mlv2023) | 0.9909 / 0.9944 / 0.9931 | 0.9993 / 0.9999 / 0.9986 |
| default (centered, minimal, invariant) | 0.0781 / 0.1047 / 0.0626 | 0.8887 / 0.9719 / 0.8359 |

Default per-drawing mean r 0.995 / 0.999 / 0.991; coverage 0.993; 0 failed drawings.

### Second review pass
The reviewer re-ran every reproducer: all five findings and four concerns fixed. The published-compatible tracer matched Octave edgelink on 238 further random skeletons (0 mismatches); the hashed cleanedgelist equalled the R reference on 160 random cases; `cpp_mirror_loop` equalled a brute-force start/direction average within 5e-15. Three minor items were then fixed:
- **N1:** isolated-loop invariance on long loops was about 2.3e-12 (1,613-point loop), above the 1e-12 claim, from running sums over the padded sequence. Loops are now scored in blocks of at most 64 points with a padding of ceiling(window) + 16 points, and `cpp_mirror_loop` uses Kahan summation: rotation and reversal now change scores by at most 7.6e-14. A test checks a loop of more than 1,500 points.
- **N2:** documentation of loop scores now states the circular-window definition above.
- **N3:** documentation no longer claims "minimal" has no 2×2 blocks; the test uses the redundant-pixel oracle on the reviewer's 400 noise masks and checks that the irreducible blocks keep the input's topology.

Full suite: 118 test blocks, 3,816 expectations, 0 failed. Non-loop branches and the compat chain are unchanged, so the receipts stand.
