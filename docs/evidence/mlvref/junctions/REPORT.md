# Track A: geometry and junctions against real MATLAB MLV outputs

Plan: `docs/plans/mlv-reference-qualification.md`, phases P1 (version drift) and P2 (geometry and junctions). Date 2026-11-19. Base commit 5fa7fa5 plus the working-tree changes listed at the end. Runtimes: R 4.5.1 (aarch64 macOS), GNU Octave 11.1.0 with image 2.20.0, pinned MLV MATLAB source f12418a (`vendor/MLV_toolbox/MLVcode`). Natural data are not redistributable and are not committed; this directory holds only aggregate receipts and the synthetic suite.

## Summary
- **No version drift.** The stored MATLAB fields are what the pinned source computes.
  - Octave re-execution of 120 stratified drawings reproduces every stored orientation, length and curvature (max |Δ| 5.7e-14 deg; relative 1.6e-14) and all 14,842 junctions.
  - An independent R transcription reproduces all 395,134 stored junctions in all 3,017 drawings.
- **MLV-compatible rule.** `detect_junctions(rule = "mlv2023")` ports the validated transcription into the package.
  - It reproduces 395,125 of 395,134 stored junctions exactly: position, type, ordered contour:segment list including MLV's appended segments, every angle, and min/max angle. Junction counts are equal (395,134).
  - All 9 differences lie in the 401 drawings that contain zero-length or rounding-residue segments, which an `mlv_drawing` cannot represent (MLV-M015/R009).
  - On those 401 drawings, `rule = "mlv2023"` equals the pinned algorithm run on the cleaned geometry for 64,208/64,208 junctions.
  - On all 41 synthetic cases it reproduces pinned Octave exactly; this is a committed test.
- **Default rule.** The default is `rule = "exit"`, `gap_metric = "euclidean"`.
  - Exact agreement with stored MATLAB: 345,018 (87.3%) at the base commit; 356,113 (90.1%) under the new defaults.
  - Every difference is attributed to one named cause: a documented MLV difference, the documented R definition, M018 (the gap metric, now a deliberate default), or one of 5 floating-point ties at the 2 px radius.
  - R equals an emulation of its documented definition except for 5 floating-point ties at the 2 px radius (389,974 junctions).
- **Segment geometry.** R matches the stored values within 1e-9 on all 2,284,239 segments, except:
  - 571 segments R drops: 552 exact zero-length (M015) and 19 rounding-residue segments (R009);
  - 99 curvatures next to those segments.
- **R-port revisions:** MLV-R005 to R009, and MLV-R014 to R017 after independent review.
- **New MLV findings:** M018 (anisotropic gap, a contract clarification) and M019–M026 (executed differences). Re-checked: M006, M007, M008, M015.
- **Synthetic suite (41 cases):**
  - default rule matches the hand-derived answer in 41/41;
  - pinned MLV under Octave in 24/41 (`gap_diag_1.2` is counted as a mismatch under the Euclidean default);
  - `rule = "mlv2023"` matches pinned MLV in 41/41.
- **Tests:** 3,355 expectations on the current working tree, 0 failures, 0 errors. `tests/testthat/test-junction-oracle.R` holds 1,049 of them.

## Method
**Inputs.** `tools/mlvref/export.m` (P0) wrote the stored vecLD fields as %.17g CSV for all 30 datasets: 3,017 drawings in IAPS (1,182), OASIS (900), Architectural (200), Snodgrass–Vanderwart (260) and Toronto scenes (475).

**Parsing caveat (harness).** R's decimal parser (`R_strtod`; `read.csv`, data.table and vroom behave the same here) is not correctly rounded for 17-digit input.
- Example: "870.99999999999977" parses one ulp from the nearest double.
- The harness re-parses every numeric field with the C library `strtod`, which is correctly rounded.
- With R's parser, two segments of length ~1e-13 looked zero-length and two orientations looked like drift. Both were parsing artifacts.

**Tolerances.** Positions ≤1e-6 px; angles and min/max angles ≤1e-6 deg; geometry ≤1e-9 (relative for curvature).
- Coordinates are integers or subpixel doubles, and both implementations evaluate the same rational expressions in IEEE double, so honest agreement is at about 1e-12.
- Junctions are paired greedily by distance within 2 px, for diagnosis only.

**Tools** (all in `tools/mlvref/`):
- `drift_pinned.m`, `drift_compare.R`: P1.
- `compare_junctions.R`: R defaults vs stored. Zero-length test: the constructor's 64-ulp rule.
- `compare_mlv2023.R`: `rule = "mlv2023"` vs stored, all drawings → `mlv2023-vs-stored.csv`.
- `mlv_junction_emulator.R`: independent transcription with one switch per cause; `arms = "exit"` and `gap = "euclidean"` emulate the package default.
- `attribute_junctions.R`: walks from the stored output through the emulator to R, one cause at a time → `attribution-*-final4.csv`, `residual-features-final4.csv`.
- `synthetic_junctions.R`, `synthetic_junctions.m`: synthetic suite. Also writes the committed fixtures `inst/fixtures/octave/junction_cases.csv` and `junction_octave.csv`: pinned-MLV outputs on synthetic geometry, with zero-length segments removed.
- `refresh_junction_witnesses.R`: surgical fixture refresh.
- `feature-spec-junction.patch`: patch for `R/feature-spec.R` (see Files changed).

Reproduce (repository root, after export):
```sh
octave-cli --quiet --no-history tools/mlvref/drift_pinned.m 4 1500
Rscript tools/mlvref/drift_compare.R <lib>
Rscript tools/mlvref/compare_junctions.R <lib> stored 7 -after
Rscript tools/mlvref/attribute_junctions.R <lib> stored 7 -final4
Rscript tools/mlvref/compare_mlv2023.R <lib> stored 7
Rscript tools/mlvref/synthetic_junctions.R <lib>
```

## P1: version drift
Receipt `drift-pinned-vs-stored.csv`: 120 drawings, 87,961 segments, 1,822 s of Octave time.

| Field | Equal within 1e-9 | Bitwise equal | Max difference |
|---|---|---|---|
| orientation | 87,961/87,961 | 75,234 | 5.7e-14 deg (last-place libm atan2d) |
| length | 87,961/87,961 | 87,961 | 0 |
| curvature | 87,961/87,961 | 70,704 | 1.6e-14 relative |
| junctions (position, type, contour:segment list, angles, min/max) | 14,842/14,842 | – | 0 |

Full-corpus corroboration: the faithful emulation (V0) reproduces 395,134/395,134 stored junctions. Every stored field can therefore support claims about the pinned code; none needs "published but drifted" status.

## MLV-compatible rule (`rule = "mlv2023"`)
`R/junctions-mlv2023.R` ports the transcription that reproduces all stored junctions, including every recorded difference:
- Chebyshev gap with contours shorter than AE skipped;
- 1e-4 endpoint snapping (M019) and discarded negative coordinates (M008);
- upper-triangular recursive merging (M006);
- start-first attachment with one adjacent segment (M022, M024);
- interior arms measured from the centroid with the y-mirrored atan2 (M007, M021), and two arms for any distant member (M025);
- no closed-contour wrap (M026); 2-arm removal.

`gap_metric = "euclidean"` together with `rule = "mlv2023"` is an error. `method` is ignored, and the members list includes MLV's appended segments.

Verification (`mlv2023-vs-stored.csv`):
- 3,017 drawings; stored 395,134, R 395,134.
- Exact: 395,125. The 9 others all lie in the 401 drawings with zero-length or near-zero segments (8 of them are MLV-M020).
- On those drawings, R equals faithful MLV run on the cleaned geometry for 64,208/64,208 junctions. The remaining gap is the input representation (M015/R009), not the rule.
- Synthetic: 41/41 exact against Octave. Zero-length cases are compared on cleaned geometry; this is asserted in `test-junction-oracle.R`.

## P2: R against stored MATLAB

### Segment geometry (`summary-stored-after.csv`)
- All compared segments agree within 1e-9.
- R drops 571 segments: 552 exactly zero-length (M015) and 19 whose endpoints differ only by rounding residue, at most 64 ulp of their own magnitude (MLV-R009).
  - Example of the latter: length 2.5e-13 at coordinate 871, where MATLAB reports curvature 2.1e13 deg/px.
- 99 curvatures differ next to a dropped segment: MLV measures the turn against the degenerate neighbour.
- There are no orientation mismatches.

### Junction attribution (`attribution-total-final5.csv`)
Each step changes one cause. "Changed" = max(n_A, n_B) − exact matches.

| Step | Cause | Changed | Junctions after |
|---|---|---|---|
| S>V0 | faithful emulation of pinned MLV | 0 | 395,134 |
| V0>V1 | interior-arm y convention changed (M007) | 24,759 | 395,134 |
| V1>V2 | transitive merging (M006) | 462 | 394,910 |
| V2>V3 | keep negative-coordinate intersections (M008) | 0 | 394,910 |
| V3>V4 | no 1e-4 snapping (M019) | 29 | 394,911 |
| V4>V5 | no short-contour exclusion | 0 | 394,911 |
| V5>V6 | zero-length and residue segments removed (M015/M020/R009) | 9 | 394,911 |
| V6>V7 | exit arm definition (R005–R008, R014–R015 vs M021–M026) | 9,554 | 394,697 |
| V7>V8 | Euclidean absolute gap (M018; package default) | 7,742 | 389,975 |
| V8>R | `detect_junctions()` defaults | 5 (ulp ties at the 2 px radius) | 389,974 |

Changed junctions by family (per-dataset sums; `attribution-by-dataset-final4.csv`; per-step totals updated in final5, per-family figures from final4):

| Family | Stored | M007 | M006 | M019 | M015/R009 | Arm definition | Euclidean gap |
|---|---|---|---|---|---|---|---|
| IAPS | 128,010 | 518 | 24 | 0 | 0 | 991 | 348 |
| OASIS | 119,248 | 221 | 148 | 0 | 0 | 3,774 | 260 |
| Architectural | 51,776 | 76 | 41 | 0 | 0 | 1,633 | 82 |
| Snodgrass | 5,917 | 3,676 | 0 | 10 | 9 | 127 | 128 |
| Toronto | 90,183 | 20,268 | 249 | 19 | 0 | 2,870 | 6,927 |

### R defaults against stored, before and after
Receipts: `summary-stored-before.csv` (base commit), `summary-stored-after.csv` (`rule = "exit"`, `gap_metric = "euclidean"`).

| Family | Exact before | Exact after | MLV-only before → after | R-only before → after |
|---|---|---|---|---|
| IAPS | 125,540 | 126,439 | 13 → 199 | 0 → 44 |
| OASIS | 114,975 | 115,115 | 85 → 327 | 0 → 39 |
| Architectural | 49,995 | 50,034 | 22 → 101 | 0 → 15 |
| Snodgrass | 1,271 | 2,084 | 0 → 32 | 0 → 1 |
| Toronto | 53,237 | 62,510 | 153 → 4,687 | 4 → 85 |
| **All** | **345,018 (87.3%)** | **356,113 (90.1%)** | **273 → 5,344** | **4 → 184** |

- Per-family rows are from the final4 run; totals are from final5. The default rule now matches stored MATLAB on 69 fewer junctions, from the R014/R015 rule changes; the split between the two was not measured.
- The large MLV-only count is mostly M018. Under the Euclidean default, diagonal gaps of 1–1.41 px are no longer joined; that is 4,725 fewer junctions (Toronto 4,438, IAPS 181, OASIS 62, Snodgrass 31, Architectural 13; step V7>V8).
- The rest are MLV junctions with fewer than three arms leaving the neighbourhood (M023/M024).
- R-only junctions involve closed contours that MLV does not wrap (M026). Gap-metric changes add a few.
- With `gap_metric = "chebyshev"`, the exact count is 362,075 (91.6%), 601 MLV-only and 164 R-only.

Arm-definition class, V6 vs V7 (from the final4 run; final5 adds R014/R015), with features in `residual-features-final5.csv`:

| Sub-class | Junctions | Oracle adjudication |
|---|---|---|
| MLV has more arms (1,255 of them include an exact duplicate direction) | 7,484 | Oracle agrees with R: M022/M023/M024/M025 |
| R has more arms (closed-contour wrap) | 68 | Oracle agrees with R: M026 |
| Same count, interior arm measured from an offset merged centroid | 960 | Oracle agrees with R: M021 |
| Same count, arm leaves via a later segment than MLV's incident one | 444 | Underdetermined at sub-radius scale; R definition documented |
| MLV junctions with <3 arms leaving the neighbourhood | 376 | Oracle agrees with R: M023/M024 |
| R junctions MLV lacks (closed contours) | 164 | Oracle agrees with R: M026 |

The commonest transition is MLV Star (5 arms) → R X (4 arms), 5,659 times. In these, one contour passes the crossing with a vertex within 2 px, and MLV counts both of its segments as forward arms (M024).

Natural effect of R007/R008 alone (final2 → final3): 419 junctions changed, 247 of them involving a closed contour (R007) and 172 removals of an additional arm (R008). R's junction count went from 394,529 to 394,697.

## Adjudications
Each class lists its oracle and minimal reproducer. Every reproducer is a case in `synthetic_junctions.R`; the R behaviour is asserted in `test-junction-oracle.R` or `test-geometry.R`.

### Re-verified, already recorded
- **MLV-M006 (merging):** 462 natural junctions; 224 fewer junctions after transitive merging.
  - Oracle: renumbering contours must not change the result.
  - `merge_chain`: MLV gives an X at (13,10) plus a Star at (10.75,10). After renumbering it gives one Star at (11.5,10). R gives one Star at (11.5,10) for every order.
- **MLV-M007 (interior y convention):** 24,759 natural junctions (6.3%). Dominant in Toronto and Snodgrass, where subpixel diagonal crossings are common.
  - `T_diagonal`: MLV min angle 0; analytical 90.
- **MLV-M008 (negative coordinates):** 0 natural junctions.
  - Synthetic: MLV finds none; R finds the analytical X at (−3,−3).
- **MLV-M015 (zero-length segments):** 552 segments, 81 curvatures, and 8 junctions via M020.

### New MLV findings (proposed register entries)
- **MLV-M018 (contract clarification only).**
  - Finding: lineIntersection.m converts AE to a parameter fraction using max(|dx|,|dy|). The absolute gap is therefore Chebyshev, up to √2 px along a diagonal, though the help text says "1 pixel".
  - Trigger: a 45° stem 1.2 px (Euclidean) short of a bar is joined; the same stem rotated to 0° is not.
  - Natural effect: a Euclidean rule would change 7,768 junctions (4,731 fewer).
  - R: default `gap_metric="chebyshev"` keeps the MLV semantics (one raster step); `"euclidean"` is rotation-invariant. Both are tested.
- **MLV-M019, endpoint snapping.** When the intersection parameter is within 1e-4 of 0 or 1, MLV returns the segment endpoint, a displacement of up to 1e-4 × length.
  - Bar (0,5)–(40000,5) crossed at x=3: MLV gives Y at (0,5); analytical X at (3,5).
  - 29 natural junctions. R computes the exact intersection.
- **MLV-M020, zero-length neighbour arm.** A zero-length neighbour adds an arm of orientation atan2d(0,0)=0.
  - MLV gives X [0 90 90 180]; analytical T [90 90 180].
  - 8 natural junctions. R drops zero-length segments explicitly.
- **MLV-M021, interior arms from the merged centroid.** MLV measures interior arms from the merged centroid, which shifts their angles by an amount set by offset/length.
  - Three lines crossing pairwise within 1.5 px: MLV gaps 43.5–46.4° and 88.1/91.9°; analytical exactly 45° ×4 and 90° ×2.
  - 960 natural junctions. R takes the arm direction from the leaving segment.
- **MLV-M022, start-first attachment.** MLV attaches a segment at its start whenever the start is within 2 px, even when the end is nearer, so reversing a contour changes the result.
  - `stub_stem_reversed`: a stem 0.2 px short of a bar with a 1 px first piece is a T drawn one way and an X drawn the other; analytical T both ways.
  - Natural: part of the 7,484.
- **MLV-M023, splitting creates junctions.** Splitting one polyline into contours with a short middle piece creates a junction.
  - (0,10)–(10,10)–(11,11)–(20,11) as three contours: MLV gives T at (10.5,10.5); analytical: a continuation with no junction.
  - Natural: part of the 7,484, plus the 376 MLV-only junctions.
- **MLV-M024, duplicate arms.** Every segment of one contour that lies inside the neighbourhood contributes an arm.
  - Straight line with vertices 1 px either side of a crossing: MLV gives a Star with a 0° gap; analytical X [90×4].
  - 1 px overshoot: MLV gives X [0 90 90 180]; analytical T.
  - Natural: 1,255 exact duplicates, and most of the 5,659 Star→X cases.
- **MLV-M025, additional arm.** A member segment with both endpoints more than 2 px away is always treated as crossing (two arms), even when it lies wholly beyond the junction.
  - Merged chain at default settings (bar; crossings at x=10, 11.9, 13.8; T stem from (15.7,10)): MLV gives 10 arms with gaps shifted by M021; analytical 9 arms, gaps [0×5, 90×4].
  - Natural: part of the 7,484.
- **MLV-M026, closed contours.** MLV treats a closed polygon as ending at its first vertex.
  - Bar plus the polygon (10,11),(10,0),(30,0),(20,20): MLV gives T started at (10,11), X from any other start. Analytical X with gaps [pb, 90, 90, 180−pb], where pb = atan2(9,10)·180/π.
  - Natural: 68 R-more-arms plus 164 R-only junctions.

### R-port revisions
- **MLV-R005, arm direction.** The arm direction was the chord from the merged junction position to the first vertex outside the radius. Straight-arm angles were biased by up to 0.4° in the offset-rays case. This was the dominant base-version class: 28,975 junctions.
  - Oracle: for straight arms, the angles between the lines do not depend on where the merged position falls.
  - Fix: the arm direction is the orientation of the segment through which the contour leaves the neighbourhood.
- **MLV-R006, arm membership.**
  - Failure: a segment was attached at its nearer endpoint, and ties went to the start, so reversing a contour could change the type. Arms toward contour ends inside the radius were counted, so overshoots and split polylines became junctions.
  - Oracle: reversal and contour-partition invariance; a straight line passing a crossing has two arms.
  - Fix: arms are the directions in which each incident contour leaves the open disc of radius `endpoint_radius`. A contour that ends inside the disc contributes nothing in that direction.
- **MLV-R007, closed contours** (found by the independent reviewer). Closed contours did not wrap, so the result depended on the start vertex.
  - Oracle: start-vertex and reversal invariance (4 rotations × 2 directions).
  - Fix: indices wrap when the last end equals the first start; a closed contour wholly inside the disc gives no arm. 247 natural junctions.
- **MLV-R008, additional arm** (found by the reviewer). A member segment that never enters the disc and lies wholly beyond the junction got an arm pointing back across it. A gapped T with gap > `endpoint_radius` read as X, and merged chains gained arms.
  - Oracle: a gapped T has exactly one stem arm.
  - Fix: if the segment does not reach the disc and the junction's projection lies at or beyond an endpoint, only the arm leading away is counted. 172 natural junctions.

- **MLV-R009, rounding-residue segments.**
  - Failure: the constructor rejected only exactly zero-length segments. 19 natural segments of length about 1e-13, the residue of earlier coordinate transforms, were accepted and produced curvatures of 1e13 deg/px.
  - Oracle: endpoints that differ by at most 64 ulp of their own magnitude are the same point to working precision.
  - Fix: `mlv_drawing()` treats such segments as zero-length, under the `zero_length` policy (error or drop). There is no absolute floor, so genuine small-scale geometry (e.g. a 1e-15 segment at the origin, already covered by `test-curvature-estimators.R`) stays valid.
  - Tests: residue segment rejected or dropped; 1e-6, 1e-15 and 1e-10-at-871 segments kept.

- **MLV-R014, merged neighbourhood** (found by the independent review).
  - Failure: the neighbourhood was one disc around the merged centroid. A merged crossing far from the centroid was handled inconsistently; a member whose projection lies inside the segment but which never enters the disc still got two arms.
  - Oracle: merging is bookkeeping. The merged junction's arms must be the union of the arms its crossings have when unmerged; per contour exit, a shared bar counts once.
  - Fix: the neighbourhood is the union of open discs of radius `endpoint_radius` around every merged intersection. The R008 "beyond the junction" test uses the segment's own intersections; the roxygen states the rule.
  - Test: five crossings 2 px apart (outer stems 4 px from the centroid) give one Star with bar 0/180 and five stems 90/270 (gaps 0×8, 90×4). With merge_radius = 1 they are five X.
- **MLV-R015, intersection acceptance** (found by the independent review).
  - Failures: the sweep box was not conservative (`sg(0,0,0,6)` + `sg(1,3,4.5,3)`: 0 junctions by sweep, 1 by exhaustive search; the gap equals the limit), and an extension exactly at the limit depended on segment direction.
  - Fix in `cpp_intersections`: each segment is evaluated in a canonical frame (lexicographically smaller endpoint, direction toward the other); a pair is evaluated in geometric order; the extension beyond either end is measured from that end with the same arithmetic, and equality is accepted; the reported point is computed from the canonical first segment; the candidate box gets a relative pad (1e-9·(1+|coordinate|)), which only filters.
  - Result: sweep and exhaustive search are bitwise identical, and acceptance is exactly reversal-symmetric. Tests: both repros; 150 half-pixel-grid drawings checking sweep == exhaustive with expect_identical (both metrics) and reversal invariance. Uncommitted check: 6,000 random half-pixel-grid drawings, 0 differences between methods, 0 under reversal, 0 under renumbering.
- **MLV-R016, residue drops** (found by the independent review).
  - Failure: dropping a residue segment (R009) left a residue gap. A dropped closing segment opened a closed contour (T vs X by start vertex); consecutive residue segments accumulated beyond 64 ulp and failed the connectivity check.
  - Fix: when segments are dropped, the next kept segment starts exactly where the previous kept one ends, cyclically for contours that were closed. Connectivity and closure stay exact.
  - Tests: the square in two start orders is exactly closed and connected with the same X at (10,15); the accumulated-residue contour cleans to (0,1)–(1,1)–(5,1).
- **MLV-R017, merge order in `rule = "mlv2023"`** (found by the independent review).
  - Failure: recursiveClasses visited merge partners in x-sorted sweep order; MATLAB `find` uses ascending index (contours `1,2,4,5,3,5,5` vs Octave's `1;2;3;5;4;5;5`, and different angles).
  - Fix: partners sorted ascending; the recursion is iterative (explicit stack, same pre-order), so long chains cannot exhaust expression depth.
  - Tests: the repro is asserted exactly against Octave and added to the synthetic suite and committed fixture (41 cases, 41/41 exact). A 3,000-crossing chain gives one junction with 6,002 arms under both rules. The same ordering also occurred in the harness emulator and was revised there; the stored corpus never exercised it (the faithful emulation still reproduces 395,134/395,134).

### Underdetermined, documented definitions
- **Sub-radius arm direction (444 junctions):** MLV uses the incident segment; R uses the first segment that leaves the disc. Both are scale choices.
- **Pieces shorter than the radius:** a contour lying wholly inside the 2 px disc contributes no arm in R. MLV counts it, but skips contours shorter than AE = 1 px. Part of the 7,484.
- **Gap metric (M018):** decided. Euclidean is the default, so the absolute gap is a rotation-invariant length. `gap_metric = "chebyshev"` (one raster step, as in MLV) stays available, and `rule = "mlv2023"` always uses it. Effect: 7,745 junctions change (4,725 fewer), 6,927 of them in Toronto.
- **Short-contour exclusion:** 0 natural effect, because every contour shorter than 1 px is zero-length and is dropped.
- **Boundary ties (5 junctions):** V8 vs R. They flip under a ±1e-9 change of the merge or endpoint radius.

## Synthetic suite (`synthetic-suite.csv`)
41 cases, each with a hand-derived answer (the 41st, `merge_order`, has analytical default answer none: every small contour lies inside the 2 px neighbourhood, so only the bar's two arms leave):
- basic types: L, T, diagonal T, X, oblique X, Y, Y from shared ends, Arrow, Star, X from shared ends, ψ;
- gaps and continuations: gaps of 0.9/1.1 px, diagonal gaps of 0.9/1.2 px, collinear continuation and overlap, T at a vertex, bent bar;
- merging: merge chain in two orders, merged chain with a T stem;
- MLV-specific triggers: negative coordinates, long-segment snapping, zero-length neighbour, offset rays, offset lines;
- the neighbourhood rule: short contour inside the neighbourhood (both directions), long crossing at a bent vertex, split polyline, overshoot, stub stem (both directions), vertices near a crossing;
- closed polygon at 4 start vertices, and reversed stem.

For each case the suite records three results:
- the default rule against the analytical answer;
- pinned MLV under Octave against the analytical answer;
- `rule = "mlv2023"` against pinned MLV, exactly: junction set, type, position, ordered angle sequence, ordered contour and segment lists.

Results:
- **Default rule:** 41/41.
- **Pinned MLV:** 24/41. The 17 mismatches (the 41st case is MLV-R017's order case, where MLV reports a junction):
  - M006 ×2, M007, M008, M018 (`gap_diag_1.2`, under the Euclidean definition), M019, M020, M021;
  - M022, M023, M024 ×4, M025, M026.
- **`rule = "mlv2023"`:** 41/41. For `zero_length_neighbour`, the reference is Octave on the cleaned geometry, the only input an `mlv_drawing` can hold. The uncleaned Octave run documents M020.

The committed fixtures `inst/fixtures/octave/junction_cases.csv` and `junction_octave.csv` pin `rule = "mlv2023"` in `test-junction-oracle.R`.

## Tests
- **`tests/testthat/test-junction-oracle.R`: 1,049 expectations.** All default-rule expectations are hand-derived. They cover R005–R009, the Euclidean-default rotation invariance and the explicit Chebyshev option (M018), and gap > `endpoint_radius`. They also cover the closed polygon in all 4 start vertices × 2 directions, and sweep == exhaustive on 60 random drawings for both metrics.
- **The `rule = "mlv2023"` block** pins all 41 synthetic cases to the Octave fixture, and checks that combining `"mlv2023"` with `"euclidean"` is an error.
- Additional random checks, not committed: sweep == exhaustive on 600 drawings for both metrics, and 200 drawings invariant under contour reversal plus renumbering. Both gave 0 differences.
- **Witness fixtures** (`statistics-0.2.0.rds`, `statistics-before-optimization.rds`): `refresh_junction_witnesses.R` replaced only `junction_angle`/`junction_type` rows (108 values per file for R005–R008). It stops if anything else differs.
  - The drawings that lost junctions are random subpixel scribbles 1–4 px across, which lie inside the 2 px disc.
  - The Euclidean default changed no witness values (0 rows).
  - Track B's `mat_refresh_witnesses.R` (embedded specifications only, values unchanged) composes with this refresh.
- **Full suite on the current working tree:** 19 files, 116 tests, **3,355 expectations, 0 failures, 0 errors**.
- **With `tools/mlvref/feature-spec-junction.patch` applied** in a scratch copy, followed by `mat_refresh_witnesses.R` and `refresh_junction_witnesses.R` (0 changes): 19 files, 100 tests, **2,182 expectations, 0 failures, 0 errors**.

## Files changed (track A)
- `R/geometry.R`: exit arm definition (R005–R008); `gap_metric` with default `"euclidean"`; `rule`; documentation.
- `R/junctions-mlv2023.R` (new): the MLV-compatible rule.
- `R/objects.R`: MLV-R009 and constructor documentation.
- `src/kernels.cpp`: `gap_norm` in `cpp_intersections` only.
- `R/RcppExports.R`, `src/RcppExports.cpp`: regenerated.
- `man/detect_junctions.Rd`, `man/mlv_drawing.Rd`: regenerated.
- `inst/fixtures/statistics-0.2.0.rds`, `statistics-before-optimization.rds`: junction rows only.
- `inst/fixtures/octave/junction_cases.csv`, `junction_octave.csv` (new, synthetic).
- `tests/testthat/test-junction-oracle.R` (new).
- `tools/mlvref/`: `drift_pinned.m`, `drift_compare.R`, `compare_junctions.R`, `compare_mlv2023.R`, `mlv_junction_emulator.R`, `attribute_junctions.R`, `synthetic_junctions.R`, `synthetic_junctions.m`, `refresh_junction_witnesses.R`, `feature-spec-junction.patch` (all new).
- `feature-spec-junction.patch` adds `junction$rule` (default `"exit"`) and `junction$gap_metric` to `R/feature-spec.R` (`"chebyshev"` for `"mlv2023"`, otherwise `"euclidean"`; `"mlv2023"` with `"euclidean"` is an error), a definition entry `junction = "exit-neighbourhood-1-<metric>"` or `"mlv2023-port-f12418a"` so specifications without it are rejected as noncanonical, and `tests/testthat/test-junction-spec.R`. `R/measurement.R` needs no change.
- **Authoritative receipts:** `drift-pinned-vs-stored.csv`, `summary-stored-before.csv`, `summary-stored-after.csv`, `attribution-total-before.csv`, `attribution-by-dataset-before.csv`, `attribution-total-final5.csv`, `attribution-by-dataset-final5.csv`, `attribution-by-dataset-final4.csv` (per-family table), `residual-features-final5.csv`, `mlv2023-vs-stored.csv`, `synthetic-suite.csv`.

## Open questions for the maintainer
1. **Zero-length segments under `rule = "mlv2023"`.** Nine stored junctions (M020) depend on zero-length segments, which an `mlv_drawing` cannot hold. Representing them would weaken the constructor invariant, so `rule = "mlv2023"` reproduces MLV on the cleaned geometry; this is the documented boundary of the compatibility rule.

Resolved: arm-rule compatibility option (added), gap default (Euclidean), near-zero segments (R009), drawing count (3,017; the plan's figure of 3,357 has been revised).
