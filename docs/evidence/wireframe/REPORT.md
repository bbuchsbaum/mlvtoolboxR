# Recovering junctions in Wireframe photographs

A learned T/X scorer now clears the issue's 50% recovery target for both types on the original panel. The follow-up reports its performance on 342 fresh Wireframe photographs and blinded judgments of which crossings are visible. [Trained photographic detector](TOPOLOGY.md)

The [follow-up experiment](FOLLOWUP.md) tests overlapping image crops, four-arm crossings and a second pretrained model. On 48 new photographs the combined rule recovered 49 more Ts and 19 more Xs than full-image M-LSD, with zero flags from either method at 17 reviewed gaps.

On the original 24-photo panel from issue #3, a pretrained M-LSD line-segment provider recovered **19 of 62 T junctions and 35 of 101 X junctions**, with **zero flags at the 15 reviewed gaps**. FDoG with 6 px gap repair recovered 3 Ts and no Xs. M-LSD recovered all four marked window crossings but none of the 33 eligible floor crossings. Boosting local contrast recovered six more Xs, lost five Ts and flagged one reviewed gap. Both versions fell below the protocol's 50% recovery threshold for each type. Keep the provider experimental. [Issue #3](https://github.com/bbuchsbaum/mlvtoolboxR/issues/3)

## Original benchmark

The public JSON matched the frozen archive for all 24 photo records, all 163 typed coordinates and all 15 reviewed near misses. Coordinates are one-based pixel centers on the resized images. The four archived CSVs matched the SHA-256 values recorded in the manifest. [Public benchmark manifest](https://github.com/bbuchsbaum/mlvtoolboxR/issues/3#issuecomment-7101349061)

All 24 resized RGB images matched the archived pixel hashes. The original references and review decisions were used unchanged. Reference counts, detected counts and matches for both FDoG settings agreed with the archive on **all 96 photo × setting × junction-type rows**. [Per-photo results](issue-scores.csv)

| Contour input | T matched / 62 | X matched / 101 | All detected T | All detected X | Reviewed gaps flagged / 15 |
| --- | ---: | ---: | ---: | ---: | ---: |
| FDoG, gap 0 | 0 | 0 | 40 | 4 | 0 |
| FDoG, gap 6 | 3 | 0 | 80 | 4 | 0 |
| M-LSD, rasterized and traced | 19 | 35 | 962 | 329 | 0 |
| M-LSD, original vector segments | 5 | 22 | 369 | 527 | 1 |
| M-LSD, local contrast then rasterized and traced | 14 | 41 | 902 | 355 | 1 |

The pretrained provider supplies straight lines as one-pixel strokes through the existing `trace_rgb()` interface. Keeping the exit-rule detector fixed, the traced strokes recovered 31% of the Ts and 35% of the Xs. The number of predicted junctions greatly exceeds the number of eligible reference sites. Overall precision remains undetermined because human annotations cover selected structural edges. The vector-only comparison used the identical M-LSD segments directly and flagged one frozen gap. [Results at all 15 reviewed gaps](issue-near-miss-scores.csv)

## Where the crossings are lost

**Photographic contrast.** We applied contrast-limited adaptive histogram equalization (CLAHE) to the lightness channel of each resized photo in CIELAB, using an 8 × 8 tile grid and clip limit 2. The M-LSD weights, line thresholds, rasterization, tracing, detector and matching rules stayed fixed. We chose these settings before scoring this comparison. On the original panel, the extra flag was a T near (206.74, 63.60) in 00358297. Neither version matched any of the 33 floor Xs. [Per-photo scores](issue-scores.csv) · [Gap scores](issue-near-miss-scores.csv)

We also compared the two photo providers on the 23 supplementary images outside the issue panel, excluding the shared photo 00358638. The contrast setting again increased matched Xs and added a gap flag. This panel uses its own reference rule. [Supplementary scores](scores.csv) · [Reviewed gaps](near-miss-scores.csv)

| 23 disjoint photos | T matched / 328 | X matched / 76 | Reviewed gaps flagged / 9 |
| --- | ---: | ---: | ---: |
| M-LSD | 136 | 31 | 1 |
| M-LSD with local contrast | 135 | 38 | 2 |

**Known input lines.** We bypassed photo extraction and supplied the human-annotated segments to the same junction detector. These controls measure the effect of the representation on known geometry. They use the source of the reference labels and cannot measure photographic accuracy.

| Human line input, original panel | T matched / 62 | X matched / 101 | Floor X matched / 33 |
| --- | ---: | ---: | ---: |
| Original subpixel vectors | 43 | 96 | 30 |
| Rasterized to one-pixel strokes, then traced | 26 | 87 | 29 |

The detector recovered 29 of the 33 floor crossings once the annotated lines were present in the mask. That places the main failure on this floor upstream of junction typing. For Ts, rasterization and tracing introduced a second loss. Of the 43 sites matched with the original vectors, 22 remained matched and 21 lost their T match. Four other sites gained a T. At 17 of the 21 lost sites, Y was the only junction type within the four-pixel scoring radius. [Issue results](issue-scores.csv) · [Detector rule](../../../R/geometry.R)

A three-segment example isolates the T-to-Y change. In photo 00190922, take rows 4, 5 and 7 of the prepared annotation CSV. They meet near (93.85, 56.91). On these three segments alone, the largest gap between adjacent arms is 160.67° with the original vectors and 157.91° after rasterization and tracing. That crosses the detector's 160° Y/T boundary. This case can be used to check a change in how arm directions are estimated. [Annotation source](https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test/metadata.jsonl)

## Supplementary 24-photo sample

The sample excludes both named witnesses. We ranked the remaining test photo IDs by the SHA-256 digest of `issue3-fixed-panel/` followed by the ID and took the first 24. Selection used neither labels nor predictions. No provider parameters were tuned. One photo (00358638) also appears in the original panel. The counts below use a separately defined reference rule and nine independently reviewed gaps. The two panels should be interpreted separately. [Panel and resized-pixel hashes](panel.csv) · [Per-photo counts](scores.csv)

| Contour input to the same junction detector | T matched / 350 | X matched / 80 | All detected T | All detected X | Reviewed gaps flagged / 9 |
| --- | ---: | ---: | ---: | ---: | ---: |
| FDoG, gap 0 | 4 | 0 | 49 | 5 | 0 |
| FDoG, gap 6 | 12 | 0 | 85 | 6 | 0 |
| M-LSD, rasterized and traced | 145 | 31 | 753 | 262 | 1 |
| M-LSD, original vector segments | 82 | 32 | 311 | 532 | 1 |

The rasterized provider recovered 41% of eligible Ts and 39% of eligible Xs. It also produced 262 X detections, of which 31 matched an eligible human reference. The annotations describe selected architectural lines. Unmatched detections require review before they can be counted as errors. The vector row bypasses rasterization and tracing and tests their effect on the same predictions. [Wireframe source paper](https://openaccess.thecvf.com/content_cvpr_2018/papers/Huang_Learning_to_Parse_CVPR_2018_paper.pdf)

## What disappeared in the two witnesses

| Photo | FDoG gap 0, detected T / X | FDoG gap 6, detected T / X | M-LSD traced, matched T | M-LSD traced, matched X |
| --- | ---: | ---: | ---: | ---: |
| 00079423, window bars | 0 / 0 | 2 / 0 | 1 / 1 | 6 / 7 |
| 00078318, floor tiles | 2 / 0 | 5 / 0 | 0 / 3 | 0 / 33 |

**Window.** The four marked sites agree with our derived coordinates to within 0.06 resized pixels. The FDoG crop retains pieces of the outside window frame. Both horizontal and vertical bars are absent at the four marked sites before tracing. The two gap settings produced identical masks throughout this crop. M-LSD supplies both bar directions, and the same detector assigns four distinct Xs there. The issue gives seven eligible Xs across the full photograph. Our independent eligibility rule gives six. [Window photo](https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test/images/00079423.png)

**Floor.** The four sites in the displayed crop lie on two families of tile seams. M-LSD follows the roughly upright seams, but misses the lines running across them. FDoG also loses the crosswise lines in this region. Its two gap settings again produce identical ink in the crop. Across the full photo the rasterized provider detected 22 Xs. None matched the 33 frozen reference sites. One matched under the supplementary rule, which admitted 43 sites. [Floor photo](https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test/images/00078318.png)

## Supplementary reviewed gaps

Before looking at local detections, we reviewed 18 candidates from the 24-photo sample. Candidates extended an annotated segment by 4–10 px to meet the interior of a second segment at a substantial angle. Each was at least 8 px from a derived positive reference. We took at most one candidate per photo, selected by hashing its photo and segment indices. Inspection of the RGB crop and the two source segments retained nine visibly separated pairs. The ledger includes the nine exclusions, the reasons and zero-based indices into the pinned annotations. [Review ledger](near-miss-review.csv)

A flag means that a T or X fell within 4 px of the projected intersection. The rasterized provider flagged the top edge of a lintel in 00102517. The vector method flagged the sloping trim in 00492185. In both cases the detected contours ran along the target edge, and the designated stem still ended short of it. The nine sites therefore yielded one proximity flag per M-LSD method and zero visually confirmed bridges across the designated gaps. Nine sites give a limited check on false bridging. [Results by site](near-miss-scores.csv)

## Method and reproduction

All methods used the same EXIF-corrected RGB image resized with Pillow LANCZOS to 256 px on its long side. The window pixel hash matched the published value. FDoG and the rasterized M-LSD provider passed one-pixel masks through `trace_rgb()` with `simplify = 1`. M-LSD internally resized the same pixels to its 512 × 512 model input. We used NAVER's large FP32 weights and its default line score (> 0.10) and length (> 20 heatmap pixels). The R package has no runtime dependency on this model. [Pinned input metadata](https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test/metadata.jsonl) · [Model and reference inference](https://github.com/navervision/mlsd/tree/453cafa09467d0272760578d35c1fda38e8895a5)

For the supplementary sample, reference candidates were human segment endpoints and strict interior intersections. Coordinates follow resized pixel centers. The scorer requires at least two source segments within half a pixel, incident rays extending at least 4 px, a 4 px border margin and a minimum 30° angular separation after merging directions within 10°. A T requires one opposite ray pair within 20° of straight, an X two such pairs. Same-type sites within 1 px are deduplicated. Matching maximizes the number of same-type pairs within 4 px, then minimizes their total distance. [Scorer](../../../tools/wireframe/score.py)

To reproduce from the repository root, use the package source at commit `38f0e26` with its R dependencies, including `png` and `pkgload`. The run used R 4.5.1 and Python 3.12.11. Keep model weights and images in a local working directory.

```sh
work=$(mktemp -d)
uv venv --python 3.12.11 "$work/venv"
uv pip install --python "$work/venv/bin/python" Pillow==12.1.1 numpy==2.5.5 \
  opencv-python-headless==5.1.0.95 ai-edge-litert==2.4.0 scipy==1.18.1
curl -fsSL https://api.github.com/repos/bbuchsbaum/mlvtoolboxR/issues/comments/7101349061 |
  jq -er '.body | split("```json\n")[1] | split("\n```")[0]' > "$work/issue-3.json"
printf '%s  %s\n' a65df9d883d62e611a6affadac0e0929202db7e56193a8879e9797a751301f14 "$work/issue-3.json" |
  shasum -a 256 -c -
"$work/venv/bin/python" tools/wireframe/prepare.py "$work" --manifest "$work/issue-3.json"
Rscript --vanilla -e 'pkgload::load_all(".", quiet = TRUE); source("tools/wireframe/compare.R")' "$work"
"$work/venv/bin/python" tools/wireframe/score.py "$work" docs/evidence/wireframe/near-miss-review.csv --manifest "$work/issue-3.json"
"$work/venv/bin/python" tools/wireframe/witnesses.py "$work"
```

The original panel results are in `issue-scores.csv`, `issue-sites.csv` and `issue-near-miss-scores.csv`. The supplementary sample uses `scores.csv` and `sites.csv`. The witness figure is `witness-inspection.png`. The package's line-extraction tests passed, including the #2 occlusion repair, near-miss controls, and filled-grid versus centerline-grid distinction. [Regression tests](../../../tests/testthat/test-lines.R)
