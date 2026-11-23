# Photographic T and X recovery with overlapping crops and line fragments

On 48 previously unused Wireframe photos, the full-image M-LSD baseline recovered 329 of 864 reference Ts and 65 of 194 Xs. Adding junctions supported by overlapping crops or four visible line arms raised the counts to **378 Ts and 84 Xs**. Both methods flagged **0 of 17** reviewed gaps. The combined rule produced 437 more detections for 68 more reference matches. On the original 24-photo issue panel it recovered **20/62 Ts and 42/101 Xs**, still below the 50% target for each type. Use this as an experimental photo detector. [Frozen issue manifest](https://github.com/bbuchsbaum/mlvtoolboxR/issues/3#issuecomment-7101349061) · [Full results](support-summary.csv)

## What improved

Each cell shows baseline → combined. A reference match requires a unique same-type detection within four resized pixels. A gap is flagged when either type falls within four pixels of its reviewed location.

| Photos | Matched T | Matched X | All T + X detections | Reviewed gaps flagged |
| --- | ---: | ---: | ---: | ---: |
| 23 used for development | 136 → 156 / 328 | 31 → 31 / 76 | 962 → 1,146 | 1 → 1 / 9 |
| 24 frozen issue photos | 19 → 20 / 62 | 35 → 42 / 101 | 1,291 → 1,540 | 0 → 0 / 15 |
| 24 additional photos used to compare alternatives | 208 → 228 / 496 | 35 → 42 / 93 | 1,399 → 1,665 | 3 → 3 / 16 |
| 48 first scored after the combined rule was fixed | 329 → 378 / 864 | 65 → 84 / 194 | 2,387 → 2,824 | 0 → 0 / 17 |

The 48-photo gain was 5.7 percentage points for Ts and 9.8 for Xs. T matches increased on 27 photos and X matches on 11. Every baseline detection was retained. Overall reference matches per detection were 16.5% before and 16.4% after. The human annotations cover selected structural lines, so these fractions cannot establish precision. On the original panel this fraction fell from 4.2% to 4.0%. [Per-photo results](support-by-photo.csv) · [Results at every reviewed gap](support-gaps.csv)

To check whether the gains follow from denser output alone, we shifted just the added detections to other positions within each photo 199 times. A common wrapped translation preserves the number of additions, their types and their relative locations. On the 48 photos, the 95th percentile of the gain was **5 Ts and 2 Xs** after shifting, compared with **49 Ts and 19 Xs** at the measured locations. On the issue panel the corresponding numbers were 1 T and 1 X after shifting, against 1 T and 7 Xs observed. A second maximum-matching implementation independently agreed on every per-photo match count for the baseline and combined rule. [Shift results](shift-check.csv) · [Independent check](../../../tools/wireframe/check_localization.py)

## What each source contributes

Both alternatives start from the previously tested M-LSD line segments rasterized and traced through the package. The **crop-only** variant uses the same contour tracer and exit-rule junction detector for every added point. On the 48 photos it recovered **383 Ts and 79 Xs**. The **combined** variant trades five of those T matches for five X matches. Its extra four-arm Xs are inferred directly from photographic line fragments. Its output is a table of typed points. It does not construct a single new contour drawing. [Method comparisons](support-summary.csv)

| Original issue panel | Matched T / 62 | Matched X / 101 | Floor X / 33 | Gaps flagged / 15 |
| --- | ---: | ---: | ---: | ---: |
| Full image | 19 | 35 | 0 | 0 |
| Add crop-supported points only | 21 | 38 | 0 | 0 |
| Add four-arm Xs only | 19 | 40 | 3 | 0 |
| Add both within one budget | 20 | 42 | 3 | 0 |

**Floor photo 00078318.** Three directly supported X detections matched reference crossings near (20.2, 160.3), (68.9, 147.2) and (201.9, 144.4). The fragment-only rule added exactly three Xs on this photo. The other 30 eligible crossings remained unresolved. In the earlier control, the same exit-rule detector recovered 29 of the 33 when given rasterized human lines. [Floor photo](https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test/images/00078318.png) · [Known-line controls](REPORT.md#where-the-crossings-are-lost)

**How additions are chosen.** The four corner crops each cover two thirds of the full width and height, rounded up. Each receives exactly the RGB pixels from the 256-pixel-long-side working image. After tracing we exclude detections within four pixels of an artificial crop edge. A crop detection needs same-type agreement within four pixels from a different crop or the contrast-enhanced full frame. We rank candidates by the number of agreeing views and then by mean separation. Four-arm X proposals use OpenCV's line segment detector on a contrast-enhanced, bilinearly enlarged grayscale image. Collinear fragments must provide at least 5.5 pixels of measured support on *each* arm between 1 and 10 pixels from the intersection. Missing intervals contribute zero support. The combined rule considers the four-arm Xs first, then the crop proposals. It accepts at most one added point per five baseline points in each photo, with more than four pixels between each addition and any retained T or X. [Crop inference](../../../tools/wireframe/predict_tiles.py) · [Measured line support](../../../tools/wireframe/fragment_crossings.py) · [Selection rule](../../../tools/wireframe/evaluate_support.py)

**Second pretrained model.** We also compared HAWPv2's native vectors and rasterized traces at four score thresholds on the disjoint 23-photo development set. We prioritized fewer gap flags, then higher recall for the less well recovered type. This selected native vectors at 0.50. They recovered only 7 Ts and 12 Xs on the issue panel. Tracing the lines at the same cutoff recovered 23 Ts and 34 Xs. Adding HAWP traces to M-LSD when supported by either HAWP vectors or contrast-enhanced M-LSD gave 21 Ts and 37 Xs on the issue panel. On the next 24 it recovered 260/496 Ts and 49/93 Xs, but flags rose from 3/16 to 4/16. HAWP's authors list Wireframe test as a validation set. All draws here come from that same dataset. [HAWP scores](hawp-summary.csv) · [Combination scores](combined-summary.csv) · [Authors' evaluation](https://github.com/cherubicXN/hawp/blob/92ae446875ef9296fc95de7313519aabd638ebb5/docs/HAWPv2.md) · [Upstream splits](https://github.com/cherubicXN/hawp/blob/92ae446875ef9296fc95de7313519aabd638ebb5/configs/hawpv2.yaml)

## Panels and review

The original panel retains the issue's exact labels. The three other panels use the supplementary geometric reference rule documented in the first report. They are disjoint from one another and from the original. The new 24 were the first unused IDs ranked by SHA-256 of `issue3-replication-v1/` plus photo ID. The final 48 used a separate salt, `issue3-confirmation-v1/`, after excluding all prior IDs. We had seen scores for the first three sets when choosing the combined rule. Only the last 48 test that choice prospectively. [Reference rules](REPORT.md#method-and-reproduction) · [24-photo roster and hashes](replication-panel.csv) · [48-photo roster and hashes](confirmation-panel.csv)

The new near-miss screens used RGB crops with the two annotated segments highlighted. Review took place before inspecting local predictions at those sites. Of 48 candidates on the 24 photos, 16 had a visible separation. Of 45 candidates on the 48 photos, 17 had a visible separation. One agent made these judgments. The CSVs include exclusions and their reasons. The counts support a comparison at those positions. Estimating photographic precision would require a separate blinded review of sampled detections. [24-photo reviews](replication-gap-review.csv) · [48-photo reviews](confirmation-gap-review.csv)

## Run on one photo

From this checkout, use the R package dependencies and the Python environment in the original reproduction instructions. `--weights` takes the pinned large FP32 model downloaded by `prepare.py`. The runner accepts a local image without any Wireframe annotations. It writes the type, one-based x and y on the 256-pixel-long-side image, and one-based centers on the EXIF-oriented original. Use `--method crops` to require contour-derived detections throughout. [Single-photo runner](../../../tools/wireframe/photo_junctions.py) · [Environment and model](REPORT.md#method-and-reproduction)

```sh
"$work/venv/bin/python" tools/wireframe/photo_junctions.py photo.jpg junctions.csv \
  --weights "$work/mlsd.tflite"
```

Both single-photo modes reproduced every T and X from the batch run on the floor image. A blank input produced an empty table with column names. Synthetic fragment controls emitted X candidates at a true crossing and none on blank pixels, parallel lines, a T, or a stem stopping six pixels short. The existing R line tests also passed. [Fragment and crop controls](../../../tools/wireframe/test_support.py) · [R line tests](../../../tests/testthat/test-lines.R)

## Reproduce the benchmark

Start with the original report's preparation and comparison commands, which produce the pinned model, manifest and original predictions in `$work`. Run the following from the repository root.

```sh
py="$work/venv/bin/python"
for panel in replication confirmation; do
  "$py" tools/wireframe/prepare.py "$work/$panel" --manifest "$work/issue-3.json" --extra-panel "$panel"
  Rscript --vanilla -e 'pkgload::load_all(".", quiet=TRUE); source("tools/wireframe/compare.R")' "$work/$panel"
done
for folder in "$work" "$work/replication" "$work/confirmation"; do
  "$py" tools/wireframe/predict_tiles.py "$folder"
  Rscript --vanilla -e 'pkgload::load_all(".", quiet=TRUE); source("tools/wireframe/trace_tiles.R")' "$folder"
  "$py" tools/wireframe/fragment_crossings.py "$folder"
done
"$py" tools/wireframe/evaluate_support.py "$work" "$work/issue-3.json" "$work/replication" --confirmation "$work/confirmation"
"$py" tools/wireframe/witnesses.py "$work" --floor-recovery "$work/floor-recovery.png"
for panel in development issue; do
  "$py" tools/wireframe/check_localization.py "$work" "$panel"
done
for panel in replication confirmation; do
  "$py" tools/wireframe/check_localization.py "$work/$panel" "$panel"
done
```

The optional HAWP comparison used the authors' supervised weights, preprocessing, 512 × 512 internal input and inference setting `use_residual=2`. The included CPU patch skips the optional CUDA training extension on a CPU host, initializes a buffer on the selected device and avoids importing the self-supervised package. The predictor verifies the checkpoint SHA-256 before loading it. [Inference code](../../../tools/wireframe/predict_hawp.py) · [CPU patch](../../../tools/wireframe/hawp-cpu.patch)

```sh
uv pip install --python "$py" torch==2.8.0 torchvision==0.23.0 yacs==0.1.8 \
  scikit-image==0.26.0 matplotlib==3.11.2 python-json-logger==4.2.0 tqdm==4.70.1
mkdir -p "$work/hawp-src"
curl -fsSL https://api.github.com/repos/cherubicXN/hawp/tarball/92ae446875ef9296fc95de7313519aabd638ebb5 \
  -o "$work/hawp.tar.gz"
tar -xzf "$work/hawp.tar.gz" --strip-components=1 -C "$work/hawp-src"
patch -s -p1 -d "$work/hawp-src" < tools/wireframe/hawp-cpu.patch
curl -fsSL https://github.com/cherubicXN/hawp-torchhub/releases/download/HAWPv2/hawpv2-edb9b23f.pth \
  -o "$work/hawpv2-edb9b23f.pth"
for folder in "$work" "$work/replication"; do
  "$py" tools/wireframe/predict_hawp.py "$folder" "$work/hawp-src" "$work/hawpv2-edb9b23f.pth"
done
"$py" tools/wireframe/evaluate_hawp.py "$work" "$work/issue-3.json"
"$py" tools/wireframe/combine_hawp.py "$work" "$work/issue-3.json" "$work/replication"
```
