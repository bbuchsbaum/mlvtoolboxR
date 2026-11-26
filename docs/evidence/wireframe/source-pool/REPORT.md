# Two new curve sources reach artist contours the current candidates miss

Structured edges and FDoG at 512 pixels both passed the frozen screen. On the 60 Toronto development photos, an artist-guided selection at each photo's FDoG ink reached **58.55%** of artist contour pixels with FDoG at 512 pixels added to the existing candidates, and **65.28%** with structured edges added, up from **49.94%**. Neither source met the crossing rule. As photo-only drawings, both sources matched the artist's MLV statistics less well than the whole-contour selector, for different reasons. The structured-edge drawing is split into very short curves. The FDoG-512 drawing kept almost no junctions. Both sources remain candidates for a selector that is scored on MLV agreement. [Frozen protocol](PLAN.md) · [Decision record](screening-decision.json)

## Screen results

| Candidate curves | Union, any ink | Artist-guided at FDoG ink | Four-arm crossings of 84 (photos) |
|---|---:|---:|---:|
| FDoG + M-LSD | 58.92% (223,878 px) | 49.94% | 5 (1) |
| + FDoG at 512 px | 80.32% (413,076 px) | **58.55%** | 7 (3) |
| + structured edges | 79.53% (388,389 px) | **65.28%** | 8 (4) |
| + both | 88.19% (535,847 px) | 68.55% | 9 (5) |

Coverage is the fraction of the 353,073 artist pixels within two working pixels of a selected curve. The artist-guided selection repeatedly takes the complete curve that adds the most artist pixels per added ink pixel, never exceeding the photo's FDoG ink. It shows what the curves make available. It is not a photo-only result. [Pooled and by category](source-pool-summary.csv) · [Per photo](source-pool-photo-scores.csv)

The rule carried a source forward for a gain of at least 2.0 points in the artist-guided selection, or at least five more four-arm crossings on at least three photos. FDoG at 512 pixels gained **+8.61 points** and structured edges **+15.34**. The artist-guided gain was positive on all 60 photos for each source. Crossings rose by two and three, below the crossing threshold. All the gained crossings are on city and office photos.

## How much of the gain is new territory

These checks were added after the screen and lie outside its rule. The structured-edge curves are short: their median curve inks two working pixels, against 13 for FDoG at 512 pixels, 25 for FDoG and 7 for M-LSD. Shorter curves let a selection spend ink more precisely. Two checks separate that effect from new contours:

- **Pixels the old pool cannot reach.** Some artist pixels lie farther than two pixels from every FDoG and M-LSD curve. These pixels contribute 14.75 points of the FDoG-512 selection and 17.92 points of the structured-edge selection.
- **Granularity control.** Each existing curve's raster pixels were split into consecutive runs of at most 2 or 13 pixels. The runs follow the rasterizer's sorted column-major order, not arc length, so a run can join pixels that are not adjacent along the curve. The existing candidates' artist-guided coverage rose to 58.74% and 53.20%. Structured edges remain 6.54 points above the 2-pixel control, and FDoG at 512 pixels 5.35 points above the 13-pixel control. [Granularity control](granularity-check.csv)

## Photo-only drawings from each source

Each source also made a drawing on its own without artist information. It took complete curves in order of decreasing mean edge strength (for FDoG at 512 pixels, decreasing length) until the next would exceed the FDoG ink. An independent recomputation from the rendered drawings reproduced each coverage and precision figure.

The MLV columns use the measures of the [agreement analysis](../mlv-agreement/REPORT.md), restricted to the 60 development photos. The artist-trained classifier is fitted by leave-one-out on the other 59 development artist drawings, and identification ranks against these 60 artist drawings, so the values differ from the 180-photo analysis. This comparison came after the screen and lies outside its rule. [MLV agreement of the source drawings](m3-mlv-agreement.csv)

| Development drawing | Coverage | Precision | F1 | Artist-trained category decoding | Mean MLV histogram distance | Contours / artist | Junctions in 60 photos |
|---|---:|---:|---:|---:|---:|---:|---:|
| FDoG | 42.65% | 75.75% | 0.546 | 0.250 | 0.344 | 0.25 | 588 |
| FDoG at 512 px | 42.88% | 80.31% | 0.559 | 0.233 | 0.492 | 0.21 | 4 |
| Structured edges | **46.68%** | **89.56%** | **0.614** | 0.283 | 0.379 | 1.23 | 1,021 |
| Whole-contour selector | 44.91% | — | 0.578 | **0.383** | **0.300** | 0.34 | 1,055 |
| Artist | — | — | — | — | — | 1 | 9,936 |

The selector's coverage and F1 are its out-of-fold [development scores](../crossing-groups/REPORT.md). Junctions are counted in the original frame with the package detector. Structured edges gave the best pixel agreement of the drawings compared here, yet their MLV statistics agree less well than the selector's. With a third of the artist's contour length, the drawing has 1.23 times as many contours as the artist's. Its length histogram is the farthest from the artist's of these methods (0.469 against 0.320 for the selector), and its T density does not track the artist's across photos (ρ −0.09). Pixel coverage alone rewards this fragmentation.

The FDoG-512 drawing fails for another reason. Its photo-only subset contains four junctions across the 60 photos, and 57 photos have none. Its junction-type and junction-angle distances are therefore near the maximum (0.968 and 0.991). In three spot-checked photos, the complete FDoG-512 drawings had junctions that the length-ordered subset lost. The ordering rule appears to drop the short branches that form them.

Linking the structured-edge fragments into longer contours is a hypothesis to test, not a result. The [reference toolbox](../../../../vendor/MLV_toolbox/MLVcode/traceLineDrawingFromRGB.m) groups structured edges with superpixel affinities before tracing. This screen used the non-maximum-suppressed edges directly, so it does not reproduce that route.

## Registration

This check was also added after the screen. Shifting each full source drawing by up to two pixels against the artist raster gave the best F1 at zero shift for structured edges, traced at 800 pixels. FDoG and M-LSD scored best at (−1, −1), and FDoG at 512 pixels at (0, −1). The largest gain from shifting was 0.004 in F1 (FDoG: 0.550 against 0.546 unshifted). The FDoG-family and M-LSD drawings may share a small coordinate offset whatever their tracing resolution. Its source has not been found. [Shift check](alignment-check.csv)

## Reproduction

Use the private Toronto inputs of the [whole-contour reproduction](../curve-rank/REPRODUCE.md). The structured-edge step needs a separate environment with `opencv-contrib-python-headless==4.14.0.94`, numpy 2.2.6 and Pillow 12.1.1. It uses the `opencv_extra` model at commit `5e3a56880fb115e757855f8e01e744c154791144` (SHA-256 `fe3586bc…bf17`). Non-maximum suppression uses OpenCV's default radii. The FDoG and M-LSD baselines, the earlier ledger and the earlier arm table were reproduced photo by photo before any source was added.

```sh
sed=.agent-work/venv-sed/bin/python
py="$work/venv/bin/python"
pool=.agent-work/source-pool
draw=.agent-work/contour-study/draw-toronto
"$sed" -I tools/wireframe/source_pool/prepare_photos.py .agent-work/toronto/CP_images.zip \
  docs/evidence/wireframe/curve-rank/development.csv "$pool/model/model.yml.gz" "$pool"
Rscript --vanilla tools/wireframe/source_pool/trace_sources.R .agent-work/toronto/staged "$draw" \
  "$pool" docs/evidence/wireframe/curve-rank/development.csv
"$py" tools/wireframe/source_pool/traced_arm_pool.py .agent-work/toronto/staged "$draw" \
  .agent-work/crossing-groups/baseline "$pool/draw" "$pool/arms"
"$py" tools/wireframe/source_pool/screen.py "$draw" .agent-work/contour-proposal "$pool" \
  "$pool/arms/traced-pool-development.csv" docs/evidence/wireframe/source-pool
"$py" tools/wireframe/source_pool/granularity_check.py "$draw" "$pool" docs/evidence/wireframe/source-pool/granularity-check.csv
"$py" tools/wireframe/source_pool/alignment_check.py "$draw" "$pool" docs/evidence/wireframe/source-pool/alignment-check.csv
"$py" -I tools/wireframe/source_pool/m3_selection.py "$pool" \
  docs/evidence/wireframe/source-pool/source-pool-photo-scores.csv "$pool/m3-selection.csv"
Rscript --vanilla tools/wireframe/source_pool/m3_statistics.R "$pool" "$pool/m3-selection.csv" "$pool/m3-stats"
"$py" -I tools/wireframe/source_pool/m3_agreement.py .agent-work/mlv-agreement "$pool/m3-stats" "$draw" \
  docs/evidence/wireframe/source-pool/m3-mlv-agreement.csv
```

The last three commands need the agreement analysis's per-drawing statistics in `.agent-work/mlv-agreement`.
