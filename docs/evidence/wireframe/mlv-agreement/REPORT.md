# Do automatic drawings reproduce the artist's MLV statistics?

Of the four automatic methods, the whole-contour selector came closest to the artist drawings, though well short of them. A classifier trained on artist-drawing statistics assigned the six Toronto scene categories correctly to **84%** of held-out artist drawings. Applied to the selector's drawings of the same photographs it was correct on **42%**, compared with **28%** for FDoG. Orientation statistics agreed closely, photograph by photograph. Length and junction-angle histograms identified the photograph only weakly. Curvature and junction-type histograms did so no better than chance. In a post hoc control, rasterizing and retracing the artist's own drawing at the 256-pixel working resolution produced disagreement of similar size. The direction was opposite, though: too many contours and junctions rather than too few. [Frozen protocol](PLAN.md)

## Category information

The analysis uses the 180 Toronto photographs that have drawings from every method: the 60 development photos with out-of-fold selector drawings and the 120 test photos with frozen-selector drawings. The 205 unscored reserve photos remain unopened. Each drawing is described by 38 MLV features: orientation, length, curvature and junction-angle histograms, the density of each junction type, and total contour length. One detector and one frame are used for artist and automatic drawings alike. Leave-one-photo-out multinomial logistic regression has a chance accuracy of 1/6. [Decoding results](category-decoding.csv)

| Drawings | Artist-trained classifier, all 180 | Same, 120 test photos | Classifier trained on the method's own drawings |
|---|---:|---:|---:|
| Artist | **0.839** (0.789–0.889) | 0.850 | — |
| FDoG | 0.278 (0.233–0.328) | 0.275 | 0.483 |
| M-LSD | 0.367 (0.311–0.422) | 0.358 | 0.611 |
| Earlier vector blend | 0.350 (0.294–0.411) | 0.350 | 0.689 |
| Whole-contour selector | **0.422** (0.372–0.478) | 0.408 | 0.550 |

The selector improved artist-scale decoding over FDoG by **+0.144** (paired 95% interval +0.083 to +0.206), and by +0.133 (+0.058 to +0.217) on the test photos alone. Its gain over the earlier blend was +0.072 (+0.006 to +0.133) pooled, and +0.058 (−0.025 to +0.142) on the test photos alone. Every automatic method carried more category information when the classifier was trained on its own drawings. Their statistics separate the categories on a scale different from the artist's. Intervals resample photographs within category 2,000 times. [Paired intervals](paired-bootstrap.csv)

## Which statistics agree

The table gives the median Jensen–Shannon distance (0 for identical histograms, 1 for disjoint ones) between each automatic drawing and the artist drawing of the same photograph. Beside it is the mean percentile at which that artist drawing ranks among all 180 artist drawings, where 0.5 is chance. The reference column gives the median distance between artist drawings of *different* photographs in the same category. The CSV also holds mean distances and the top-1 identification rate. The top-1 rate was at most 0.200 for orientation and horizontal/vertical balance, and at most 0.028 for the other properties. Two pairs were excluded because both histograms were empty: one FDoG drawing and its artist drawing for junction type and for junction angle. [Distribution agreement](distribution-agreement.csv) · [Per photo](photo-distances.csv)

| Property | Different photos, same category | FDoG | M-LSD | Earlier blend | Selector |
|---|---:|---:|---:|---:|---:|
| Orientation | 0.226 | 0.174 / 0.90 | 0.310 / 0.82 | 0.235 / 0.90 | **0.145 / 0.95** |
| Horizontal/vertical | 0.207 | 0.165 / 0.90 | 0.314 / 0.80 | 0.224 / 0.89 | **0.146 / 0.93** |
| Length | 0.259 | 0.336 / 0.57 | 0.357 / 0.61 | 0.372 / 0.61 | 0.306 / 0.58 |
| Curvature | 0.149 | **0.252** / 0.54 | 0.463 / 0.50 | 0.367 / 0.50 | 0.295 / 0.51 |
| Junction type | 0.192 | 0.468 / 0.52 | **0.305** / 0.51 | 0.401 / 0.50 | 0.386 / 0.53 |
| Junction angle | 0.263 | 0.565 / 0.54 | 0.440 / 0.59 | 0.600 / 0.58 | **0.414** / 0.58 |

For orientation and horizontal/vertical balance, FDoG and the selector each lie closer to the photograph's artist drawing than two artist drawings from the same category lie to each other. The blend lies slightly farther. The photograph's own artist drawing ranks, on average, in the top tenth of the 180. For length, curvature and junctions, every automatic drawing is farther from its own artist drawing than different artist drawings are from one another. Length and junction-angle histograms still rank the right photograph above chance, though weakly (mean percentiles 0.54–0.61; for the selector, 0.580 and 0.584). Curvature and junction-type histograms do not (0.50–0.54). Comparing means with paired intervals, the selector was closer to the artist than FDoG on every property except curvature, which worsened by 0.034 (0.027 to 0.041).

## Scalar summaries

| Summary, original frame | Artist median | FDoG | Selector |
|---|---:|---:|---:|
| Total contour length (px) | 19,856 | ρ 0.74, ratio 0.36 | ρ 0.73, ratio 0.36 |
| Contours | 260.5 | 0.61, 0.25 | 0.74, 0.34 |
| T junctions per 10,000 px | 50.2 | −0.12, 0.09 | **0.48, 0.20** |
| Y junctions per 10,000 px | 5.4 | 0.10, 0.72 | 0.23, 0.72 |
| X junctions per 10,000 px | 5.9 | 0.07, 0.00 | 0.26, 0.78 |
| Length-weighted mean curvature | 2.18 | 0.36, 0.39 | 0.22, 0.23 |
| Fraction near horizontal or vertical | 0.77 | 0.81, 0.98 | 0.80, 1.03 |

ρ is the Spearman correlation with the artist value across the 180 photos. The ratio is the median per-photo ratio of method to artist. The automatic drawings contain about a third of the artist's contour length. The selector's T density is a fifth of the artist's, but it now ranks photographs by T density (ρ 0.48, against −0.12 for FDoG; paired difference +0.60, 0.43 to 0.78). Mean curvature is 0.23 of the artist's for the selector (about four times smoother) and 0.39 for FDoG. The blend and M-LSD, dominated by straight segments, reach 0.09 and 0.03. Orientation balance matches. [Scalar agreement](scalar-agreement.csv) · [By category](by-category.csv)

Both sensitivity analyses gave the same pattern. They were the 256-pixel working frame, and the published MLV curvature and junction definitions in the original frame. Under the published definitions, artist-trained decoding was 0.350 for FDoG and 0.439 for the selector.

## How much comes from the drawing representation

This control was added after the primary results were read and is not part of the frozen plan. Each artist drawing was rasterized with 256, 512 or 800 pixels on its long side, retraced with `trace_contours(simplify = 1)`, and scored exactly like a method drawing. The content is the artist's own; only the representation changes. Every automatic source in this study was traced at 256 pixels. [Control tables](control/)

| Drawing | Artist-trained decoding, original frame | Same, working frame | Mean histogram distance | Identification percentile | Contours / artist | T density / artist | Y density / artist | Mean curvature / artist |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Artist, retraced at 256 px | 0.267 | 0.300 | 0.282 | 0.68 | 2.93 | 2.91 | 14.6 | 0.29 |
| Artist, retraced at 512 px | 0.328 | 0.539 | 0.207 | 0.74 | 2.64 | 2.28 | 7.86 | 0.49 |
| Artist, retraced at 800 px | **0.544** | **0.672** | **0.165** | **0.80** | 2.27 | 1.46 | 3.72 | 0.61 |
| Whole-contour selector | 0.422 | 0.456 | 0.303 | 0.68 | 0.34 | 0.20 | 0.72 | 0.23 |

Distances, percentiles and ratios are in the original frame; the ratios are medians over photos. Retraced at the working resolution, the artist's own contours disagreed with the artist drawing about as much as the selector's drawings did. They gave lower artist-scale decoding in both frames. The scalar signature differs: the retraced drawings have about three times as many contours and T junctions as the artist's and fourteen times as many Ys, whereas the automatic drawings have too few.

Higher tracing resolution reduced every distance. At 800 pixels, decoding rose to 0.544 (0.483–0.611) in the original frame and 0.672 in the working frame, and the curvature histogram again identified its photograph (percentile 0.71). This remains short of the 0.839 for the artist's vector drawings. The representation can therefore limit MLV agreement whatever contours a method finds. The control does not separate the steps of that representation:

- rasterization;
- thinning;
- splitting at raster junctions;
- simplification, with a tolerance of one raster pixel that shrinks in the original frame as resolution rises.

Which of them causes the excess contours and Y junctions remains a hypothesis. Decoding also depends on the measurement frame, because length bins and junction tolerances are in pixels. The published MLV definitions give the same ordering in the original frame (0.283, 0.328 and 0.428).

## Reproduction

The commands read the private Toronto inputs listed in the [whole-contour reproduction](../curve-rank/REPRODUCE.md). They write the per-drawing histograms to a private directory and the aggregate tables here. The recorded run used R 4.5.1 and the Python 3.12.11 environment of the earlier comparisons.

```sh
py="$work/venv/bin/python"
Rscript --vanilla tools/wireframe/mlv_agreement/compute.R .agent-work/mlv-agreement
"$py" -I tools/wireframe/mlv_agreement/analyze.py .agent-work/mlv-agreement docs/evidence/wireframe/mlv-agreement
for res in 256 512 800; do
  Rscript --vanilla tools/wireframe/mlv_agreement/control.R .agent-work/mlv-agreement "$res"
done
"$py" -I tools/wireframe/mlv_agreement/analyze.py .agent-work/mlv-agreement \
  docs/evidence/wireframe/mlv-agreement/control --control
```

Two steps deviate from the written protocol. The features use log(1 + total contour length) instead of the log, because some M-LSD drawings are empty. Test-only decoding restricts the 180 leave-one-out predictions to the test photos, so its training rows include development photos.
