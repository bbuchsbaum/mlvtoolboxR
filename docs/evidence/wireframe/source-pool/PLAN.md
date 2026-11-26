# Screening additional candidate curve sources

On the 60 Toronto development photos, the union of all FDoG and M-LSD curves reaches 58.92% of artist contour pixels within two working pixels, even with no ink limit. An artist-guided greedy selection at the per-photo FDoG ink count reaches 49.94%, and the photo-only whole-contour selector 44.91%. This development screen asks whether additional photo-derived curve sources contain artist contours that the current candidates lack. It uses artist drawings to measure what is available and cannot establish a photo-only gain. The 120 test photos and the 205 unscored reserve photos are not used. [Earlier headroom ledger, local](../../../../.agent-work/contour-proposal/headroom.csv) · [Development roster](../curve-rank/development.csv)

## Sources

- **FDoG at 512 pixels.** Apply the package's FDoG drawing (`xdog_lines()` defaults, as in the existing `fdog` source) to the photo resized with Pillow LANCZOS to 512 pixels on its long side. Trace with `trace_contours(simplify = 1)`, then map the coordinates to the working frame with the staging transform.
- **Structured edges.** Use the structured edge detector of Dollár and Zitnick, the original photo route of the reference MLV toolbox, through OpenCV's `ximgproc` module. Use the pretrained BSDS model from `opencv_extra`, pinned by commit and SHA-256. Run it on the original 800 × 600 photo, apply the module's non-maximum suppression, keep pixels with edge strength at least 0.1, trace them with `trace_contours(simplify = 1)`, and map the coordinates to the working frame. A source curve's strength is the mean edge strength along its pixels.

Each source is rasterized curve by curve with the package rasterizer at the working resolution, in the same pixel-index form as the existing per-curve exports.

## Measurements

- **M1.** Coverage and distinct ink of the union of candidate curves with no ink limit, with and without each new source.
- **M2.** Artist-guided greedy coverage at each photo's FDoG ink count, using the rule in the earlier headroom diagnostic (largest additional artist pixels per additional ink pixel, complete curves only, never above the budget), with and without each new source.
- **M3.** A photo-only drawing from the source alone. Order its curves by decreasing strength (for FDoG at 512 pixels, by decreasing length) and take complete curves until the next would exceed the photo's FDoG ink. Report coverage, precision and F1 against the artist ink, alongside FDoG's 42.65%, 75.75% and 0.5457.
- **M4.** The number of the 84 measurable development crossings with four supported arms when the source joins the FDoG and M-LSD curves, using the traced-branch rule of the [arm opportunity check](../arm-opportunity/REPORT.md). The baseline is five, on one photo.

Report the totals pooled over the 60 photos and by category.

## Screening rule

Carry a source forward to refitting the selector if adding it raises M2 by at least 2.0 percentage points, or raises M4 by at least five crossings spread over at least three photos. Report sources that fail and drop them. A carried source enters a refitted selector under a separately frozen protocol before any confirmation photos are opened.
