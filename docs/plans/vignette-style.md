# Vignette plan and style sheet (2026-11-20)

## Readers and job

Vision scientists and their students who use R but are not computational
geometers. They want to measure line drawings, understand what each number
means well enough to defend it in a methods section, and choose between the
corrected defaults and published-compatible MLV values. Voice (PRODUCT.md):
precise, technical, peer-to-peer; units always shown; no hype. Follow the
`write-with-judgment` technical guidance: begin with the thing itself, name real
functions and values, explain cause and consequence, define terms at first use,
put an example where uncertainty first appears, separate guarantees from
defaults and recommendations, complete sentences, no "simply/just/powerful/
seamless/robust" as praise.

## Vignette set

| File | Title | Job |
|---|---|---|
| `getting-started.Rmd` | From a line drawing to contour features | Build a drawing, measure orientation, length and curvature, read histograms and feature tables, use a feature specification. |
| `junctions.Rmd` | Junctions: how they are found and what they measure | Junction types, the gap allowance, the arm rule, merging; where the corrected rule and `rule = "mlv2023"` differ, shown on small drawings. |
| `medial-axis.Rmd` | Medial axis and symmetry scores | Distance map → average outward flux → skeleton → branches → parallelism, separation, mirror, taper → mapped contour scores. |
| `mlv-compatibility.Rmd` | Reproducing published MLV values | The two modes, the exact compatible recipe, what each correction changes and why, cache/spec rejection, and the qualification evidence. |
| `workflows.Rmd` | (revise) File input, extraction and feature recipes | Keep the I/O, extraction and checked-prediction material; remove anything now covered elsewhere and correct stale claims (it still says the Pratt backend is used and MATLAB morphology equivalence is unestablished). |
| `stimulus-design.Rmd` | unchanged | Cross-link only. |

Cross-link with `vignette("junctions", package = "mlvtoolbox")` style references.
Each new vignette roughly 1,500–2,500 words of prose; whole set renders in well
under a minute.

## Shared vocabulary (use exactly these terms)

drawing (`mlv_drawing`), contour, segment, vertex; orientation (degrees, axial
modulo 180 for histograms), length (pixels), curvature (degrees per pixel);
junction, arm, junction neighbourhood, gap allowance; distance map, average
outward flux (AOF, "flux"), skeleton, branch, branch score; parallelism,
separation, mirror, taper; mapped (contour) score; feature specification;
"corrected defaults" vs "published-compatible mode" (`mlv2023` options).
Image coordinates: x rightward, y downward, pixel centres at integers.

## Figures

- Base graphics only (no new hard dependencies). `fig.retina = 2`,
  `dev = "png"`, sizes chosen per figure (typically 7 × 3.5 in for paired panels).
- Stimuli on a white ground with black ink; colour only when it encodes a
  quantity, always with a legend that names the quantity and unit.
- Palettes: orientation is cyclic → a cyclic palette (e.g. `hcl.colors(n,
  "Hue")`-style or a custom cyclic HCL wheel) over 0–180°; magnitudes →
  sequential (`"Viridis"` or `"Blues 3"`); signed flux → diverging
  (`"Blue-Red 3"`, symmetric limits); categories → Okabe–Ito.
- Equal aspect ratio, y axis downward to match image coordinates, minimal
  axes (light, or none on stimulus panels), panel titles that say what is shown.
- Annotate the point of each figure (an arrow, a label, a highlighted junction)
  rather than relying on the caption. Every figure gets a one-sentence
  `fig.cap` stating what to notice.
- Restore graphics state (`op <- par(...); on.exit(par(op))` inside helpers or
  restore at chunk end).
- Prefer public plotting methods. Cosmetic helper code goes in a hidden chunk
  (`echo = FALSE`) only when the lesson is not plotting.

## Computation rules

- Every number stated in prose comes from an evaluated chunk (inline `r`) or is
  explicitly cited to a committed evidence document (`docs/PARITY.md`) as a
  qualification result not recomputed in the vignette. Natural reference data
  are not available to vignettes; use synthetic drawings only.
- Add quiet `stopifnot()` checks where a degenerate output could mislead.
- Inspect every rendered figure (open the PNGs) before declaring done.
