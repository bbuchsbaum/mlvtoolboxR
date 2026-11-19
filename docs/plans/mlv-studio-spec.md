# MLV Studio — specification (draft)

Status: proposal, 2026-10-26. Visual design approved 2026-10-26 (see §8);
implementation not started.

## 1. Goal

A real, polished desktop-grade interface for mlvtoolbox that ships inside the
package and launches from R:

```r
library(mlvtoolbox)
s <- mlv_studio(drawing)
```

It is not Shiny. It is a prebuilt single-page app served from the package's
`inst/` directory by a minimal local server. R remains the only compute
engine; the browser renders, inspects, and drives R. The live R session and
the UI share objects in both directions.

Design principles:

1. **R is the source of truth for every number.** No measurement is
   reimplemented in JavaScript. This preserves the parity guarantees in
   `docs/PARITY.md` and avoids two drifting implementations.
2. **Every UI action is an R call.** The UI keeps a provenance log that can
   be exported as a runnable R script that reproduces the session state.
3. **Geometry first.** The drawing canvas is the centre of the product; tables
   and charts are linked views of it, not separate pages.
4. **Zero cost to non-UI users.** UI dependencies live in `Suggests`; the core
   package does not change.

## 2. What it contains

### 2.1 Workspace layout

```
┌────────┬──────────────────────────────────────────┬──────────────┐
│Library │  Canvas (drawing / image / MAT layers)   │  Inspector   │
│        │                                          │              │
│drawings│  pan · zoom · hover · select · lasso     │ selection    │
│images  │                                          │ properties   │
│MATs    │                                          │ parameters   │
│splits  ├──────────────────────────────────────────┤              │
│        │  Drawer: histograms · feature table ·    │              │
│        │  provenance / R code                     │              │
└────────┴──────────────────────────────────────────┴──────────────┘
```

Command palette (⌘K) for every action; keyboard shortcuts for layer toggles,
property colouring, and selection tools. Light and dark themes, with the
drawing surface itself always rendered as "paper" (white ground) unless the
user switches it, because stimuli are judged on white. *Superseded by the
approved design (§8): the chrome is a calibrated mid-grey field.*

### 2.2 Modules

**A. Library / import**
- Drag-and-drop PNG/JPEG, SVG, MLV `.mat`, and `.rds` drawings; folder import
  for batches.
- Objects pushed from R appear here live (`s$push()`).
- Thumbnails, canvas size, contour/segment counts, and derived-state badges
  (junctions computed, MAT computed, stale).

**B. Tracing workbench** (`read_mlv_image`, `trace_contours`, `trace_rgb`,
`labels_to_edges`, `simplify_contours`)
- Source image underlay with the traced vector overlay.
- Parameter controls with debounced re-trace and a before/after split slider.
- Diagnostics: dropped zero-length segments, degenerate geometry, contour
  count deltas as parameters change.

**C. Drawing canvas** (core view)
- Vector rendering of contours with colour-by: orientation (cyclic map),
  length, curvature (sequential), contour ID (categorical), junction
  membership, or any per-contour model score (`predict_contours`).
- Hover tooltip: contour ID, segment index, length, orientation, curvature.
- Junction glyphs by type (L, T, X, Y, arrow), with angle annotation on hover
  (`detect_junctions`, `segment_at_junctions`, `split_junction_regions`).
- Selection: click, shift-click, lasso, "select by predicate" (e.g.
  curvature > x), producing a contour-ID set shared with R
  (`s$selection()`).
- Transform tools as R-backed operations: rotate, shift, clip to circle,
  subset (`rotate_drawing`, `shift_contours`, `clip_circle`,
  `subset_contours`). Each produces a new library object with lineage.

**D. Medial-axis view** (`compute_mat`, `average_outward_flux`,
`skeletonize`, `trace_skeleton`, `compute_mat_properties`,
`map_mat_properties`, `analyse_mat`)
- Layer stack: distance map, AOF map, skeleton, traced branches, score images
  for parallelism / separation / mirror / taper.
- Branch inspector: per-point score profile along the branch with radius
  plot; clicking a point highlights the tangent-mapped contour points.
- Parameter panel (thresholds, `window`, `power`, minimum branch points) with
  explicit "recompute" (MAT is the expensive path; no auto-recompute).
- Mapping coverage display: contour pixels with no MAT assignment shown
  distinctly from low scores.

**E. Statistics** (`contour_statistics`, `feature_table`, `feature_density`,
`average_property`, `batch_statistics`, `plot_properties`)
- Histograms per property with bin edges and units visible, not just
  `length_4`. Linked brushing: brushing a bin selects its contours on the
  canvas and vice versa.
- Feature table grid (sortable, column filters), per drawing or per contour.
- Compare mode: two to four drawings side by side with overlaid histograms
  and difference bars.
- Dataset view for batches: small multiples and a drawings × features matrix.
- When `feature_spec` lands (roadmap item 3), the active spec is shown in the
  header and incompatible comparisons are blocked in the UI just as in R.

**F. Stimulus lab** (`split_contours`, `split_by_histogram`, later the
constrained designer)
- Choose criterion/criteria, weights, and length fraction; see top / bottom /
  unassigned drawings side by side, updating live for cheap criteria.
- Matching report: for each nuisance property, the residual difference
  between conditions (histogram distance, total length, coverage) with
  pass/fail against user tolerances.
- Set-level mode across a batch: condition separation versus nuisance
  residuals for the whole stimulus set.
- Export: `write_mlv_png` at exact pixel dimensions, SVG, CSV of features,
  plus a manifest (seed, parameters, package version, contour IDs retained).
- The constrained designer (roadmap item 2) slots into this module
  unchanged: its feasibility, achieved separation, and residuals map onto the
  same report.

**G. Provenance and R code**
- Every operation recorded as `{op, args, inputs, output, pkg_version}`.
- "Copy as R" renders the log as a script using exported functions only;
  CI tests that the generated script reproduces the same objects.
- Session save/load as JSON plus RDS objects.

**H. Estimator lab** (after roadmap item 1)
- Overlay curvature estimators (published endpoint convention vs turning
  angle vs arc-length scale) on the same contour; invariance demo
  (subdivide / reverse / rotate) with live numbers. Useful for teaching and
  for the methods paper.

## 3. R interface

```r
mlv_studio(x = NULL, ..., browser = getOption("mlvtoolbox.browser", TRUE),
           port = NULL, block = !interactive())
```

- `x`: a drawing, a named list of drawings/MATs, an image path, or a folder.
- `browser`: `TRUE` (system browser), `"viewer"` (RStudio/Positron viewer
  pane when available), or `FALSE` (print URL only).
- `block = FALSE` in interactive sessions: the server runs on the `later`
  event loop while the console is idle, so the user keeps working in R.
  `block = TRUE` for `Rscript`.

Returns an `mlv_studio` handle (environment-backed, no R6 dependency):

| Method | Purpose |
|---|---|
| `s$push(obj, name)` | Send an object to the library |
| `s$get(name)` | Pull a UI-created object (split, transform result) into R |
| `s$list()` | Library contents with lineage |
| `s$selection()` | Currently selected contour IDs on the active drawing |
| `s$code()` | Provenance log rendered as R code |
| `s$url()`, `s$stop()` | URL, shutdown |

`print(s)` shows URL, object count, and status.

Later: `mlv_snapshot(x, file)` writes a self-contained read-only HTML viewer
with embedded data, for sharing with collaborators or embedding in
vignettes/pkgdown.

## 4. Architecture

### 4.1 Components

```
R process                                   Browser
─────────                                   ───────
mlv_studio()  ──► httpuv server (127.0.0.1)  ◄── static SPA (inst/studio/)
                  ├─ GET  /           static assets
                  └─ WS   /rpc        JSON-RPC + binary frames
session registry (env): id → object, lineage
RPC handlers: whitelisted wrappers over exported functions
```

- **Server:** `httpuv` (static files + WebSocket). `jsonlite` for messages.
  Both in `Suggests`; `mlv_studio()` checks with `requireNamespace()` and
  gives an install hint.
- **RPC:** a fixed table of handler functions (`.studio_handlers`), each a
  thin wrapper that validates arguments, calls an exported function, stores
  the result in the registry, appends to provenance, and returns an object
  ID plus a summary. **There is no eval endpoint.**
- **Payloads:** drawings go as binary frames: `Float32Array` of segment
  endpoints (n × 4), `Int32Array` contour offsets, plus a JSON header (IDs,
  width, height, units). Raster maps (distance, AOF, scores) go as
  `Float32Array` with dims; the client applies colour maps on the GPU. This
  keeps a 50k-segment drawing to about 1 MB.
- **Protocol versioning:** a `protocol` integer in the handshake; client and
  server refuse a mismatch with a clear message (prevents stale cached
  bundles after a package upgrade).
- **Long operations:** handlers run on the R main thread. Rcpp kernels
  already call `checkUserInterrupt`; add an optional progress callback to
  the expensive stages (`compute_mat`, tracing, batch) that sends progress
  frames. A cancel button sends an interrupt request. Out-of-process workers
  (mirai/callr) are explicitly deferred: they complicate object sharing and
  are not needed until batch sizes demand it.

### 4.2 Front end

Stack aligned with psycloud so tooling and habits transfer:

- TypeScript, React 19, Vite, Tailwind 4.
- Canvas rendering: a custom WebGL2 layer (instanced line segments with
  per-segment attributes; textures for raster maps) with a Canvas2D fallback.
  Hit-testing through a static spatial index (`flatbush`) over segments, not
  GPU picking, so hover stays exact.
- Charts: Observable Plot for histograms and small multiples; a virtualised
  grid for the feature table.
- State: a small store (Zustand) mirroring the R registry by ID; the client
  never mutates drawings, it requests operations.
- Colour maps: perceptually uniform, with a cyclic map for orientation
  (orientation is mod 180°), sequential for magnitudes, colour-blind-safe
  categorical palette for junction types. Every colour encoding has a legend
  with units.

### 4.3 Repository layout and packaging

```
studio/                 # front-end source (in .Rbuildignore)
  src/ package.json vite.config.ts
inst/studio/            # built, minified bundle committed to git
R/studio.R              # mlv_studio(), handle, server
R/studio-handlers.R     # RPC table
R/studio-protocol.R     # binary encoders, protocol version
```

- Bundle budget: **≤ 1.5 MB** gzipped-on-disk total in `inst/studio/`, to
  stay well under CRAN's installed-size NOTE threshold. Fonts subset or
  system stack; no icon fonts (inline SVG icons).
- CRAN JS policy: include a pointer to the unminified source (`studio/` on
  GitHub, exact commit recorded in `inst/studio/BUILD.json`) and list bundled
  third-party licences in `inst/COPYRIGHTS`. All bundled JS must be
  GPL-3-compatible (MIT/ISC/BSD fine).
- A CI job rebuilds `studio/` and fails if the result differs from the
  committed `inst/studio/` (no stale bundles).

### 4.4 Security

A local server that drives R is a code-execution surface. Requirements:

- Bind to `127.0.0.1` only; random port by default.
- Random per-session token in the launch URL; required on the WebSocket
  handshake. Check the `Origin` header.
- Whitelisted handlers only; argument validation in R before any call.
- File access only to paths the user supplied from R or bytes uploaded
  through the browser (written to `tempdir()`); exports go to a directory the
  user picks in R or via browser download.
- Server stops on `s$stop()`, on R session end, and after an idle timeout
  when `block = TRUE`.

## 5. Testing and evidence

- **R handlers:** testthat, calling handlers directly without a browser
  (argument validation, registry, provenance, binary encoding round-trip).
- **Provenance round-trip:** generated R script reproduces identical objects
  (`identical()` on drawings, tolerance on doubles) for a scripted session.
- **Protocol:** golden binary payloads shared between R tests and vitest.
- **Front end:** vitest for state and encoders.
- **End to end:** Playwright against a real `mlv_studio()` server using the
  project-installed Playwright browser; screenshot regression for canvas
  rendering of fixture drawings. Skipped on CRAN, required in CI.
- **Performance:** a 50k-segment drawing pans and zooms at 60 fps; hover
  latency under 16 ms; push of that drawing under 200 ms locally.
- **Accessibility:** keyboard reachability for every action, WCAG AA
  contrast for chrome (not for stimulus rendering).

## 6. Milestones

| Milestone | Contents | Exit criterion |
|---|---|---|
| M0 — skeleton | server, token auth, protocol, `mlv_studio()` handle, push/get, canvas with colour-by and hover | Playwright test opens a fixture drawing and reads a hover value that equals R's |
| M1 — inspect | library, junction glyphs, selection ↔ R, statistics with linked brushing, provenance → R code | generated script reproduces the session |
| M2 — trace | tracing workbench with live re-trace and diagnostics | image → drawing entirely in UI, reproduced by script |
| M3 — stimulus lab | splits, matching report, transforms, export with manifest | exported PNG byte-identical to `write_mlv_png` from the script |
| M4 — medial axis | MAT layers, branch inspector, tangent mapping, progress/cancel | MAT session on fixture drawing reproduced by script |
| M5 — share | `mlv_snapshot()` static HTML; dataset view | snapshot opens offline with no server |
| Later | estimator lab, constrained designer, `feature_spec` integration, optional WebR demo on pkgdown | tied to roadmap items 1–3 |

M0–M1 is the first useful release: a beautiful inspector with R round-trip.

## 7. Open decisions

1. **Name:** `mlv_studio()` proposed.
2. **Front-end framework:** React (matches psycloud) vs Svelte/Solid
   (smaller bundle). Recommendation: React, for shared tooling; the bundle
   budget is achievable either way.
3. **`httpuv`/`jsonlite` in `Suggests` vs `Imports`:** recommendation
   `Suggests`, keeping the core package dependency-light.
4. **WebR:** r-universe builds WebAssembly binaries, so a browser-only demo
   (no local R) is possible later; it would need an in-browser transport
   alongside the WebSocket one. Deferred; design the client transport as an
   interface so it can be added.
5. **Ordering against the measurement roadmap:** the UI should not wait for
   roadmap items 1–3, but modules E/F/H should be built against the
   `feature_spec` and estimator-naming interfaces once they are fixed, to
   avoid relabelling bins in the UI later.

## 8. Approved design (2026-10-26)

The visual design is approved. The reference is the claude.ai Design canvas
"MLV Studio" (private; https://claude.ai/artifact/JZZcNovsqhvRksphdbNSUj,
version 9, "Hill climb v4f"), with artboards Main, MedialAxis, StimulusLab,
Tracing and FirstRun. Direction and product context are in `PRODUCT.md` and
`.impeccable/surfaces/studio.md`; the approved composition is
`.impeccable/mocks/comp-2.png` (canvas-dominant).

Design decisions that change §2:

- **Field:** calibrated grey (#808080) chrome around the stimulus; red
  (#D2342B) is reserved for the fixation/selection accent. Type is Atkinson
  Hyperlegible Next and Atkinson Hyperlegible Mono.
- **Every displayed number comes from R.** The sketch was filled from real
  mlvtoolbox 0.2.1 output for a synthetic drawing; the implementation keeps
  that rule (principle 1).
- **Medial axis:** a Coverage view showing which contour samples receive a
  score, per contour and pooled. Open-background branches are hidden by
  default and excluded from ranking and coverage, labelled as a screen rule
  until the package defines one (see below).
- **Stimulus lab:** the verdict leads with the outcome at default tolerances;
  edited tolerances never hide it. Invalid tolerances block export. The export
  manifest records measured values, the tolerances used and the
  default-tolerance outcome.
- **Copy** actions always report success or failure; selections are announced
  to screen readers.
- **Global ⌘K** was not demonstrable in the canvas format and remains a
  requirement for the implementation.

Known gaps carried into implementation (none blocking): focus is lost when
the contour table opens; the contour-ID legend is approximate; the lab does
not say when no fraction passes at default tolerances, nor mark per-row which
checks fail at defaults; legends lack tick marks; "Accept and inspect" in
Tracing should open the traced drawing.

### 8.1 Package questions raised by the sketch

The sketch had to compute some things in the screen layer. Principle 1 says
they belong in the package, or not in the UI:

1. **Open-background MAT branches.** The sketch classifies a branch as
   background when more than half its points lie in free space 4-connected to
   the border of `rasterize(x)`. On house1 this hides 99 of 195 branches but
   also removes open parts of objects (the walkway; the trunk, whose edge
   stops 5 px short of the crown). The package should either provide a
   documented classifier or the UI should not filter.
2. **MAT mapping coverage.** Coverage per contour now exists inside
   measurement metadata (`R/measurement.R`); the UI needs it as a public,
   documented result for a MAT, both for all branches and for a branch subset.
3. **Split feasibility.** On house1 no `split_contours()` fraction passes all
   four matching checks at default tolerances. The UI needs a package-level
   answer to "is any fraction feasible?" (`design_subsets()` already reports
   `infeasible`; `split_contours()` does not).
4. **Export manifest schema.** Fields: inputs and provenance, criterion,
   fraction, tolerances used, measured check values, default-tolerance outcome,
   package version. Should be defined alongside `feature_spec` so the UI and
   scripted exports write the same record.
