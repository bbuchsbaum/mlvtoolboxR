# Purpose, lineage and differences from MLV

`mlvtoolbox` is an independently maintained R implementation and extension of
the Mid-Level Vision Toolbox. Its purpose is to let vision researchers measure
contour structure, inspect the resulting geometry, and construct stimulus sets
whose chosen properties and matching constraints can be reproduced in R.

The scientific starting point is the original MLV Toolbox. This port keeps its
published estimators available, documents numerical conventions, fixes
demonstrated defects, adds explicitly named estimators and workflow tools,
and measures performance while preserving tested numerical behavior.
Bradley Buchsbaum maintains the R implementation. Upstream credits do not imply
that those authors maintain or endorse this port.

## Sources and people

| Contribution | Credit and source |
|---|---|
| Original MLV Toolbox and 2023 methods paper | Dirk Bernhardt-Walther (Dirk B. Walther), Delaram Farzanfar, Seohee Han, Morteza Rezanejad, and the Bernhardt-Walther Lab, University of Toronto. [MATLAB repository](https://github.com/bwlabToronto/MLV_toolbox), [paper](https://doi.org/10.3389/fcomp.2023.1140723). |
| Python reference implementation | [pyMLV](https://github.com/bwlabToronto/pyMLV), Bernhardt-Walther Lab. Aravind Narayanan is credited specifically in its AOF implementation. |
| AOF and medial-axis methods | Morteza Rezanejad's [AOFSkeletons](https://github.com/mrezanejad/AOFSkeletons) and MLV work. [Rezanejad, Downs, Wilder, Walther, Jepson, Dickinson and Siddiqi (2019)](https://openaccess.thecvf.com/content_CVPR_2019/papers/Rezanejad_Scene_Categorization_From_Contours_Medial_Axis_Based_Salience_Measures_CVPR_2019_paper.pdf). |
| Pratt morphology lookup tables | Josep Mones i Teixidor and Carnë Draug; adapted from `bwmorph.m` in [GNU Octave image 2.20.0](https://gnu-octave.github.io/packages/image/). |
| Retained MATLAB edge-linking references | Peter Kovesi, University of Western Australia. [Original notices](../inst/LICENSES/Kovesi-notices.txt) remain attached to these references; they are not R runtime code. |
| R implementation and extensions | Bradley Buchsbaum: R APIs, numerical safeguards, native kernels, feature schemas, alternative estimators, constrained subset design, graph descriptors and qualification infrastructure. |

The original toolbox also cites Walther and Shen (2014), *Nonaccidental
properties underlie human categorization of complex natural scenes*, and
Rezanejad and Siddiqi (2013), *Flux graphs for 2D shape analysis*. Those works
supply scientific context; the R port does not claim authorship of these ideas.

Pinned references: MATLAB `f12418a5e03c2f38e222fed1dffccf42b227cd29` and Python
`14b3183b91146b1021eb7ef97db75ad703b06832`. The [source manifest](../vendor/SOURCES.json)
records retained paths and hashes. The Python port shares the MATLAB lineage,
so agreement between them is differential evidence, not independent scientific
validation.

## What stays recognizable and what changes

| Area | Original/reference behavior | This R implementation |
|---|---|---|
| Core measurements | Length, orientation, curvature, junctions and medial-axis salience | Keeps these concepts and the `mlv2023` estimator family; exposes R drawing objects and tables with explicit parameter metadata. |
| Numerical conventions | Some routines truncate coordinates, mishandle boundaries or conflate missing and zero scores | Preserves subpixel geometry, validates degeneracies and missingness, and corrects documented bins and geometry. Evidence is in the [bug ledger](../UPSTREAM_BUGS.md). |
| Alternative estimators | Published contour and branch formulas | Uses vertex-turning curvature by default (published curvature via `estimator = "mlv2023"`) and adds an opt-in `arclength` estimator, calibrated physical units and fixed feature schemas. See the [measurement contract](FEATURE_SPEC.md). |
| Execution | MATLAB toolboxes and the Python stack | Uses R and Rcpp/C++17 for the ordinary pipeline. Compiled geometry and fused histograms have [workload-specific performance and equality evidence](FUSED_STATISTICS.md). |
| Raster and MAT behavior | MATLAB morphology, rasterization and boundary choices | Uses Octave-derived Pratt morphology, Bresenham rasterization, deterministic EDT ties and a pixel-graph tracer. These choices can change flux and skeletons despite equal distances. |
| Stimulus and graph workflows | Original contour splitting/manipulation | Adds [bounded constrained subset design](SUBSET_DESIGNER.md), with feasibility/search diagnostics, and [graph descriptors](CONTOUR_GRAPHS.md) for components, closure, continuation and enclosure. |
| Photographs and models | MATLAB SAM/structured-edge tooling and serialized models | Requires an explicit mask-producing provider for `trace_rgb()`. No SAM weights/backend are bundled. R models must match their feature schema; MATLAB serialized regressors are not imported. |

## Qualification boundaries

This is a port with documented departures, not a drop-in MATLAB replacement.
Actual MATLAB was not available for qualification. Executed references use
pinned MATLAB sources through Octave and selected Python functions with
explicitly recorded adapters.

[Parity evidence](PARITY.md) separates matches, upstream defects and method
differences. Nearest-boundary ties change default MAT skeletons on 4/6 test
masks; matching boundary labels resolves those six cases. Tests establish the
stated computational behavior on their fixtures, not perceptual validity,
broad natural-image equivalence, all-platform qualification or identical
performance everywhere. The tools-only Voronoi probe is a research experiment,
not a new default.

## Citation and licenses

Use `citation("mlvtoolbox")` for the R software and original MLV paper, plus
the medial-axis paper when appropriate. [CITATION](../inst/CITATION) contains
all three entries. Notices and credits are installed with the package:

~~~r
system.file("COPYRIGHTS", package = "mlvtoolbox")
system.file("LICENSES", package = "mlvtoolbox")
~~~

The R implementation is declared GPL-3.0-or-later. The [complete text](../COPYING)
and [third-party notices](../inst/LICENSES) accompany the source.
**The vendored MATLAB core has unresolved license grants:** the
[per-file audit](../vendor/LICENSING.md) records the permission/lineage review
needed before public redistribution. Our declaration does not relicense
upstream sources or datasets.
