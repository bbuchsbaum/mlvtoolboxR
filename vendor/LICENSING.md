# Licenses in the pinned reference snapshots

Audit: 2026-10-26. **These snapshots have mixed terms; this is not blanket
redistribution clearance.** All 159 retained files match [SOURCES.json](SOURCES.json).
Their original headers and copyright notices are unchanged. The entire vendor
directory is excluded from the R package.

| Snapshot | Pinned revision | Observed terms |
|---|---|---|
| [MLV_toolbox](https://github.com/bwlabToronto/MLV_toolbox) | f12418a5e03c2f38e222fed1dffccf42b227cd29 | 5 GPL-3.0-or-later files; 6 Peter Kovesi files with explicit permissive notices; 82 files without an explicit grant identified |
| [pyMLV](https://github.com/bwlabToronto/pyMLV) | 14b3183b91146b1021eb7ef97db75ad703b06832 | MIT root license is the baseline for 63 retained files; 3 files retain explicit GPL-3.0-or-later notices |

The [per-file inventory](../docs/evidence/attribution/vendor-license-inventory.csv)
records paths, hashes, observed terms and attribution lines. `NOASSERTION` means
the audit did not establish a license grant; it does not mean public domain.
The scanner recognizes notices, not ownership or complete derivative lineage.
Run `python3 tools/audit_licenses.py` to verify the snapshots and regenerate it.

## Notices preserved

- [pyMLV MIT license](pyMLV/LICENSE), copyright 2024 Bernhardt-Walther Lab,
  also included verbatim in [the installed bundle](../inst/LICENSES/pyMLV-MIT.txt).
- [Explicit MLV GPL notices](../inst/LICENSES/MLV-GPL-notices.txt), including
  Morteza Rezanejad's AOF work and Aravind Narayanan's Python AOF credit.
  The MATLAB files are `computeAOF.m`, `mapMATtoContour.m`,
  `extract2DSkeletonFromBinaryImage.m`, `InitializeNeighborhoods.m` and
  `computeGradientVectorField.m`. The Python files are `computeAOF.py`,
  `mapMATtoContour.py` and `extract2DSkeletonFromBinaryImage.py`.
- [Peter Kovesi notices](../inst/LICENSES/Kovesi-notices.txt), retained verbatim
  for `edgelink.m`, `findendsjunctions.m`, `drawedgelist.m`, `lineseg.m`,
  `cleanedgelist.m` and `edgelist2image.m`. Their abbreviated permissive wording
  is recorded as `LicenseRef-Kovesi-permissive`, not silently replaced by MIT.
- [Other copyright and author lines](../inst/LICENSES/Vendored-copyrights.txt).
  A copyright credit alone does not supply a grant.
- [Complete GPL v3 text](../inst/LICENSES/GPL-3.0.txt), also at [COPYING](../COPYING).
  The R kernels additionally adapt GNU Octave image Pratt tables:
  [original notice](../inst/LICENSES/Octave-image-notice.txt),
  [package credits](../inst/COPYRIGHTS).

## Unresolved upstream permissions

**The public repository does not contain these 82 files.** They were removed
from the entire Git history before first publication and are ignored by Git.
Their paths and hashes remain in [SOURCES.json](SOURCES.json) and the
inventory. Parity, benchmark and license-audit tools need them; run
`python3 tools/fetch_ungranted_vendor.py` to fetch them from the pinned
upstream commit, verified against the recorded SHA-256 hashes.


The pinned MATLAB repository has **no root license file**. Its README directs
readers to individual file licenses. Of 82 retained files without an explicit
grant, 77 are `.m` files. This includes much of `MLVcode/`, `maxlinedev.m`,
and helpers under `3rdParty/contour_tracing/main/`. The remainder is a README,
SVG, HTML and two small MAT reference files. See the complete inventory rather
than assuming a neighboring file's notice covers them.

The [upstream tree audit](../docs/evidence/attribution/upstream-license-tree.json)
checked license-named files across both complete pinned trees. Other MATLAB
license files belong to omitted m2html, redblue and structured-edge subtrees.
Those subtrees are not retained here and their licenses do not supply a grant
for the retained MATLAB core. Python's root license is retained.

Before public redistribution of these unresolved reference files, or of code
derived from them, obtain the applicable copyright holder's permission or
complete a documented replacement/lineage review. That review must also consider
any R routines translated from copyright-only MATLAB sources; excluding
`vendor/` from the tarball does not resolve this question. Attribution, the
paper's open-access status and our GPL declaration are not a substitute for a
grant. [GitHub's licensing guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository)
explains the distinction between public source and a license.

No upstream permission request or public redistribution was performed by this
audit. It preserves the evidence and identifies the permission gap.

## Data and optional tools

Datasets, historical MAT data and weights have separate provenance and
permissions. Their presence upstream does not establish redistribution rights.
They are excluded from package archives; package tests use synthetic fixtures.
The local reference collection remains at `reference-data/mlv` under the custody
recorded in [SOURCES.json](SOURCES.json).

The Boost.Polygon probe under `tools/voronoi/` uses separately installed Boost/BH
headers and is excluded from the R package. No Boost headers, MATLAB runtime or
SAM weights are vendored in the package.
