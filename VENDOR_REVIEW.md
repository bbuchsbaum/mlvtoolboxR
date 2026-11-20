# Vendored repository size and redundancy review

Date: 2026-10-25. This is a review; no repository files were removed,
moved, rewritten, or repacked.

## Findings

The two local reference checkouts total **3.281 GiB** in logical file
sizes. **1.917 GiB** is in their Git stores; the checked-out files total
**1.364 GiB**.

| Checkout        | Working files (MiB) | Git store (MiB) | Total (MiB) |
|-----------------|--------------------:|----------------:|------------:|
| `MLV_toolbox/`  |              712.18 |         1240.01 |     1952.18 |
| `vendor/pyMLV/` |              684.92 |          722.89 |     1407.81 |

Both are complete, non-shallow Git clones. The MATLAB checkout has 459
commits reachable from HEAD; Python has 16. Their tracked and untracked
status checks were clean. Their exact revisions and origins are recorded
in
[VENDOR_DUPLICATES.json](https://bbuchsbaum.github.io/mlvtoolboxR/VENDOR_DUPLICATES.json).

The largest individual files are the pack files: approximately 1,239.27
MiB for MATLAB and 722.77 MiB for Python. They hold Git objects for the
checkouts and their histories; they are not disposable cache files.
Deleting a pack file would damage its repository.

The MATLAB checkout is currently at `MLV_toolbox/`; Python is at
`vendor/pyMLV/`. These are the paths actually inspected.

## Confirmed duplicates

SHA-256 checks of same-size regular-file candidates at least 1 MiB found
**33 exact duplicate groups**, representing **661.04 MiB** of repeated
logical working-tree data. All copies have distinct file identities;
none are existing hard links.

| Repeated material | Duplicate groups | One-copy reduction (MiB) | Assessment |
|----|---:|---:|----|
| All dataset MAT files across both repositories | 30 | 650.54 | Strongest consolidation candidate: retain one complete shared dataset copy |
| `demos/Forest_DBW.mat` across both repositories | 1 | 7.74 | Keep one historical reference copy and preserve provenance |
| `images/example_output.png` and `images/medial_axis_transform.png` | 2 | 2.76 | Repeated illustrations; small savings |
| Python `MLVcode/MAT.mat` versus `demos/MAT.mat` | 1 payload match | 7.78 | Full hashes differ, but all bytes after their 128-byte MATLAB-v5 headers match |

The dataset duplication includes IAPS, OASIS, TorontoScenes,
ArchitecturalScenes, and Snodgrass/Vanderwart. This conclusion comes
from checksums, not just matching names or sizes. The largest pair is
`dataSets/TorontoScenes/cities_vecLD.mat`, at 41.87 MiB per copy.

The two Python `MAT.mat` files are each 8,160,320 bytes. Their
complete-file SHA-256 values differ, but both have MATLAB-v5 headers and
the payload SHA-256 is identical:
`683382d03752b8b7259ab34f7203dec7a4f51df726c0c57c637c67ee045898ac`. Thus
header differences explain why ordinary full-file checksums did not
group them. The code explicitly writes `MAT.mat` from `computeMAT()`;
retain at least one copy until its provenance and regeneration behavior
are documented.

Keeping one instance of every confirmed exact/payload duplicate would
reduce logical working-tree content by about **668.82 MiB**. Removing
tracked paths from a checkout does not remove their existing objects
from Git packs.

The exact-duplicate scan was limited to files at least 1 MiB. Smaller
repeated documentation, source assets and metadata were not exhaustively
checksum-tested.

## Additional material to slim or exclude

| Material | Observed size | Role and recommendation |
|----|---:|----|
| MATLAB `MLVcode/edges-master/` | 33.99 MiB | Legacy structured-edge backend, including a trained forest model and a large bundled toolbox. Defer from a minimal analysis snapshot; preserve as an optional backend reference |
| MATLAB `demos/` | 21.54 MiB | Mix of useful inputs, live scripts and generated figures; keep scripts/inputs selectively after identifying dependencies |
| Python `demos/` | 19.96 MiB | Includes the repeated MAT/Forest artifacts and notebooks; retain notebooks as code references |
| Python notebook embedded outputs | About 4.38 MiB serialized | Almost all notebook bulk is saved output/images; clear outputs in a derived snapshot, retaining original inputs and cells |
| Python `docs/` and MATLAB `doc/` | 1.69 + 0.53 MiB | Generated documentation; optional in a compact vendor tree |
| MATLAB `resources/` | 0.24 MiB | MATLAB project metadata; low-value storage reduction |

Python’s contour-feature notebook is 3.77 MiB and the medial-axis
notebook is 0.63 MiB. Their serialized output records occupy
approximately 4.38 MiB combined; this is an estimate from JSON
serialization rather than a measured output-cleared file-size
difference.

The first-party MATLAB `MLVcode/*.m` source is about 0.16 MiB, and the
Python code is similarly small once embedded MAT artifacts are excluded.
Most disk usage is reference data, object stores and bundled
dependencies.

Do not classify every MAT file as redundant. For example,
`vendor/pyMLV/MLVcode/imgLD.mat` is a distinct raster input, and
Python’s README says the medial-axis demo currently depends on
MATLAB-rendered input. The demonstration forest MAT file also preserves
a historical processed result; regeneration from the current source has
not been established.

## Recommended cleanup design

1.  **Keep one shared reference-data collection.** Consolidate the
    repeated dataset files and Forest example using the checksum
    manifest. Record the canonical destination and the source
    repositories/commits for each retained item. Update reference
    harness paths before removing any duplicate copy. Avoid writable
    hard links between independently editable checkouts.

2.  **Vendor pinned source snapshots for the R port.** Preserve the
    exact commit, origin, licenses and required helper files in a source
    manifest. Export only the needed code/reference artifacts into the
    eventual compact layout. If full Git history is still wanted
    locally, a shared Git object store or separately retained clone is
    another option. Shallow/sparse clones can help, but their actual
    footprint must be measured; neither guarantees removal of all large
    current-tree blobs.

3.  **Keep expensive or historical fixtures selectively.** Retain one
    copy of the MAT payload, useful rendered input, tiny synthetic
    fixtures and selected demonstrations. Preserve original specimens
    until replacements are independently reproducible. Keep legacy RGB
    weights/backends separate from ordinary contour-analysis
    development.

4.  **Exclude local references from R distribution.** Once package
    scaffolding exists, use explicit build exclusions for the vendor
    checkouts, reference data, object stores and model weights. Confirm
    the tarball contents. Exclusion affects packaging; it does not
    reclaim local disk space.

A compact source-only snapshot could remove most of the current checkout
footprint from the working project, while one shared dataset copy
remains available for optional scientific evaluation. Exact final
savings require a staged replacement, verification of its needed files
and data, and a decision about custody of the full clones. The 1.917 GiB
Git-store total is not a promise that a history-preserving cleanup can
recover that amount.

## Verification and limitations

- Read every regular file’s metadata inside the two scoped checkouts; no
  permission errors or symlinks were encountered.
- Hashed 68 candidate files, totaling about 1337.64 MiB of I/O, and
  checked file size/modification time before and after hashing.
- Examined source references and notebook output sizes without
  displaying image or dataset contents. No MAT datasets were loaded into
  the model or a scientific analysis.
- Recorded Git heads, remotes, shallow status, commit counts,
  object-store summaries and clean worktree status.
- Filesystem allocated-block totals are approximately 3.295 GiB, but
  APFS shared/cloned extents mean sums of file sizes or allocated blocks
  do not establish uniquely reclaimable space.
- The volume reported approximately 20 GiB available during review.
- No deletion, history rewrite, clone conversion, notebook clearing,
  source change, or data relocation was performed.
