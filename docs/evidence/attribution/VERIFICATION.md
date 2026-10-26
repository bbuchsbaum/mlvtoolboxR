# Attribution audit verification

The audit verifies notice preservation and package contents, not blanket
permission to redistribute the port or the vendored files.

- All 159 retained vendor files match the original snapshot manifest.
- Pinned upstream tree inspection found no MATLAB root license; the Python
  root MIT license was retained. All 8 explicit GPL notices and 6 Kovesi
  permission notices are preserved in the installed notice bundle.
- Complete GPL and MIT texts, original notice extracts, author lines and
  citations are included in the source tarball and installed package.
  Exact SHA-256 comparisons passed for all 8 license/credit/citation files.
- Three installed citation entries parse successfully; the README drawing,
  feature-table and plotting example executes successfully.
- Full R CMD check --as-cran --no-manual: 0 errors, 0 warnings, 2 notes
  (New submission and inability to verify current time). 1,332 assertions pass
  without test failures, warnings or skips.
- Runtime R/native code, tests, numerical definitions and vendor bytes were
  unchanged. Coverage and benchmarks were not rerun for this documentation and
  metadata change.

Archive: `mlvtoolbox_0.2.1-attribution.tar.gz`, SHA-256
`7a4497651076d244319c1a83486148cc349ca8806c59d0a1c650efddc00038ca`.
All package inputs are bound by [the source manifest](package-source-manifest.json).
The check used a frozen package over `8a1452b3af4d8953c1edc6d7f53a49f771f2ab54`; that parent
alone does not include these metadata changes. The later checkpoint receipt
binds the tested input files to their committed blobs. Concurrent studio
housekeeping in the live .Rbuildignore was preserved separately.

[Package inspection](packaging.json), [check log](00check.log),
[test output](testthat.Rout), [citation output](citation.txt),
[vendor inventory](vendor-license-inventory.csv),
[notice provenance](notice-sources.json), [upstream tree check](upstream-license-tree.json).

Unresolved: 82 MATLAB reference files, including 77 .m files, have no identified
explicit grant. Obtain appropriate permission or complete a documented
replacement/lineage review before publicly redistributing those files or
derived code. See [licensing audit](../../../vendor/LICENSING.md).
