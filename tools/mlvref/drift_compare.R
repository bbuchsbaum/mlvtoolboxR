# P1 version drift: compare the stored MATLAB fields with the pinned MLV
# source re-executed under Octave (tools/mlvref/drift_pinned.m) on the same
# stored contours. No mlvtoolbox code is involved.
# Usage: Rscript tools/mlvref/drift_compare.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else ".agent-work/track-a/lib"
commandArgs <- function(trailingOnly = TRUE) c(lib, "pinned-vs-stored", "1")
src <- readLines("tools/mlvref/compare_junctions.R"); n <- grep("^sets <- complete_sets", src)
eval(parse(text = src[1:(n - 1)]))
rm(commandArgs)
sets <- intersect(complete_sets(export_root), complete_sets(pinned_root))
rows <- lapply(sets, function(ds) {
  sample <- utils::read.csv(file.path(pinned_root, ds, "sample.csv"))
  ps <- read_gz(file.path(pinned_root, ds, "segments.csv"))
  pj <- read_gz(file.path(pinned_root, ds, "junctions.csv"))
  ss <- read_gz(file.path(export_root, ds, "segments.csv"))
  sj <- read_gz(file.path(export_root, ds, "junctions.csv"))
  ss <- ss[ss$drawing %in% sample$drawing, ]
  i <- match(paste(ps$drawing, ps$contour, ps$segment), paste(ss$drawing, ss$contour, ss$segment))
  stopifnot(!anyNA(i), length(i) == nrow(ss))
  s <- ss[i, ]
  same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & (a == b | abs(a - b) <= tol_geom * pmax(1, abs(b))))
  od <- same(ps$orientation, s$orientation); ld <- same(ps$length, s$length); cd <- same(ps$curvature, s$curvature)
  jstats <- lapply(sample$drawing, function(k) {
    A <- with_minmax(sj[sj$drawing == k, ]); B <- with_minmax(pj[pj$drawing == k, ])
    m <- compare_one(A, B, k)$matched
    ex <- if (is.null(m)) 0 else sum(m$pos_ok & m$type_a == m$type_b & m$pairs_ok & m$angles_ok & m$minmax_ok)
    c(nA = nrow(A), nB = nrow(B), exact = ex)
  })
  js <- colSums(do.call(rbind, jstats))
  data.frame(dataset = ds, drawings = nrow(sample), segments = nrow(ps),
    orientation_equal = sum(od), length_equal = sum(ld), curvature_equal = sum(cd),
    orientation_bitwise = sum(ps$orientation == s$orientation, na.rm = TRUE),
    length_bitwise = sum(ps$length == s$length, na.rm = TRUE),
    curvature_bitwise = sum(ps$curvature == s$curvature, na.rm = TRUE),
    orientation_max_abs = max(abs(ps$orientation - s$orientation), na.rm = TRUE),
    curvature_max_rel = max(abs(ps$curvature - s$curvature) / pmax(1, abs(s$curvature)), na.rm = TRUE),
    stored_junctions = js[["nA"]], pinned_junctions = js[["nB"]], junctions_identical = js[["exact"]],
    octave_seconds = sum(sample$seconds))
})
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(receipt_root, "drift-pinned-vs-stored.csv"), row.names = FALSE)
print(colSums(res[, -1])); print(sapply(res[, c("orientation_max_abs", "curvature_max_rel")], max))
