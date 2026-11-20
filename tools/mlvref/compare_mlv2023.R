# Verify detect_junctions(rule = "mlv2023") against the stored MATLAB junctions
# on every exported drawing: position (<= 1e-6 px), type, contour set,
# contour:segment list (including the adjacent segments MLV appends), all
# angles and min/max angle (<= 1e-6 deg), and junction counts.
# Usage: Rscript tools/mlvref/compare_mlv2023.R <lib> stored <cores>
# Writes docs/evidence/mlvref/junctions/mlv2023-vs-stored.csv.
src <- readLines("tools/mlvref/compare_junctions.R"); n <- grep("^sets <- complete_sets", src)
eval(parse(text = src[1:(n - 1)]))
rows <- lapply(complete_sets(export_root), function(ds) {
  segs <- read_gz(file.path(export_root, ds, "segments.csv"))
  draws <- read_gz(file.path(export_root, ds, "drawings.csv"))
  sj <- read_gz(file.path(export_root, ds, "junctions.csv"))
  out <- parallel::mclapply(draws$drawing, function(k) {
    i <- which(draws$drawing == k); s <- segs[segs$drawing == k, ]
    b <- build_drawing(s, draws$width[i], draws$height[i])
    zero_drawing <- any(zero_seg(s))
    A <- with_minmax(sj[sj$drawing == k, ]); B <- with_minmax(r_junctions(b, rule = "mlv2023"))
    m <- compare_one(A, B, k)$matched
    ex <- if (is.null(m)) logical() else m$pos_ok & m$type_a == m$type_b & m$contours_ok & m$pairs_ok & m$angles_ok & m$minmax_ok
    c(drawings = 1, zero_length_drawings = zero_drawing, stored = nrow(A), r = nrow(B), exact = sum(ex),
      differing_in_zero_length_drawings = if (zero_drawing) max(nrow(A), nrow(B)) - sum(ex) else 0)
  }, mc.cores = cores)
  bad <- vapply(out, inherits, NA, "try-error"); if (any(bad)) stop(ds, ": ", as.character(out[bad][[1]]))
  r <- as.data.frame(t(colSums(do.call(rbind, out)))); r$dataset <- ds
  message(sprintf("%s: %d/%d exact", ds, r$exact, r$stored)); r
})
res <- do.call(rbind, rows); res <- res[, c("dataset", setdiff(names(res), "dataset"))]
utils::write.csv(res, file.path(receipt_root, "mlv2023-vs-stored.csv"), row.names = FALSE)
print(colSums(res[, -1]))
