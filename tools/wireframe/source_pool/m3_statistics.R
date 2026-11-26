# MLV statistics of the photo-only M3 drawings from the source screen
# (development photos only; exploratory, outside the frozen screening rule).
#
# Usage: Rscript --vanilla tools/wireframe/source_pool/m3_statistics.R \
#   .agent-work/source-pool .agent-work/source-pool/m3-selection.csv OUT_DIR
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Expected the source-pool folder, the selection CSV and an output folder")
pkgload::load_all(".", quiet = TRUE)
src <- args[1]; out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sel <- read.csv(args[2], colClasses = c(id = "character", curves = "character"))
dev <- read.csv("docs/evidence/wireframe/curve-rank/development.csv", colClasses = c(id = "character"))
stopifnot(all(sel$id %in% dev$id))
to_native <- function(m, W, H, w, h) {
  m[, c(1, 3)] <- (m[, c(1, 3)] - 0.5) * W / w + 0.5
  m[, c(2, 4)] <- (m[, c(2, 4)] - 0.5) * H / h + 0.5
  m
}
scalars <- function(drawing, stats) {
  seg <- compute_contour_properties(drawing)$segment_properties
  len <- seg$length; total <- sum(len)
  jt <- stats[stats$property == "junction_type", ]
  dens <- setNames(if (total > 0) jt$value / total * 1e4 else rep(NA_real_, nrow(jt)), paste0("density_", jt$label))
  o <- seg$orientation %% 90
  hv <- if (total > 0) sum(len[pmin(o, 90 - o) <= 22.5]) / total else NA_real_
  curv <- if (total > 0) sum(len * abs(seg$curvature), na.rm = TRUE) / sum(len[is.finite(seg$curvature)]) else NA_real_
  c(total_length = total, n_contours = length(drawing$contours), dens, mean_abs_curvature = curv, hv_fraction = hv)
}
hist_rows <- list(); scalar_rows <- list()
for (k in seq_len(nrow(sel))) {
  id <- sel$id[k]; m <- sel$source[k]
  refs <- jsonlite::fromJSON(file.path(".agent-work/toronto/staged", paste0(id, ".json")), simplifyVector = FALSE)
  W <- refs$source$original_width; H <- refs$source$original_height
  w <- refs$working$width; h <- refs$working$height
  js <- jsonlite::fromJSON(file.path(src, "draw", paste(id, m, "contours.json", sep = ".")), simplifyVector = FALSE)
  take <- as.integer(strsplit(sel$curves[k], " ", fixed = TRUE)[[1]])
  segs <- lapply(js$contours[take], function(x) matrix(unlist(x$segments), ncol = 4, byrow = TRUE))
  working <- mlv_drawing(segs, w, h, zero_length = "drop")
  png::writePNG(1.0 * rasterize(working), file.path(out, paste(id, m, "m3.png", sep = ".")))
  d <- mlv_drawing(lapply(segs, to_native, W = W, H = H, w = w, h = h), W, H, zero_length = "drop")
  st <- contour_statistics(d)
  label <- paste0(m, "_m3")
  hist_rows[[k]] <- data.frame(id = id, category = refs$category, set = "development", method = label,
    config = "native", property = st$property, bin = st$bin, label = st$label, value = st$value)
  sc <- scalars(d, st)
  scalar_rows[[k]] <- data.frame(id = id, category = refs$category, set = "development", method = label,
    config = "native", name = names(sc), value = unname(sc))
}
write.csv(do.call(rbind, hist_rows), file.path(out, "histograms-sources.csv"), row.names = FALSE)
write.csv(do.call(rbind, scalar_rows), file.path(out, "scalars-sources.csv"), row.names = FALSE)
message("done")
