# Post hoc representation control for the MLV agreement analysis.
# Rasterize each artist drawing with RES pixels on the long side and retrace it
# with the same call used for traced sources (trace_contours(simplify = 1)),
# then measure it exactly as compute.R measures the method drawings.
#
# Usage: Rscript --vanilla tools/wireframe/mlv_agreement/control.R WORK_DIR [RES]
# RES defaults to the 256-pixel working resolution. Writes method
# "artist_retraced" (RES = 256) or "artist_retraced_RES" rows to
# WORK_DIR/histograms-control[RES].csv and WORK_DIR/scalars-control[RES].csv.
args <- commandArgs(TRUE)
if (!length(args) %in% 1:2) stop("Expected the work directory and an optional resolution")
out <- args[1]
res <- if (length(args) == 2L) as.integer(args[2]) else 256L
label <- if (res == 256L) "artist_retraced" else paste0("artist_retraced_", res)
suffix <- if (res == 256L) "" else as.character(res)
pkgload::load_all(".", quiet = TRUE)
dev <- read.csv("docs/evidence/wireframe/curve-rank/development.csv", colClasses = c(id = "character"))
test <- read.csv("docs/evidence/wireframe/curve-rank/reserve-test.csv", colClasses = c(id = "character"))
roster <- rbind(data.frame(id = dev$id, category = dev$category, set = "development"),
                data.frame(id = test$id, category = test$category, set = "test"))
staged_file <- function(id, set) {
  dir <- if (set == "development") ".agent-work/toronto/staged" else ".agent-work/curve-rank/reserve/staged"
  file.path(dir, paste0(id, ".json"))
}
to_native <- function(m, W, H, w, h) {
  m[, c(1, 3)] <- (m[, c(1, 3)] - 0.5) * W / w + 0.5
  m[, c(2, 4)] <- (m[, c(2, 4)] - 0.5) * H / h + 0.5
  m
}
published <- function(W, H) feature_spec(
  properties = c("orientation", "length", "curvature", "junction_type", "junction_angle"),
  bins = 8, normalize = FALSE, length_range = c(2, W + H), estimator = "mlv2023",
  junction = list(rule = "mlv2023"))
scalars <- function(drawing, stats) {
  seg <- compute_contour_properties(drawing)$segment_properties
  len <- seg$length; total <- sum(len)
  jt <- stats[stats$property == "junction_type", ]
  dens <- setNames(if (total > 0) jt$value / total * 1e4 else rep(NA_real_, nrow(jt)), paste0("density_", jt$label))
  o <- seg$orientation %% 90
  hv <- if (total > 0) sum(len[pmin(o, 90 - o) <= 22.5]) / total else NA_real_
  curv <- if (total > 0) sum(len * abs(seg$curvature), na.rm = TRUE) / sum(len[is.finite(seg$curvature)]) else NA_real_
  c(total_length = total, n_contours = length(drawing$contours), dens,
    mean_abs_curvature = curv, hv_fraction = hv)
}
hist_rows <- list(); scalar_rows <- list()
for (i in seq_len(nrow(roster))) {
  id <- roster$id[i]; set <- roster$set[i]
  refs <- jsonlite::fromJSON(staged_file(id, set), simplifyVector = FALSE)
  W <- refs$source$original_width; H <- refs$source$original_height
  w <- refs$working$width; h <- refs$working$height
  art <- lapply(refs$polylines, function(p) {
    xy <- matrix(unlist(p$points), ncol = 2, byrow = TRUE); cbind(head(xy, -1), tail(xy, -1))
  })
  rw <- round(w * res / 256); rh <- round(h * res / 256)
  art <- lapply(art, to_native, W = rw, H = rh, w = w, h = h)
  retraced <- trace_contours(rasterize(mlv_drawing(art, rw, rh, zero_length = "drop")), simplify = 1)
  segs <- lapply(retraced$contours, as.matrix)
  for (config in c("native", "working", "published")) {
    s <- if (config == "working") lapply(segs, to_native, W = w, H = h, w = rw, h = rh)
         else lapply(segs, to_native, W = W, H = H, w = rw, h = rh)
    dims <- if (config == "working") c(w, h) else c(W, H)
    d <- mlv_drawing(s, dims[1], dims[2], zero_length = "drop")
    st <- if (config == "published") contour_statistics(d, spec = published(W, H)) else contour_statistics(d)
    hist_rows[[length(hist_rows) + 1L]] <- data.frame(id = id, category = roster$category[i], set = set,
      method = label, config = config, property = st$property, bin = st$bin, label = st$label, value = st$value)
    if (config != "published") {
      sc <- scalars(d, st)
      scalar_rows[[length(scalar_rows) + 1L]] <- data.frame(id = id, category = roster$category[i], set = set,
        method = label, config = config, name = names(sc), value = unname(sc))
    }
  }
}
write.csv(do.call(rbind, hist_rows), file.path(out, paste0("histograms-control", suffix, ".csv")), row.names = FALSE)
write.csv(do.call(rbind, scalar_rows), file.path(out, paste0("scalars-control", suffix, ".csv")), row.names = FALSE)
message("done")
