# Natural-data receipt for two 0.3.0 changes: the specification-free curvature
# default (estimator = "turning" instead of "mlv2023") and bin-edge snapping.
# For every stored drawing it records whether the turning profile is defined
# (exact retracing is accepted with an unsigned 180-degree turn), the
# integrated absolute curvature under both estimators, the L1 distance between
# their normalized 8-bin curvature histograms, and how many orientation,
# junction-angle, length and curvature values snapping moves to another bin.
# Usage: Rscript tools/mlvref/compare_curvature_default.R <lib> [cores]
args <- commandArgs(TRUE)
.libPaths(c(normalizePath(args[1]), .libPaths()))
library(mlvtoolbox)
cores <- if (length(args) > 1) as.integer(args[2]) else 4L
export_root <- ".agent-work/mlvref"
receipt_root <- "docs/evidence/mlvref/curvature-default"
dir.create(receipt_root, recursive = TRUE, showWarnings = FALSE)

exact_num <- Rcpp::cppFunction('NumericVector exact_num(CharacterVector s) {
  NumericVector out(s.size());
  for (R_xlen_t i = 0; i < s.size(); ++i) {
    if (CharacterVector::is_na(s[i]) || !std::strlen(s[i])) { out[i] = NA_REAL; continue; }
    out[i] = std::strtod(s[i], nullptr);
  }
  return out;
}', includes = "#include <cstdlib>\n#include <cstring>")
read_gz <- function(path) {
  z <- utils::read.csv(gzfile(paste0(path, ".gz")), colClasses = "character")
  for (v in intersect(names(z), c("drawing", "contour", "segment", "width", "height"))) z[[v]] <- as.integer(z[[v]])
  for (v in intersect(names(z), c("x1", "y1", "x2", "y2"))) z[[v]] <- exact_num(z[[v]])
  z
}
zero_seg <- function(p) {
  tol <- function(a, b) 64 * .Machine$double.eps * pmax(abs(a), abs(b))
  abs(p$x2 - p$x1) <= tol(p$x1, p$x2) & abs(p$y2 - p$y1) <= tol(p$y1, p$y2)
}
build_drawing <- function(s, width, height) {
  s <- s[order(s$contour, s$segment), ]
  keep <- lapply(split(s, s$contour), function(p) p[!zero_seg(p), ])
  keep <- keep[vapply(keep, nrow, 1L) > 0L]
  mlv_drawing(unname(lapply(keep, function(p) as.matrix(p[, c("x1", "y1", "x2", "y2")]))), width, height, ids = names(keep))
}
moved <- function(values, edges) {
  values <- values[is.finite(values)]
  sum(findInterval(values, edges, all.inside = TRUE) != mlvtoolbox:::.bin_index(values, edges))
}
one <- function(d) {
  bins <- 8L
  st <- mlvtoolbox:::.segment_table(d)
  p <- tryCatch(as.data.frame(d), error = function(e) conditionMessage(e))
  if (is.character(p)) return(data.frame(turning_error = p))
  rev <- sum(is.na(p$signed_curvature) & p$curvature > 0)
  h <- function(e) contour_statistics(d, "curvature", bins = bins, normalize = TRUE, estimator = e)$value
  axial <- st$orientation %% 180
  old_o <- (floor(axial * bins / 180 + .5) %% bins) + 1
  j <- detect_junctions(d)
  data.frame(turning_error = NA_character_, segments = nrow(st), reversal_pieces = rev,
    integrated_mlv2023 = sum(st$curvature * st$length), integrated_turning = sum(p$curvature * p$length),
    hist_l1 = sum(abs(h("mlv2023") - h("turning"))),
    orientation_moved = sum(old_o != mlvtoolbox:::.orientation_index(st$orientation, bins)),
    junctions = nrow(j), junction_angle_moved = moved(j$min_angle, seq(0, 120, length.out = bins + 1L)),
    curvature_moved = moved(log10(st$curvature + 1), seq(0, log10(91), length.out = bins + 1L)),
    turning_curvature_moved = moved(p$curvature, 10^seq(0, log10(91), length.out = bins + 1L) - 1))
}
sets <- basename(dirname(Sys.glob(file.path(export_root, "*", "COMPLETE"))))
rows <- list()
for (ds in sets) {
  segs <- read_gz(file.path(export_root, ds, "segments.csv"))
  draws <- read_gz(file.path(export_root, ds, "drawings.csv"))
  out <- parallel::mclapply(draws$drawing, function(k) {
    r <- one(build_drawing(segs[segs$drawing == k, ], draws$width[draws$drawing == k], draws$height[draws$drawing == k]))
    cbind(dataset = ds, drawing = k, r)
  }, mc.cores = cores)
  bad <- vapply(out, inherits, NA, "try-error")
  if (any(bad)) stop(ds, ": ", paste(unique(vapply(out[bad], as.character, "")), collapse = "\n"))
  rows[[ds]] <- do.call(rbind, lapply(out, function(o) {for (v in setdiff(c("segments", "reversal_pieces",
    "integrated_mlv2023", "integrated_turning", "hist_l1", "orientation_moved", "junctions", "junction_angle_moved",
    "curvature_moved", "turning_curvature_moved"), names(o))) o[[v]] <- NA; o}))
  message(ds, ": ", nrow(draws), " drawings")
}
all <- do.call(rbind, rows)
by_set <- do.call(rbind, lapply(split(all, all$dataset), function(z) data.frame(dataset = z$dataset[1],
  drawings = nrow(z), turning_errors = sum(!is.na(z$turning_error)), segments = sum(z$segments, na.rm = TRUE),
  drawings_with_reversal = sum(z$reversal_pieces > 0, na.rm = TRUE),
  integrated_ratio_median = stats::median(z$integrated_mlv2023 / z$integrated_turning, na.rm = TRUE),
  hist_l1_median = stats::median(z$hist_l1, na.rm = TRUE),
  orientation_moved = sum(z$orientation_moved, na.rm = TRUE), junctions = sum(z$junctions, na.rm = TRUE),
  junction_angle_moved = sum(z$junction_angle_moved, na.rm = TRUE),
  curvature_moved = sum(z$curvature_moved, na.rm = TRUE), turning_curvature_moved = sum(z$turning_curvature_moved, na.rm = TRUE))))
total <- data.frame(dataset = "ALL", drawings = nrow(all), turning_errors = sum(!is.na(all$turning_error)),
  segments = sum(all$segments, na.rm = TRUE), drawings_with_reversal = sum(all$reversal_pieces > 0, na.rm = TRUE),
  integrated_ratio_median = stats::median(all$integrated_mlv2023 / all$integrated_turning, na.rm = TRUE),
  hist_l1_median = stats::median(all$hist_l1, na.rm = TRUE),
  orientation_moved = sum(all$orientation_moved, na.rm = TRUE), junctions = sum(all$junctions, na.rm = TRUE),
  junction_angle_moved = sum(all$junction_angle_moved, na.rm = TRUE),
  curvature_moved = sum(all$curvature_moved, na.rm = TRUE), turning_curvature_moved = sum(all$turning_curvature_moved, na.rm = TRUE))
summary <- rbind(by_set, total)
utils::write.csv(summary, file.path(receipt_root, "summary.csv"), row.names = FALSE)
print(summary)
