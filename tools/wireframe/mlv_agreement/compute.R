# MLV statistics for artist and automatic Toronto drawings (frozen plan:
# docs/evidence/wireframe/mlv-agreement/PLAN.md).
#
# Usage: Rscript --vanilla tools/wireframe/mlv_agreement/compute.R OUT_DIR
# Run from the package root. Writes histograms.csv and scalars.csv.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Expected the output directory")
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(".", quiet = TRUE)

dev <- read.csv("docs/evidence/wireframe/curve-rank/development.csv", colClasses = c(id = "character"))
test <- read.csv("docs/evidence/wireframe/curve-rank/reserve-test.csv", colClasses = c(id = "character"))
roster <- rbind(data.frame(id = dev$id, category = dev$category, set = "development"),
                data.frame(id = test$id, category = test$category, set = "test"))
stopifnot(!anyDuplicated(roster$id), nrow(roster) == 180L)

methods <- c("fdog", "mlsd", "union_vector", "curve_rank")
source_file <- function(id, set, method) {
  if (set == "development") {
    if (method == "curve_rank") file.path(".agent-work/curve-rank/development", paste0(id, ".curve_rank.contours.json"))
    else file.path(".agent-work/contour-study/draw-toronto", paste0(id, ".", method, ".contours.json"))
  } else {
    if (method == "curve_rank") file.path(".agent-work/curve-rank/reserve/predicted", paste0(id, ".curve_rank.contours.json"))
    else file.path(".agent-work/curve-rank/reserve/base", paste0(id, ".", method, ".contours.json"))
  }
}
staged_file <- function(id, set) {
  dir <- if (set == "development") ".agent-work/toronto/staged" else ".agent-work/curve-rank/reserve/staged"
  file.path(dir, paste0(id, ".json"))
}

# Inverse of scaled_point() in tools/wireframe/toronto_scenes.py.
to_native <- function(m, W, H, w, h) {
  m[, c(1, 3)] <- (m[, c(1, 3)] - 0.5) * W / w + 0.5
  m[, c(2, 4)] <- (m[, c(2, 4)] - 0.5) * H / h + 0.5
  m
}
artist_segments <- function(refs) {
  lapply(refs$polylines, function(p) {
    xy <- matrix(unlist(p$points), ncol = 2, byrow = TRUE)
    cbind(head(xy, -1), tail(xy, -1))
  })
}
method_segments <- function(path) {
  js <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  keep <- Filter(function(k) length(k$segments) > 0L, js$contours)
  lapply(keep, function(k) matrix(unlist(k$segments), ncol = 4, byrow = TRUE))
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

hist_rows <- list(); scalar_rows <- list(); missing <- character()
for (i in seq_len(nrow(roster))) {
  id <- roster$id[i]; set <- roster$set[i]
  refs <- jsonlite::fromJSON(staged_file(id, set), simplifyVector = FALSE)
  stopifnot(identical(refs$id, id))
  W <- refs$source$original_width; H <- refs$source$original_height
  w <- refs$working$width; h <- refs$working$height
  segs <- list(artist = artist_segments(refs))
  for (m in methods) {
    f <- source_file(id, set, m)
    if (!file.exists(f)) { missing <- c(missing, f); next }
    segs[[m]] <- method_segments(f)
  }
  for (m in names(segs)) {
    for (config in c("native", "working", "published")) {
      s <- segs[[m]]
      if (config != "working") s <- lapply(s, to_native, W = W, H = H, w = w, h = h)
      dims <- if (config == "working") c(w, h) else c(W, H)
      d <- mlv_drawing(s, dims[1], dims[2], zero_length = "drop")
      st <- if (config == "published") contour_statistics(d, spec = published(W, H)) else contour_statistics(d)
      hist_rows[[length(hist_rows) + 1L]] <- data.frame(id = id, category = roster$category[i], set = set,
        method = m, config = config, property = st$property, bin = st$bin, label = st$label, value = st$value)
      if (config != "published") {
        sc <- scalars(d, st)
        scalar_rows[[length(scalar_rows) + 1L]] <- data.frame(id = id, category = roster$category[i], set = set,
          method = m, config = config, name = names(sc), value = unname(sc))
      }
    }
  }
  if (i %% 20L == 0L) message(i, " photos")
}
if (length(missing)) stop("missing drawings: ", paste(head(missing), collapse = ", "))
write.csv(do.call(rbind, hist_rows), file.path(out, "histograms.csv"), row.names = FALSE)
write.csv(do.call(rbind, scalar_rows), file.path(out, "scalars.csv"), row.names = FALSE)
message("done")
