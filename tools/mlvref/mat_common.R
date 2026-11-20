# Shared helpers for the track-B (medial axis) reference qualification.
# Natural data are read from reference-data/ and never written to the repository.
suppressPackageStartupMessages({
  lib <- Sys.getenv("MLVREF_LIB", "")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  library(mlvtoolbox)
})
mlvref_root <- function() {
  r <- Sys.getenv("MLVREF_ROOT", "")
  if (nzchar(r)) return(r)
  normalizePath(".")  # run from the repository root
}
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("mlvtoolbox")

# MLV signed-distance construction on an already-thinned ink raster (no re-thinning).
# Mirrors computeGradientVectorField.m and R/mat.R lines after skeletonize().
mat_field <- function(thin_ink) {
  h <- nrow(thin_ink); w <- ncol(thin_ink); white <- !thin_ink
  outer <- matrix(FALSE, h, w); ink_nb <- matrix(FALSE, h, w)
  rows <- 2:(h - 1L); cols <- 2:(w - 1L)
  for (dc in -1:1) for (dr in -1:1) if (dr != 0L || dc != 0L) {
    outer[rows, cols] <- outer[rows, cols] | (thin_ink[rows, cols] & white[rows + dr, cols + dc])
    ink_nb[rows, cols] <- ink_nb[rows, cols] | thin_ink[rows + dr, cols + dc]
  }
  outer <- outer & ink_nb
  d1 <- distance_transform(thin_ink); d2 <- distance_transform(white | outer)
  D <- d1$distance - d2$distance
  idx <- d1$index; on_ink <- d1$distance == 0
  idx[on_ink] <- d2$index[on_ink]; idx[outer] <- which(outer)
  list(D = D, index = idx, outer = outer, ink = thin_ink,
       sq = ifelse(on_ink, d2$distance^2, d1$distance^2))
}
mat_skeleton <- function(flux, threshold = 0.3, samples = 60L, min_component) {
  # MATLAB: values < t set 0, > t set 1, == t keep (nonzero) -> true.
  skeletonize(flux >= threshold * samples, min_component)
}
read_forest <- function(root) {
  x <- R.matlab::readMat(file.path(root, "reference-data/mlv/examples/Forest_DBW.mat"), fixNames = FALSE)
  M <- x$MAT
  list(raw = x, skeleton = M[["skeleton", 1, 1]] != 0, D = M[["distance_map", 1, 1]],
       AOF = M[["AOF", 1, 1]])
}
fround <- function(x) { y <- readBin(writeBin(as.vector(x), raw(), size = 4), "double", size = 4, n = length(x)); dim(y) <- dim(x); y }
load_tie_kernels <- function(root = mlvref_root()) {
  if (!exists("tie_sets", mode = "function"))
    Rcpp::sourceCpp(file.path(root, "tools/mlvref/mat_ties.cpp"), env = globalenv(),
                    cacheDir = file.path(root, ".agent-work/track-b/cpp-cache"))
}
field_ties <- function(fld) {
  cls <- matrix(ifelse(fld$ink, 1L, 0L), nrow(fld$ink)); cls[fld$outer] <- 2L
  tie_sets(fld$ink, !fld$ink | fld$outer, fld$sq, cls)
}
# Choose one label per pixel from tie sets by a rule.
tie_rule_labels <- function(ts, dims, rule = c("lowest", "highest")) {
  rule <- match.arg(rule); n <- length(ts$ptr) - 1L
  len <- diff(ts$ptr); pos <- if (rule == "lowest") ts$ptr[-length(ts$ptr)] + 1L else ts$ptr[-1L]
  out <- rep(NA_integer_, n); ok <- len > 0; out[ok] <- ts$lab[pos[ok]]
  matrix(out, dims[1], dims[2])
}
# Build an mlv_mat from externally supplied stage fields (controlled substitution).
make_mat <- function(D, AOF, skeleton, mask, min_branch_points = 10L) {
  structure(list(distance_map = D, nearest_boundary = NULL, AOF = AOF, skeleton = skeleton,
    mask = mask, branches = list(), parameters = list(threshold_angle = NA, threshold = 0.3,
    samples = 60L, epsilon = 1, min_component = max(floor(.025 * max(dim(mask))), 1L),
    min_branch_points = min_branch_points, morphology = "substituted", tracing = "pixel_graph_no_corner_cut")),
    class = "mlv_mat")
}
# Stored MATLAB per-pixel samples for one drawing: list(prop -> data.frame(x,y,score)).
stored_samples_forest <- function(raw) {
  v <- raw$vecLD; nm <- dimnames(v)[[1]]; g <- function(f) v[[which(nm == f)]]
  out <- list()
  for (p in c("parallelism", "separation", "mirror"))
    out[[p]] <- data.frame(x = as.vector(g(paste0(p, "_allX"))), y = as.vector(g(paste0(p, "_allY"))),
                           score = as.vector(g(paste0(p, "_allScores"))))
  out
}
compare_samples <- function(img, st, tol = 1e-6) {
  v <- img[cbind(st$y, st$x)]
  e <- abs(v - st$score)
  c(n = nrow(st), covered = mean(is.finite(v)), within_tol = mean(e <= tol, na.rm = TRUE),
    median_abs = median(e, na.rm = TRUE), q90_abs = unname(quantile(e, .9, na.rm = TRUE)),
    max_abs = max(e, na.rm = TRUE), cor = suppressWarnings(cor(v, st$score, use = "complete.obs")))
}
# R mapping rule (map_mat_properties) evaluated at arbitrary pixel positions (row, col).
# matlab_image = TRUE emulates mapMATtoContour.m: the score is read from the skeleton rating
# image, in which later branches overwrite shared pixels (last writer wins).
map_at <- function(mat, rc, radius = 4, matlab_image = FALSE, matlab_order = matlab_image) {
  props <- intersect(c("parallelism", "separation", "mirror", "taper"), names(mat$branches[[1]]))
  out <- lapply(props, function(p) rep(NA_real_, nrow(rc))); names(out) <- props
  h <- nrow(mat$skeleton); simg <- list()
  if (matlab_image) for (p in props) { im <- matrix(0, h, ncol(mat$skeleton))
    for (b in mat$branches) im[cbind(b$row, b$column)] <- b[[p]]; simg[[p]] <- im }
  for (branch in mat$branches) {
    fp <- if (matlab_order) tangent_matlab(branch) else ns$.tangent_points(branch); if (!nrow(fp$points)) next
    lo <- apply(fp$points, 2, min) - radius; hi <- apply(fp$points, 2, max) + radius
    cand <- which(rc[, 1] >= lo[1] & rc[, 1] <= hi[1] & rc[, 2] >= lo[2] & rc[, 2] <= hi[2]); if (!length(cand)) next
    nn <- ns$cpp_nearest(fp$points, rc[cand, , drop = FALSE]); w <- which(nn$distance < radius); if (!length(w)) next
    ok <- cand[w]; src <- fp$source[nn$index[w]]
    for (p in props) { v <- if (matlab_image) simg[[p]][cbind(branch$row[src], branch$column[src])] else branch[[p]][src]
      good <- is.finite(v); at <- ok[good]; v <- v[good]
      old <- out[[p]][at]; out[[p]][at] <- ifelse(is.na(old), v, pmax(old, v)) }
  }
  out
}
# Tangent points in exactly the order of mapMATtoContour.m/getTangentPointsContour
# ([all FP1; all FP2], loop i -> r1 -> r2), including d == 0 centre points.
tangent_matlab <- function(branch) {
  X <- branch$row; Y <- branch$column; R <- branch$radius; n <- length(X)
  if (n < 2) return(list(points = matrix(numeric(), 0, 2), source = integer()))
  P1 <- P2 <- matrix(NA_real_, 9 * (n - 1), 2); src <- integer(9 * (n - 1)); k <- 0
  for (i in seq_len(n - 1)) { x1 <- X[i]; x2 <- X[i + 1]; y1 <- Y[i]; y2 <- Y[i + 1]
    for (r1 in c(R[i] - .1, R[i], R[i] + .1)) for (r2 in c(R[i + 1] - .1, R[i + 1], R[i + 1] + .1)) {
      k <- k + 1; src[k] <- i; d <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
      if (d == 0) { P1[k, ] <- P2[k, ] <- c(x1, y1); next }
      ca <- (r2 - r1) / d; if (abs(ca) > 1) ca <- ca / (abs(ca) + .00001); sa <- sqrt(1 - ca^2)
      mx <- .5 * x1 + .5 * x2; my <- .5 * y1 + .5 * y2; mr <- .5 * r1 + .5 * r2; vx <- x1 - mx; vy <- y1 - my
      f1 <- c(ca * vx + sa * vy, -sa * vx + ca * vy); f2 <- c(ca * vx - sa * vy, sa * vx + ca * vy)
      P1[k, ] <- c(mx, my) + f1 * mr / sqrt(sum(f1^2)); P2[k, ] <- c(mx, my) + f2 * mr / sqrt(sum(f2^2)) } }
  list(points = rbind(P1, P2), source = c(src, src))
}
