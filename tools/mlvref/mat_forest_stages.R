# Track B: controlled-substitution stage attribution on Forest_DBW (aggregates only).
# Each row replaces stages of the R pipeline by MATLAB's stored fields or pinned-Octave branches
# and compares mapped per-pixel scores with the stored MATLAB samples at |error| <= 1e-4.
source("tools/mlvref/mat_common.R")
f <- read_forest(mlvref_root()); ink <- f$D <= 0; st <- stored_samples_forest(f$raw)
x <- read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat", zero_length = "drop")
od <- ".agent-work/track-b/oct_forest_single"; allpix <- readRDS(file.path(od, "allpix.rds")); rc <- cbind(allpix$y, allpix$x)
b <- read.csv(file.path(od, "branches.csv"), header = FALSE, col.names = c("b", "row", "col", "radius"))
edgelink <- unname(lapply(split(b, b$b), function(z) data.frame(row = z$row, column = z$col, x = z$col, y = z$row,
  radius = f$D[cbind(z$row, z$col)], AOF = f$AOF[cbind(z$row, z$col)])))
score <- function(name, m, emulate) {
  m <- compute_mat_properties(m)
  if (emulate) { v <- map_at(m, rc, matlab_image = TRUE); imgs <- lapply(v, function(z) { im <- matrix(NA_real_, 768, 1024); im[rc] <- z; im }) }
  else imgs <- map_mat_properties(x, m)$mat_images
  do.call(rbind, lapply(names(st), function(p) { r <- compare_samples(imgs[[p]], st[[p]], 1e-4)
    data.frame(stage = name, property = p, n = r[["n"]], covered = r[["covered"]], within_1e4 = r[["within_tol"]],
               median_abs = r[["median_abs"]], q90_abs = r[["q90_abs"]], cor = r[["cor"]]) }))
}
mk <- function(D, A, S, br = NULL) { m <- make_mat(D, A, S, ink); m$branches <- if (is.null(br)) trace_skeleton(m, 10L) else br; m }
thr <- asin(.3 * pi / 2) * 180 / pi
rows <- list(
  score("A pinned-equivalent: MATLAB S + edgelink branches, MLV mapping order", mk(f$D, f$AOF, f$skeleton, edgelink), TRUE),
  score("B MATLAB S + edgelink branches, R mapping", mk(f$D, f$AOF, f$skeleton, edgelink), FALSE),
  score("C MATLAB S + R tracing, MLV mapping order", mk(f$D, f$AOF, f$skeleton), TRUE),
  score("D MATLAB S + R tracing, R mapping", mk(f$D, f$AOF, f$skeleton), FALSE),
  score("E R hmt8 + mlv2023 AOF (published-compatible MAT) + R tracing/mapping", compute_mat(x, threshold_angle = thr, morphology = "hmt8", aof = "mlv2023"), FALSE),
  score("F R pratt skeleton (R raster) + R tracing/mapping", compute_mat(x, threshold_angle = thr, morphology = "pratt", aof = "mlv2023"), FALSE),
  score("G R defaults: hmt8 + centred tie-mean AOF", compute_mat(x, threshold_angle = thr, morphology = "hmt8", aof = "centered"), FALSE))
r <- do.call(rbind, rows); print(r[, c("stage", "property", "covered", "within_1e4", "median_abs", "cor")], digits = 3, row.names = FALSE)
write.csv(r, "docs/evidence/mlvref/mat/forest_stage_attribution.csv", row.names = FALSE)
