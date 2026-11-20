# Track B: is the published-compatible residual on Forest purely MLV's single precision?
# Compare the R compat chain on MATLAB's stored MAT fields with the pinned Octave chain run in
# double (oct_forest) and in single (oct_forest_single) precision, at the stored sample pixels.
source("tools/mlvref/mat_common.R")
f <- read_forest(mlvref_root()); x <- read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat", zero_length = "drop")
m <- make_mat(f$D, f$AOF, f$skeleton, f$D <= 0); m$parameters$branch_rule <- "mlv2023"; m$parameters$min_branch_points <- 10L
m$branches <- trace_skeleton(m, 10L); y <- map_mat_properties(x, compute_mat_properties(m))
allpix <- readRDS(".agent-work/track-b/oct_forest/allpix.rds"); rc <- cbind(allpix$y, allpix$x)
for (p in c("parallelism", "separation", "mirror")) {
  r <- y$mat_images[[p]][rc]
  for (od in c("oct_forest", "oct_forest_single")) { o <- readBin(file.path(".agent-work/track-b", od, paste0("mapped_", p, ".bin")), "double", nrow(allpix)); o[o == 0] <- NA
    e <- abs(r - o); cat(sprintf("%-11s R compat vs Octave %-17s: both NA %d, NA mismatch %d, |e|>1e-9 %d, |e|>1e-4 %d of %d\n", p, od,
      sum(is.na(r) & is.na(o)), sum(is.na(r) != is.na(o)), sum(e > 1e-9, na.rm = TRUE), sum(e > 1e-4, na.rm = TRUE), nrow(rc))) } }
