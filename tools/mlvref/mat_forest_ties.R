# Track B: effect of nearest-boundary tie policy and AOF sampling on Forest_DBW (natural data,
# read locally; only aggregates are written).
source("tools/mlvref/mat_common.R"); load_tie_kernels()
f <- read_forest(mlvref_root()); ink <- f$D <= 0; fl <- mat_field(ink); ts <- field_ties(fl)
pol <- list(
  mlv_lowest = aof_matlab(f$D, tie_rule_labels(ts, dim(ink), "lowest"), 60L, 1),
  mlv_highest = aof_matlab(f$D, tie_rule_labels(ts, dim(ink), "highest"), 60L, 1),
  mlv_tie_mean = aof_tie_mean(f$D, ts$ptr, ts$lab, 60L, 1, FALSE, 1e-9),
  centered_tie_mean = aof_general(f$D, ts$ptr, ts$lab, 60L, 1, 0, 1, 1, 1e-9))
ref <- f$AOF; S0 <- ref >= 18
rows <- lapply(names(pol), function(p) { A <- pol[[p]]; e <- abs(A - ref); S <- A >= 18
  data.frame(policy = p, max_abs_vs_matlab = max(e), n_diff_1e9 = sum(e > 1e-9),
             thr_mask_n = sum(S), thr_mask_changed = sum(S != S0)) })
r <- do.call(rbind, rows); print(r, digits = 4)
write.csv(r, "docs/evidence/mlvref/mat/forest_tie_policy.csv", row.names = FALSE)
tsz <- table(pmin(diff(ts$ptr), 5)); cat("tie set sizes:", paste(names(tsz), tsz, sep = ":"), "\n")
# Downstream effect of the tie policy: hmt8 skeleton (MATLAB-exact on this image), then the
# R tracing/scoring/mapping pipeline; compare with stored MATLAB samples and with the lowest-label run.
x <- read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat", zero_length = "drop"); st <- stored_samples_forest(f$raw)
run <- function(A) { S <- skeletonize(A >= 18, 25L, method = "hmt8"); m <- make_mat(f$D, A, S, ink)
  m$branches <- trace_skeleton(m, 10L); m <- compute_mat_properties(m); list(S = S, img = map_mat_properties(x, m)$mat_images) }
base <- run(pol$mlv_lowest)
down <- do.call(rbind, lapply(names(pol), function(p) { r <- if (p == "mlv_lowest") base else run(pol[[p]])
  do.call(rbind, lapply(names(st), function(q) { cs <- compare_samples(r$img[[q]], st[[q]], 1e-4)
    a <- r$img[[q]]; b <- base$img[[q]]; both <- is.finite(a) & is.finite(b)
    data.frame(policy = p, property = q, skeleton_changed = sum(r$S != base$S), skeleton_vs_matlab = sum(r$S != f$skeleton),
      within_1e4_vs_matlab = cs[["within_tol"]], changed_vs_lowest = mean(abs(a[both] - b[both]) > 1e-4),
      median_abs_vs_lowest = median(abs(a[both] - b[both]))) })) }))
print(down, digits = 3, row.names = FALSE)
write.csv(down, "docs/evidence/mlvref/mat/forest_tie_downstream.csv", row.names = FALSE)
