# Track B: Forest_DBW end to end under the published-compatible chain and the new defaults.
source("tools/mlvref/mat_common.R")
f <- read_forest(mlvref_root()); st <- stored_samples_forest(f$raw)
x <- read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat", zero_length = "drop")
thr <- asin(.3 * pi / 2) * 180 / pi
run <- function(label, ...) { t0 <- Sys.time()
  m <- compute_mat(x, threshold_angle = thr, ...); m <- compute_mat_properties(m); y <- map_mat_properties(x, m)
  sec <- as.numeric(Sys.time() - t0, units = "secs")
  do.call(rbind, lapply(names(st), function(p) { r <- compare_samples(y$mat_images[[p]], st[[p]], 1e-4)
    data.frame(config = label, property = p, covered = r[["covered"]], within_1e4 = r[["within_tol"]],
               median_abs = r[["median_abs"]], max_abs = r[["max_abs"]], branches = length(m$branches), seconds = sec) })) }
r <- rbind(run("compat: aof=mlv2023, hmt8, branch_rule=mlv2023", aof = "mlv2023", morphology = "hmt8", branch_rule = "mlv2023"),
           run("defaults: centered, minimal, invariant"),
           run("previous default: centered, hmt8, invariant", morphology = "hmt8"),
           run("mixed: aof=mlv2023, hmt8, invariant", aof = "mlv2023", morphology = "hmt8"))
print(r, digits = 4, row.names = FALSE)
write.csv(r, "docs/evidence/mlvref/mat/forest_compat.csv", row.names = FALSE)
