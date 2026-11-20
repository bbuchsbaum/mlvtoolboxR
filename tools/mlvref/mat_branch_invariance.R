# Track B: sensitivity of mlv2023 branch scores to traversal direction and to the duplicated
# merge points emitted by edgelink/cleanedgelist (Forest_DBW, MATLAB skeleton, MLV mapping order).
source("tools/mlvref/mat_common.R")
f <- read_forest(mlvref_root()); ink <- f$D <= 0; st <- stored_samples_forest(f$raw)
od <- ".agent-work/track-b/oct_forest_single"; allpix <- readRDS(file.path(od, "allpix.rds")); rc <- cbind(allpix$y, allpix$x)
b <- read.csv(file.path(od, "branches.csv"), header = FALSE, col.names = c("b", "row", "col", "radius"))
base <- unname(lapply(split(b, b$b), function(z) data.frame(row = z$row, column = z$col, x = z$col, y = z$row, radius = f$D[cbind(z$row, z$col)], AOF = 0)))
variants <- list(edgelink = base,
  reversed = lapply(base, function(z) z[nrow(z):1, ]),
  deduplicated = lapply(base, function(z) z[c(TRUE, diff(z$row) != 0 | diff(z$column) != 0), ]))
rows <- list()
for (nm in names(variants)) { m <- make_mat(f$D, f$AOF, f$skeleton, ink); m$branches <- variants[[nm]]; m <- compute_mat_properties(m)
  v <- map_at(m, rc, matlab_image = TRUE)
  for (p in names(st)) { im <- matrix(NA_real_, 768, 1024); im[rc] <- v[[p]]; r <- compare_samples(im, st[[p]], 1e-4)
    rows[[length(rows) + 1]] <- data.frame(branches = nm, property = p, within_1e4 = r[["within_tol"]], median_abs = r[["median_abs"]], q90_abs = r[["q90_abs"]]) } }
r <- do.call(rbind, rows); print(r, digits = 3, row.names = FALSE)
write.csv(r, "docs/evidence/mlvref/mat/forest_branch_invariance.csv", row.names = FALSE)
