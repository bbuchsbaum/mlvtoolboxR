# Track B: compare pinned-Octave chain output (on MATLAB's stored MAT fields) with stored per-pixel scores.
source("tools/mlvref/mat_common.R")
d <- Sys.getenv("OD", ".agent-work/track-b/oct_forest_single"); f <- read_forest(mlvref_root()); st <- stored_samples_forest(f$raw)
allpix <- readRDS(file.path(d, "allpix.rds"))
res <- list()
for (p in names(st)) { fn <- file.path(d, paste0("mapped_", p, ".bin")); if (!file.exists(fn)) next
  v <- readBin(fn, "double", nrow(allpix)); img <- matrix(NA_real_, 768, 1024); img[cbind(allpix$y, allpix$x)] <- ifelse(v == 0, NA, v)
  res[[p]] <- compare_samples(img, st[[p]], 1e-12) }
print(do.call(rbind, res), digits = 4)
