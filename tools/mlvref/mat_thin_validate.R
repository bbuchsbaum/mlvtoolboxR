# Track B: validate candidate MATLAB-skel orders on stages not used to find them:
# (1) Forest ink thinning: R raster -> MATLAB thin ink (D <= 0); (2) imgLD -> MAT.mat thin ink.
source("tools/mlvref/mat_common.R"); Rcpp::sourceCpp("tools/mlvref/mat_thin.cpp", cacheDir = ".agent-work/track-b/cpp-cache")
f <- read_forest(mlvref_root()); x <- read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat", zero_length = "drop")
R <- rasterize(x); T1 <- f$D <= 0
M <- R.matlab::readMat("reference-data/mlv/examples/MAT.mat", fixNames = FALSE); I <- R.matlab::readMat("reference-data/mlv/examples/imgLD.mat")$imgLD
ink2 <- I[, , 1] == 0; T2 <- M$distance_map <= 0
pre <- f$AOF >= 18; S <- f$skeleton
orders <- list(clockwise_from_E = c(2:7, 0:1), S_SW_E_SE_W_NW_NE_N = c(4,5,2,3,6,7,1,0), S_SW_E_SE_W_NE_N_NW = c(4,5,2,3,6,1,0,7), S_SW_E_SE_W_NE_NW_N = c(4,5,2,3,6,1,7,0))
for (nm in names(orders)) { o <- as.integer(orders[[nm]])
  a <- thin_hmt_order(R, o, 0L); b <- thin_hmt_order(ink2, o, 0L); k <- ns$cpp_remove_small(thin_hmt_order(pre, o, 0L), 25L)
  cat(sprintf("%-22s flux skel diff %4d | Forest ink thin diff %3d | imgLD ink thin diff %3d\n", nm, sum(k != S), sum(a != T1), sum(b != T2))) }
