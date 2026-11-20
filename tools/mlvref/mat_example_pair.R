# Track B: does examples/MAT.mat pair with examples/imgLD.mat? If so, it is an independent
# test bed for the ink-thinning and flux-thinning stages (aggregates only).
source("tools/mlvref/mat_common.R"); Rcpp::sourceCpp("tools/mlvref/mat_thin.cpp", cacheDir = ".agent-work/track-b/cpp-cache")
root <- mlvref_root()
M <- R.matlab::readMat(file.path(root, "reference-data/mlv/examples/MAT.mat"), fixNames = FALSE)
g <- function(f) M[[f]]   # top-level variables skeleton, distance_map, AOF
D <- g("distance_map"); A <- g("AOF"); S <- g("skeleton") != 0
I <- R.matlab::readMat(file.path(root, "reference-data/mlv/examples/imgLD.mat"))$imgLD
ink <- I[, , 1] == 0   # computeMAT: rgb2gray then ~binaryImage; mapMATtoContour uses imgLD == 0
cat("dims", dim(D), " ink pixels", sum(ink), " stored thin ink", sum(D <= 0), "\n")
fld <- mat_field(D <= 0); cat("stored D reproduced from D<=0:", sum(fround(fld$D) == D), "of", length(D), "\n")
cat("AOF from stored D + lowest-index labels: max |diff|", max(abs(average_outward_flux(D, fld$index) - A)), "\n")
mc <- max(floor(.025 * max(dim(D))), 1)
for (meth in c("pratt", "hmt8")) { t <- skeletonize(ink, method = meth)
  cat(sprintf("ink thinning %-5s: extra %d missing %d (MATLAB thin ink %d)\n", meth, sum(t & !(D <= 0)), sum((D <= 0) & !t), sum(D <= 0))) }
for (thr in c(.3, 2 / pi * sin(28 * pi / 180))) { pre <- A >= thr * 60
  for (meth in c("pratt", "hmt8")) { k <- skeletonize(pre, mc, method = meth)
    cat(sprintf("thr %.5f flux thinning %-5s: skel %d extra %d missing %d\n", thr, meth, sum(k), sum(k & !S), sum(S & !k))) }
  k <- ns$cpp_remove_small(thin_hmt_order(pre, as.integer(c(4,5,2,3,7,6,1,0)), 0L), as.integer(mc))
  cat(sprintf("thr %.5f flux thinning order 4,5,2,3,7,6,1,0: extra %d missing %d\n", thr, sum(k & !S), sum(S & !k))) }
# MAT.mat stores AOF as a binary mask containing its skeleton: test the thinning stage on it.
Ab <- A != 0
cat("binary AOF field:", all(A %in% c(0, 1)), " skeleton subset of it:", all(S <= Ab), "\n")
for (meth in c("pratt", "hmt8")) { k <- skeletonize(Ab, mc, method = meth); cat(sprintf("thin(AOF mask) %-5s: extra %d missing %d of %d\n", meth, sum(k & !S), sum(S & !k), sum(S))) }
k <- ns$cpp_remove_small(thin_hmt_order(Ab, as.integer(c(4,5,2,3,7,6,1,0)), 0L), as.integer(mc)); cat("thin(AOF mask) order 4,5,2,3,7,6,1,0: extra", sum(k & !S), "missing", sum(S & !k), "\n")
