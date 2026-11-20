# Track B: synthetic diagnostic figures (no natural data).
source("tools/mlvref/mat_common.R")
out <- "docs/evidence/mlvref/mat"
img <- function(M, main, col = c("white", "black")) { image(t(M[nrow(M):1, ]), col = col, axes = FALSE, main = main, asp = nrow(M) / ncol(M)); box() }
png(file.path(out, "thinning_synthetic.png"), 900, 330, res = 110)
par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))
m <- matrix(FALSE, 16, 30); m[3:12, 3:28] <- TRUE; m[13:14, 12:16] <- TRUE
img(m, "input (height-10 bar + stub)"); img(skeletonize(m, method = "pratt"), "pratt (Octave skel)"); img(skeletonize(m, method = "hmt8"), "hmt8")
dev.off()
png(file.path(out, "aof_registration_synthetic.png"), 900, 360, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
lines <- matrix(FALSE, 41, 41); lines[10, 4:38] <- lines[22, 4:38] <- TRUE
a <- compute_mat(lines, trace_branches = FALSE, aof = "mlv2023")$AOF[, 20]; c <- compute_mat(lines, trace_branches = FALSE, aof = "centered")$AOF[, 20]
thr <- 2 / pi * sin(28 * pi / 180) * 60
for (z in list(list(a, "mlv2023 AOF, column 20"), list(c, "centred tie-mean AOF, column 20"))) {
  plot(11:21, z[[1]][11:21], type = "h", lwd = 6, xlab = "row (ink rows 10 and 22; medial row 16)", ylab = "flux", main = z[[2]], ylim = c(-5, 40))
  abline(h = thr, lty = 2); abline(v = 16, col = "grey50") }
dev.off(); cat("figures written\n")
