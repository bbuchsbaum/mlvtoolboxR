# Track B: D4 equivariance of AOF tie policies on synthetic thin-ink masks.
# Usage: Rscript tools/mlvref/mat_equivariance.R  (from repository root; MLVREF_LIB optional)
source("tools/mlvref/mat_common.R"); load_tie_kernels()
synth <- function() {
  m <- list()
  a <- matrix(FALSE, 41, 41); a[8, 8:34] <- a[34, 8:34] <- a[8:34, 8] <- a[8:34, 34] <- TRUE; m$square <- a
  b <- matrix(FALSE, 41, 41); th <- seq(0, 2*pi, length.out = 2000)
  b[cbind(round(21 + 13*cos(th)), round(21 + 13*sin(th)))] <- TRUE; m$circle <- b
  p <- matrix(FALSE, 41, 41); p[10, 6:36] <- p[22, 6:36] <- TRUE; m$parallel_even <- p   # gap 11 (odd interior)
  q <- matrix(FALSE, 41, 41); q[10, 6:36] <- q[21, 6:36] <- TRUE; m$parallel_odd <- q    # gap 10 (even interior)
  r <- matrix(FALSE, 41, 41); r[cbind(6:34, 6:34)] <- TRUE; r[cbind(6:34, 40 - (6:34))] <- TRUE; m$cross_diag <- r
  s <- matrix(FALSE, 41, 41); s[6:34, 6] <- TRUE; s[34, 6:34] <- TRUE; s[cbind(6:34, 6:34)] <- TRUE; m$triangle <- s
  m
}
# MLV's AOF at pixel (i,j) is the flux through a circle centred on the pixel corner
# (i+0.5, j+0.5) (computeAOF.m adds 0.5 to the sample points). Equivariance is therefore
# tested on that dual lattice: a transform g maps corner P to g(P), i.e. pixel g(P)-0.5.
D4 <- list(
  rot90 = list(f = function(x) t(x[nrow(x):1, ]), p = function(r, c, n) cbind(c, n + 1 - r)),
  flip_v = list(f = function(x) x[nrow(x):1, ], p = function(r, c, n) cbind(n + 1 - r, c)),
  flip_h = list(f = function(x) x[, ncol(x):1], p = function(r, c, n) cbind(r, n + 1 - c)),
  transpose = list(f = function(x) t(x), p = function(r, c, n) cbind(c, r)))
policies <- list(
  lowest = function(fl, ts) aof_matlab(fl$D, tie_rule_labels(ts, dim(fl$D), "lowest"), 60L, 1),
  highest = function(fl, ts) aof_matlab(fl$D, tie_rule_labels(ts, dim(fl$D), "highest"), 60L, 1),
  tie_mean = function(fl, ts) aof_tie_mean(fl$D, ts$ptr, ts$lab, 60L, 1, FALSE, 1e-9),
  tie_mean_sym = function(fl, ts) aof_tie_mean(fl$D, ts$ptr, ts$lab, 60L, 1, TRUE, 1e-9),
  centered_lowest = function(fl, ts) aof_general(fl$D, ts$ptr, ts$lab, 60L, 1, 0, 1, 0, 1e-9),
  centered_tie_mean = function(fl, ts) aof_general(fl$D, ts$ptr, ts$lab, 60L, 1, 0, 1, 1, 1e-9))
lattice_offset <- c(lowest = .5, highest = .5, tie_mean = .5, tie_mean_sym = .5, centered_lowest = 0, centered_tie_mean = 0)
aof_of <- function(mask, pol) { fl <- mat_field(mask); policies[[pol]](fl, field_ties(fl)) }
run_equivariance <- function() {
  rows <- list()
  for (nm in names(synth())) { mask <- synth()[[nm]]
    for (pol in names(policies)) { A <- aof_of(mask, pol)
      for (g in names(D4)) { G <- aof_of(D4[[g]]$f(mask), pol); n <- nrow(A)
        # interior corners only: the epsilon border exclusion (i > eps, i < m - eps) is itself asymmetric
        ij <- as.matrix(expand.grid(3:(n - 3), 3:(n - 3)))
        o <- lattice_offset[[pol]]; q <- D4[[g]]$p(ij[, 1] + o, ij[, 2] + o, n) - o
        a <- A[ij]; b <- G[q]; e <- abs(a - b)
        rows[[length(rows) + 1]] <- data.frame(mask = nm, policy = pol, transform = g, max_abs = max(e),
          n_bad = sum(e > 1e-9), skel_diff = sum((a >= 18) != (b >= 18))) } } }
  do.call(rbind, rows)
}
if (sys.nframe() == 0L) {
  r <- run_equivariance()
  print(aggregate(cbind(max_abs, n_bad, skel_diff) ~ policy, r, max))
  out <- "docs/evidence/mlvref/mat/equivariance.csv"; write.csv(r, out, row.names = FALSE); cat("wrote", out, "\n")
}
