# Independent R transcription of the pinned MLV junction pipeline
# (detectJunctions.m, lineIntersection.m, cleanupJunctions.m,
# computeJunctionAnglesTypes.m at f12418a). Each upstream heuristic or defect
# is a switch, so that a difference between mlvtoolbox and MATLAB can be
# attributed to one mechanism at a time:
#
#   snap      lineIntersection returns a segment endpoint when the
#             intersection parameter is within 1e-4 of 0 or 1.
#   negative  lineIntersection discards a computed (non-snapped) intersection
#             with any negative coordinate (MLV-M008).
#   short     detectJunctions skips contours of total length < AE.
#   merge     "mlv" = upper-triangular recursive classes (MLV-M006);
#             "transitive" = connected components.
#   arms      "mlv" = computeJunctionAnglesTypes (start-first attachment,
#             incident-segment orientation, one adjacent segment);
#             "exit" = mlvtoolbox definition: one arm per direction in which
#             each member's contour leaves the open disc of radius arm_thresh,
#             directed along the leaving segment, deduplicated per contour;
#             a segment wholly outside the disc and beyond the junction
#             leaves only away from it; closed contours wrap.
#   gap       "chebyshev" = AE/max(|dx|,|dy|) parameter threshold (upstream);
#             "euclidean" = AE/length, a rotation-invariant absolute gap.
#   angles    "mlv" = interior arms via atan2d(py-y1, px-x1) (MLV-M007);
#             "fixed" = interior arms in the image-y-down orientation
#             convention used by computeOrientation.m.
#
# The emulator works on the original MATLAB numbering, including zero-length
# segments (MLV-M015). It is harness code only and is validated against the
# stored and pinned-Octave outputs before being used for attribution.
mlv_emulate <- function(contours, AE = 1, RE = 0.3, merge_thresh = 2, arm_thresh = 2,
                        snap = TRUE, negative = TRUE, short = TRUE,
                        merge = c("mlv", "transitive"), angles = c("mlv", "fixed"),
                        gap = c("chebyshev", "euclidean"), arms = c("mlv", "exit")) {
  merge <- match.arg(merge); angles <- match.arg(angles); gap <- match.arg(gap); arms <- match.arg(arms)
  nc <- length(contours)
  sizes <- vapply(contours, nrow, 1L)
  Z <- do.call(rbind, contours)
  cid <- rep.int(seq_len(nc), sizes); sid <- sequence(sizes)
  ori_all <- (atan2(Z[, 2] - Z[, 4], Z[, 3] - Z[, 1]) * 180 / pi) %% 360
  len_all <- sqrt((Z[, 3] - Z[, 1])^2 + (Z[, 4] - Z[, 2])^2)
  clen <- vapply(split(len_all, factor(cid, levels = seq_len(nc))), sum, 0)
  eligible <- if (short) clen[cid] >= AE else rep(TRUE, nrow(Z))
  # lineIntersection, vectorized over reference segments for one query.
  Ay <- Z[, 3] - Z[, 1]; Ax <- Z[, 4] - Z[, 2] # (sic) upstream naming
  thr <- pmin(RE, AE / if (gap == "chebyshev") pmax(abs(Ax), abs(Ay)) else sqrt(Ax^2 + Ay^2))
  lo_x <- pmin(Z[, 1], Z[, 3]); hi_x <- pmax(Z[, 1], Z[, 3])
  lo_y <- pmin(Z[, 2], Z[, 4]); hi_y <- pmax(Z[, 2], Z[, 4])
  pad <- 2 * AE + 1
  raw <- vector("list", nrow(Z)); nr <- 0L
  for (q in which(eligible)) {
    r <- which(cid > cid[q] & eligible & lo_x <= hi_x[q] + pad & hi_x >= lo_x[q] - pad &
                 lo_y <= hi_y[q] + pad & hi_y >= lo_y[q] - pad)
    if (!length(r)) next
    By <- Ay[r]; Bx <- Ax[r]
    Cy <- Z[r, 1] - Z[q, 1]; Cx <- Z[r, 2] - Z[q, 2]
    D <- Ay[q] * Bx - Ax[q] * By
    a <- (Bx * Cy - By * Cx) / D; b <- (Ax[q] * Cy - Ay[q] * Cx) / D
    at <- thr[q]; bt <- thr[r]
    hit <- !is.na(a) & !is.na(b) & -at <= a & a <= 1 + at & -bt <= b & b <= 1 + bt
    if (!any(hit)) next
    r <- r[hit]; a <- a[hit]; b <- b[hit]
    px <- py <- numeric(length(r)); keep <- rep(TRUE, length(r))
    for (k in seq_along(r)) {
      rr <- r[k]
      if (snap && abs(a[k]) < 1e-4) { px[k] <- Z[q, 1]; py[k] <- Z[q, 2] }
      else if (snap && abs(a[k] - 1) < 1e-4) { px[k] <- Z[q, 3]; py[k] <- Z[q, 4] }
      else if (snap && abs(b[k]) < 1e-4) { px[k] <- Z[rr, 1]; py[k] <- Z[rr, 2] }
      else if (snap && abs(b[k] - 1) < 1e-4) { px[k] <- Z[rr, 3]; py[k] <- Z[rr, 4] }
      else {
        A1 <- Z[q, 2] - Z[q, 4]; B1 <- Z[q, 3] - Z[q, 1]; C1 <- Z[q, 1] * Z[q, 4] - Z[q, 2] * Z[q, 3]
        A2 <- Z[rr, 2] - Z[rr, 4]; B2 <- Z[rr, 3] - Z[rr, 1]; C2 <- Z[rr, 1] * Z[rr, 4] - Z[rr, 2] * Z[rr, 3]
        DD <- A1 * B2 - A2 * B1
        px[k] <- (B1 * C2 - B2 * C1) / DD; py[k] <- (A2 * C1 - A1 * C2) / DD
        if (negative && (px[k] < 0 || py[k] < 0)) keep[k] <- FALSE
      }
    }
    if (!any(keep)) next
    nr <- nr + 1L
    raw[[nr]] <- data.frame(q = q, r = r[keep], x = px[keep], y = py[keep])
  }
  raw <- do.call(rbind, raw[seq_len(nr)])
  if (is.null(raw)) return(empty_emulated())
  # detectJunctions loop order: queryC, queryS, refC, refS == global order of (q, r).
  raw <- raw[order(raw$q, raw$r), ]
  n <- nrow(raw)
  J <- lapply(seq_len(n), function(i) list(pos = c(raw$x[i], raw$y[i]),
    c = c(cid[raw$q[i]], cid[raw$r[i]]), s = c(sid[raw$q[i]], sid[raw$r[i]])))
  # Merge candidates (i < j, squared distance <= thresh^2).
  ox <- order(raw$x); edges <- list(); ne <- 0L
  for (t in seq_len(n)) {
    i <- ox[t]; u <- t + 1L
    while (u <= n && raw$x[ox[u]] - raw$x[i] <= merge_thresh) {
      j <- ox[u]
      if ((raw$x[i] - raw$x[j])^2 + (raw$y[i] - raw$y[j])^2 <= merge_thresh^2) {
        ne <- ne + 1L; edges[[ne]] <- c(min(i, j), max(i, j))
      }
      u <- u + 1L
    }
  }
  E <- if (ne) do.call(rbind, edges) else matrix(integer(), 0, 2)
  classes <- list()
  if (nrow(E)) {
    if (merge == "transitive") {
      g <- seq_len(n)
      find <- function(x) { while (g[x] != x) x <- g[x]; x }
      for (k in seq_len(nrow(E))) { a <- find(E[k, 1]); b <- find(E[k, 2]); if (a != b) g[max(a, b)] <- min(a, b) }
      root <- vapply(seq_len(n), find, 1L)
      in_class <- root %in% root[duplicated(root)]
      # emulate output ordering: classes ordered by first member's discovery
      for (rt in unique(root[in_class])) classes[[length(classes) + 1L]] <- which(root == rt)
      avail <- !in_class
    } else {
      active <- rep(TRUE, nrow(E)); avail <- rep(TRUE, n)
      rec <- function(Jx) {
        if (!avail[Jx]) return(integer())
        out <- Jx; avail[Jx] <<- FALSE
        row <- which(active & E[, 1] == Jx); newJs <- sort(E[row, 2]) # MATLAB find: ascending
        active[E[, 1] == Jx | E[, 2] == Jx] <<- FALSE
        for (v in newJs) out <- c(out, rec(v))
        out
      }
      repeat {
        if (!any(active)) break
        # first true of mergeMatrix in column-major order: smallest column, then row
        k <- which(active); k <- k[order(E[k, 2], E[k, 1])][1]
        classes[[length(classes) + 1L]] <- rec(E[k, 1])
      }
    }
  } else avail <- rep(TRUE, n)
  for (i in seq_len(n)) J[[i]]$idx <- i
  cleaned <- J[avail]
  for (cl in classes) {
    pos <- colMeans(do.call(rbind, lapply(J[cl], `[[`, "pos")))
    cc <- integer(); ss <- integer()
    for (j in cl) for (t in 1:2) {
      if (!any(cc == J[[j]]$c[t] & ss == J[[j]]$s[t])) { cc <- c(cc, J[[j]]$c[t]); ss <- c(ss, J[[j]]$s[t]) }
    }
    cleaned[[length(cleaned) + 1L]] <- list(pos = pos, c = cc, s = ss, idx = cl)
  }
  # computeJunctionAnglesTypes
  start <- c(0L, cumsum(sizes))
  gidx <- function(ci, si) start[ci] + si
  out <- list()
  # Neighbourhood: union of open discs around every merged intersection
  # (MLV-R014); see detect_junctions().
  exit_arms <- function(p, cs, ss, idx) {
    pts <- cbind(raw$x[idx], raw$y[idx]); r2 <- max(arm_thresh^2, .Machine$double.eps)
    inside <- function(g, end) { v <- Z[g, if (end) 3:4 else 1:2]; any((pts[, 1] - v[1])^2 + (pts[, 2] - v[2])^2 < r2) }
    keys <- character(); oris <- numeric()
    for (t in seq_along(cs)) for (dir in c(1L, -1L)) {
      ci <- cs[t]; j <- ss[t]; ok <- TRUE; n <- sizes[ci]
      g0 <- gidx(ci, j); a <- Z[g0, 1:2]; u <- Z[g0, 3:4] - a
      ownr <- idx[raw$q[idx] == g0 | raw$r[idx] == g0]
      tt <- ((raw$x[ownr] - a[1]) * u[1] + (raw$y[ownr] - a[2]) * u[2]) / sum(u^2)
      tc <- pmin(pmax(((pts[, 1] - a[1]) * u[1] + (pts[, 2] - a[2]) * u[2]) / sum(u^2), 0), 1)
      reaches <- any((a[1] + tc * u[1] - pts[, 1])^2 + (a[2] + tc * u[2] - pts[, 2])^2 < r2)
      fwd <- reaches || all(tt <= 0) || any(tt > 0 & tt < 1); bwd <- reaches || all(tt >= 1) || any(tt > 0 & tt < 1)
      if ((dir > 0 && !fwd) || (dir < 0 && !bwd)) next
      closed <- n > 1 && all(Z[gidx(ci, 1), 1:2] == Z[gidx(ci, n), 3:4])
      while (inside(gidx(ci, j), dir > 0)) {
        if (!closed && ((dir > 0 && j == n) || (dir < 0 && j == 1L))) { ok <- FALSE; break }
        j <- if (closed) (j - 1L + dir) %% n + 1L else j + dir
        if (j == ss[t]) { ok <- FALSE; break }
      }
      key <- paste(ci, j, dir)
      if (ok && !key %in% keys) {
        keys <- c(keys, key); g <- gidx(ci, j)
        oris <- c(oris, if (dir > 0) ori_all[g] else (ori_all[g] + 180) %% 360)
      }
    }
    oris
  }
  for (jj in cleaned) {
    p <- jj$pos; oris <- numeric(); cc <- jj$c; ss <- jj$s
    if (arms == "exit") oris <- exit_arms(p, jj$c, jj$s, jj$idx)
    for (t in seq_along(jj$c)) {
      if (arms == "exit") break
      ci <- jj$c[t]; si <- jj$s[t]; g <- gidx(ci, si)
      d1 <- sum((p - Z[g, 1:2])^2); d2 <- sum((p - Z[g, 3:4])^2)
      if (d1 < arm_thresh^2) {
        oris <- c(oris, ori_all[g])
        if (si > 1 && !any(cc == ci & ss == si - 1)) {
          cc <- c(cc, ci); ss <- c(ss, si - 1); oris <- c(oris, (ori_all[g - 1] + 180) %% 360)
        }
      } else if (d2 < arm_thresh^2) {
        oris <- c(oris, (ori_all[g] + 180) %% 360)
        if (si < sizes[ci] && !any(cc == ci & ss == si + 1)) {
          cc <- c(cc, ci); ss <- c(ss, si + 1); oris <- c(oris, ori_all[g + 1])
        }
      } else if (angles == "mlv") {
        oris <- c(oris, (atan2(p[2] - Z[g, 2], p[1] - Z[g, 1]) * 180 / pi) %% 360,
                  (atan2(p[2] - Z[g, 4], p[1] - Z[g, 3]) * 180 / pi) %% 360)
      } else {
        oris <- c(oris, (atan2(p[2] - Z[g, 2], Z[g, 1] - p[1]) * 180 / pi) %% 360,
                  (atan2(p[2] - Z[g, 4], Z[g, 3] - p[1]) * 180 / pi) %% 360)
      }
    }
    oris <- sort(oris)
    ang <- (c(oris[-1], oris[1]) - oris) %% 360
    nA <- length(ang)
    if (nA == 2 || (arms == "exit" && nA < 3)) next
    mx <- max(ang)
    type <- if (nA == 3) { if (mx < 160) "Y" else if (mx <= 200) "T" else "Arrow" } else if (nA == 4) "X" else "Star"
    out[[length(out) + 1L]] <- data.frame(x = p[1], y = p[2], type = type,
      contours = paste(cc, collapse = ";"), segments = paste(ss, collapse = ";"),
      angles = paste(sprintf("%.17g", ang), collapse = ";"), min_angle = min(ang), max_angle = mx)
  }
  if (!length(out)) empty_emulated() else do.call(rbind, out)
}
empty_emulated <- function() data.frame(x = double(), y = double(), type = character(),
  contours = character(), segments = character(), angles = character(),
  min_angle = double(), max_angle = double())
