# MLV-compatible junction rule (detect_junctions(rule = "mlv2023")).
# A port of detectJunctions.m, lineIntersection.m, cleanupJunctions.m and
# computeJunctionAnglesTypes.m at MLV f12418a, including their documented
# defects (MLV-M006, M007, M008, M019, M021-M026; see UPSTREAM_BUGS.md).
# It exists to reproduce published MLV junction values; the default rule
# ("exit") is the corrected definition. Zero-length segments cannot occur in
# an mlv_drawing, so MLV-M020 has no counterpart here.
.junctions_mlv2023 <- function(x, AE, RE, merge_radius, endpoint_radius) {
  contours <- x$contours; nc <- length(contours)
  if (!nc) return(.empty_junctions())
  sizes <- vapply(contours, nrow, integer(1))
  Z <- unname(do.call(rbind, contours))
  cid <- rep.int(seq_len(nc), sizes); sid <- sequence(sizes)
  ori <- (atan2(Z[, 2] - Z[, 4], Z[, 3] - Z[, 1]) * 180 / pi) %% 360
  len <- sqrt((Z[, 3] - Z[, 1])^2 + (Z[, 4] - Z[, 2])^2)
  eligible <- (vapply(split(len, factor(cid, levels = seq_len(nc))), sum, double(1)) >= AE)[cid]
  # lineIntersection.m, vectorized over reference segments of later contours.
  Ay <- Z[, 3] - Z[, 1]; Ax <- Z[, 4] - Z[, 2] # upstream naming
  thr <- pmin(RE, AE / pmax(abs(Ax), abs(Ay)))
  lo_x <- pmin(Z[, 1], Z[, 3]); hi_x <- pmax(Z[, 1], Z[, 3])
  lo_y <- pmin(Z[, 2], Z[, 4]); hi_y <- pmax(Z[, 2], Z[, 4])
  pad <- 2 * AE + 1 # bounds any parameter extension in each coordinate
  raw <- vector("list", nrow(Z)); nr <- 0L
  for (q in which(eligible)) {
    r <- which(cid > cid[q] & eligible & lo_x <= hi_x[q] + pad & hi_x >= lo_x[q] - pad &
               lo_y <= hi_y[q] + pad & hi_y >= lo_y[q] - pad)
    if (!length(r)) next
    Cy <- Z[r, 1] - Z[q, 1]; Cx <- Z[r, 2] - Z[q, 2]
    D <- Ay[q] * Ax[r] - Ax[q] * Ay[r]
    a <- (Ax[r] * Cy - Ay[r] * Cx) / D; b <- (Ax[q] * Cy - Ay[q] * Cx) / D
    hit <- !is.na(a) & !is.na(b) & -thr[q] <= a & a <= 1 + thr[q] & -thr[r] <= b & b <= 1 + thr[r]
    if (!any(hit)) next
    r <- r[hit]; a <- a[hit]; b <- b[hit]
    s1 <- abs(a) < 1e-4; s2 <- !s1 & abs(a - 1) < 1e-4
    s3 <- !s1 & !s2 & abs(b) < 1e-4; s4 <- !s1 & !s2 & !s3 & abs(b - 1) < 1e-4
    free <- !(s1 | s2 | s3 | s4)
    A1 <- Z[q, 2] - Z[q, 4]; B1 <- Z[q, 3] - Z[q, 1]; C1 <- Z[q, 1] * Z[q, 4] - Z[q, 2] * Z[q, 3]
    A2 <- Z[r, 2] - Z[r, 4]; B2 <- Z[r, 3] - Z[r, 1]; C2 <- Z[r, 1] * Z[r, 4] - Z[r, 2] * Z[r, 3]
    DD <- A1 * B2 - A2 * B1
    px <- (B1 * C2 - B2 * C1) / DD; py <- (A2 * C1 - A1 * C2) / DD
    px[s1] <- Z[q, 1]; py[s1] <- Z[q, 2]; px[s2] <- Z[q, 3]; py[s2] <- Z[q, 4]
    px[s3] <- Z[r[s3], 1]; py[s3] <- Z[r[s3], 2]; px[s4] <- Z[r[s4], 3]; py[s4] <- Z[r[s4], 4]
    keep <- !(free & (px < 0 | py < 0)) # MLV-M008
    if (!any(keep)) next
    nr <- nr + 1L; raw[[nr]] <- cbind(q, r[keep], px[keep], py[keep])
  }
  if (!nr) return(.empty_junctions())
  raw <- do.call(rbind, raw[seq_len(nr)])
  raw <- raw[order(raw[, 1], raw[, 2]), , drop = FALSE] # queryC, queryS, refC, refS order
  n <- nrow(raw); P <- unname(raw[, 3:4, drop = FALSE])
  # cleanupJunctions.m: upper-triangular recursive classes (MLV-M006).
  ox <- order(P[, 1]); E <- list(); ne <- 0L
  for (t in seq_len(n)) {
    i <- ox[t]; u <- t + 1L
    while (u <= n && P[ox[u], 1] - P[i, 1] <= merge_radius) {
      j <- ox[u]
      if (sum((P[i, ] - P[j, ])^2) <= merge_radius^2) { ne <- ne + 1L; E[[ne]] <- c(min(i, j), max(i, j)) }
      u <- u + 1L
    }
  }
  avail <- rep(TRUE, n); classes <- list()
  if (ne) {
    E <- do.call(rbind, E); active <- rep(TRUE, nrow(E))
    # recursiveClasses, iteratively: pre-order depth-first visit in which each
    # junction's merge partners are taken in ascending index order (MATLAB
    # find on the matrix row), read when the junction is visited, after which
    # its row and column are cleared.
    rec <- function(J0) {
      if (!avail[J0]) return(integer())
      out <- integer(); stack <- list()
      visit <- function(J) {
        avail[J] <<- FALSE; out <<- c(out, J)
        kids <- sort(E[active & E[, 1] == J, 2])
        active[E[, 1] == J | E[, 2] == J] <<- FALSE
        stack[[length(stack) + 1L]] <<- list(kids = kids, pos = 1L)
      }
      visit(J0)
      while (length(stack)) {
        top <- length(stack); fr <- stack[[top]]
        if (fr$pos > length(fr$kids)) { stack[[top]] <- NULL; next }
        v <- fr$kids[fr$pos]; stack[[top]]$pos <- fr$pos + 1L
        if (avail[v]) visit(v)
      }
      out
    }
    while (any(active)) {
      k <- which(active); k <- k[order(E[k, 2], E[k, 1])][1] # first true, column-major
      classes[[length(classes) + 1L]] <- rec(E[k, 1])
    }
  }
  jl <- lapply(which(avail), function(i) list(p = P[i, ], c = cid[raw[i, 1:2]], s = sid[raw[i, 1:2]]))
  for (cl in classes) {
    cc <- integer(); ss <- integer()
    for (j in cl) for (g in raw[j, 1:2]) if (!any(cc == cid[g] & ss == sid[g])) { cc <- c(cc, cid[g]); ss <- c(ss, sid[g]) }
    jl[[length(jl) + 1L]] <- list(p = colMeans(P[cl, , drop = FALSE]), c = cc, s = ss)
  }
  # computeJunctionAnglesTypes.m (MLV-M007, M021, M022, M024, M025, M026).
  start <- c(0L, cumsum(sizes)); r2 <- endpoint_radius^2
  result <- list()
  for (jj in jl) {
    p <- jj$p; cc <- jj$c; ss <- jj$s; o <- numeric()
    for (t in seq_along(jj$c)) {
      ci <- jj$c[t]; si <- jj$s[t]; g <- start[ci] + si
      if (sum((p - Z[g, 1:2])^2) < r2) {
        o <- c(o, ori[g])
        if (si > 1L && !any(cc == ci & ss == si - 1L)) { cc <- c(cc, ci); ss <- c(ss, si - 1L); o <- c(o, (ori[g - 1L] + 180) %% 360) }
      } else if (sum((p - Z[g, 3:4])^2) < r2) {
        o <- c(o, (ori[g] + 180) %% 360)
        if (si < sizes[ci] && !any(cc == ci & ss == si + 1L)) { cc <- c(cc, ci); ss <- c(ss, si + 1L); o <- c(o, ori[g + 1L]) }
      } else {
        o <- c(o, (atan2(p[2] - Z[g, 2], p[1] - Z[g, 1]) * 180 / pi) %% 360,
               (atan2(p[2] - Z[g, 4], p[1] - Z[g, 3]) * 180 / pi) %% 360)
      }
    }
    o <- sort(unname(o)); angles <- (c(o[-1L], o[1L]) - o) %% 360; k <- length(angles)
    if (k == 2L) next
    mx <- max(angles)
    type <- if (k == 3L) { if (mx < 160) "Y" else if (mx <= 200) "T" else "Arrow" } else if (k == 4L) "X" else "Star"
    z <- data.frame(junction_id = length(result) + 1L, x = unname(p[1]), y = unname(p[2]), type = type,
                    min_angle = min(angles), max_angle = mx)
    z$angles <- I(list(angles)); z$arms <- I(list(o))
    z$members <- I(list(data.frame(contour_index = cc, segment_index = ss, contour_id = x$ids[cc])))
    result[[length(result) + 1L]] <- z
  }
  if (!length(result)) .empty_junctions() else do.call(rbind, result)
}
