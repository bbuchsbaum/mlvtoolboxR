# Track B: identify the rasterization rule behind MATLAB insertShape('Line', LineWidth 1,
# SmoothEdges false) using a proxy of MATLAB's raster (thin ink from distance_map plus every
# stored per-pixel score position). Natural data stay local; only counts are printed.
source("tools/mlvref/mat_common.R")
rast_rule <- function(z, h, w, rule) {
  pts <- list()
  for (i in seq_len(nrow(z))) {
    x0 <- z[i, 1]; y0 <- z[i, 2]; x1 <- z[i, 3]; y1 <- z[i, 4]
    if (rule == "dda_round") {            # sample the major axis at integers, round the minor
      n <- max(abs(round(x1) - round(x0)), abs(round(y1) - round(y0)))
      t <- if (n == 0) 0 else (0:n) / n
      xs <- x0 + t * (x1 - x0); ys <- y0 + t * (y1 - y0)
      pts[[i]] <- cbind(round(xs), round(ys))
    } else if (rule == "dda_major") {       # step major axis over rounded endpoints, interpolate exact line
      ax <- round(x0); bx <- round(x1); ay <- round(y0); by <- round(y1)
      if (abs(bx - ax) >= abs(by - ay)) { xs <- if (ax == bx) ax else ax:bx
        ys <- if (x1 == x0) rep(ay, length(xs)) else round(y0 + (xs - x0) * (y1 - y0) / (x1 - x0)); pts[[i]] <- cbind(xs, ys)
      } else { ys <- ay:by; xs <- if (y1 == y0) rep(ax, length(ys)) else round(x0 + (ys - y0) * (x1 - x0) / (y1 - y0)); pts[[i]] <- cbind(xs, ys) }
    } else if (rule == "dda_rounded_endpoints") { # integer DDA between rounded endpoints
      ax <- round(x0); bx <- round(x1); ay <- round(y0); by <- round(y1); n <- max(abs(bx - ax), abs(by - ay))
      t <- if (n == 0) 0 else (0:n) / n; pts[[i]] <- cbind(round(ax + t * (bx - ax)), round(ay + t * (by - ay)))
    } else if (rule == "slope_half") {     # major-axis DDA; exact half ties go up for positive slope, down for negative
      dx <- x1 - x0; dy <- y1 - y0; up <- sign(dx) == sign(dy)
      rnd <- function(v) if (up) floor(v + .5) else ceiling(v - .5)
      if (abs(dx) >= abs(dy)) { xs <- if (dx == 0) x0 else seq(x0, x1); ys <- if (dx == 0) y0 else rnd(y0 + (xs - x0) * dy / dx); pts[[i]] <- cbind(xs, ys) }
      else { ys <- seq(y0, y1); xs <- rnd(x0 + (ys - y0) * dx / dy); pts[[i]] <- cbind(xs, ys) }
    } else if (rule == "floor_half") {     # DDA with floor(x + 0.5) (MATLAB round differs at negative halves only)
      n <- max(abs(floor(x1 + .5) - floor(x0 + .5)), abs(floor(y1 + .5) - floor(y0 + .5)))
      t <- if (n == 0) 0 else (0:n) / n
      pts[[i]] <- cbind(floor(x0 + t * (x1 - x0) + .5), floor(y0 + t * (y1 - y0) + .5))
    }
  }
  p <- do.call(rbind, pts); p <- p[p[, 1] >= 1 & p[, 1] <= w & p[, 2] >= 1 & p[, 2] <= h, , drop = FALSE]
  m <- matrix(FALSE, h, w); m[cbind(p[, 2], p[, 1])] <- TRUE; m
}
if (sys.nframe() == 0L) {
  f <- read_forest(mlvref_root()); st <- stored_samples_forest(f$raw)
  proxy <- f$D <= 0; for (s in st) proxy[cbind(s$y, s$x)] <- TRUE
  x <- read_mlv_mat("reference-data/mlv/examples/Forest_DBW.mat", zero_length = "drop")
  z <- do.call(rbind, x$contours); h <- x$height; w <- x$width
  cat("segments", nrow(z), " integer endpoints:", mean(z == round(z)), " proxy pixels", sum(proxy), "\n")
  res <- list(r_package = rasterize(x))
  for (rule in c("dda_round", "dda_major", "dda_rounded_endpoints", "floor_half", "slope_half")) res[[rule]] <- rast_rule(z, h, w, rule)
  tab <- do.call(rbind, lapply(names(res), function(n) data.frame(rule = n, pixels = sum(res[[n]]),
    proxy_missing = sum(proxy & !res[[n]]), not_in_proxy = sum(res[[n]] & !proxy))))
  print(tab); write.csv(tab, "docs/evidence/mlvref/mat/forest_raster_rules.csv", row.names = FALSE)
}
# The pre-qualification package rasterizer (Bresenham with error-term ties) for reference:
# 23,754 pixels, 427 proxy pixels missing, 430 outside the proxy; reversing all segments
# changed 1,984 pixels. Recorded from the 0.2.1 build before MLV-R010 was fixed.
