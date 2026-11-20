# Line extraction (extract_lines, xdog_lines, line_extraction_stages).
# Oracles are independent of the implementation: analytic shapes rendered
# with 8 x 8 supersampled anti-aliasing and known boundaries, lattice
# symmetries, and documented contracts.

lx_render <- function(h, w, inside, fg = 0.2, bg = 0.9, ss = 8) {
 o <- (seq_len(ss) - 0.5) / ss - 0.5
 cover <- matrix(0, h, w)
 for (a in o) for (b in o) {
  xy <- expand.grid(y = seq_len(h) + b, x = seq_len(w) + a)
  cover <- cover + matrix(inside(xy$x, xy$y), h, w)
 }
 bg + (fg - bg) * cover / ss^2
}
lx_vertices <- function(d) {
 z <- do.call(rbind, d$contours)
 rbind(z[, 1:2, drop = FALSE], z[, 3:4, drop = FALSE])
}
lx_lengths <- function(d) vapply(d$contours, function(z)
 sum(sqrt(rowSums((z[, 3:4, drop = FALSE] - z[, 1:2, drop = FALSE])^2))), 0)
lx_closed <- function(z) all(z[1, 1:2] == z[nrow(z), 3:4])
# Distance to the boundary of a square of half-side `half` centred at (cx, cy)
# and rotated by `angle`.
lx_square_distance <- function(x, y, cx, cy, half, angle) {
 u <- (x - cx) * cos(angle) + (y - cy) * sin(angle)
 v <- -(x - cx) * sin(angle) + (y - cy) * cos(angle)
 du <- abs(u) - half; dv <- abs(v) - half
 out <- sqrt(pmax(du, 0)^2 + pmax(dv, 0)^2)
 inside <- du <= 0 & dv <= 0
 out[inside] <- pmin(-du, -dv)[inside]
 out
}
# Short strokes of random orientation and lightness offset, like brush texture.
lx_brush <- function(h, w, n, len = 7, amp = 0.2, seed = 1) {
 set.seed(seed); out <- matrix(0, h, w)
 for (i in seq_len(n)) {
  x0 <- runif(1, 1, w); y0 <- runif(1, 1, h); th <- runif(1, 0, pi); a <- runif(1, -amp, amp)
  t <- seq(-len / 2, len / 2, by = 0.5)
  xs <- round(x0 + t * cos(th)); ys <- round(y0 + t * sin(th))
  ok <- xs >= 1 & xs <= w & ys >= 1 & ys <= h
  out[cbind(ys[ok], xs[ok])] <- a
 }
 out
}
lx_plain <- function(m) { attributes(m) <- list(dim = dim(m)); m }
lx_rot90 <- function(a) if (length(dim(a)) == 3L) aperm(a[dim(a)[1]:1, , , drop = FALSE], c(2, 1, 3)) else t(a[nrow(a):1, , drop = FALSE])

test_that("a disc is traced as one closed contour on its boundary", {
 for (fg in c(0.2, 0.9)) {
  img <- lx_render(96, 96, function(x, y) (x - 48.3)^2 + (y - 47.6)^2 < 30^2, fg = fg, bg = 1.1 - fg)
  ink <- extract_lines(img)
  d <- trace_contours(ink)
  expect_length(d$contours, 1L)
  expect_true(lx_closed(d$contours[[1]]))
  v <- lx_vertices(d)
  expect_lte(max(abs(sqrt((v[, 1] - 48.3)^2 + (v[, 2] - 47.6)^2) - 30)), 1)
  s <- trace_contours(ink, simplify = 1)
  expect_length(s$contours, 1L)
  expect_lt(abs(lx_lengths(s) / (2 * pi * 30) - 1), 0.02)
 }
})

test_that("squares at several orientations are traced within 1 px of their sides", {
 square <- function(cx, cy, half, angle) lx_render(96, 96, function(x, y) {
  u <- (x - cx) * cos(angle) + (y - cy) * sin(angle)
  v <- -(x - cx) * sin(angle) + (y - cy) * cos(angle)
  abs(u) < half & abs(v) < half
 })
 for (angle in seq(0, pi / 4, length.out = 10)) {
  img <- square(48.2, 47.7, 25, angle)
  # Without the pre-filter every vertex is within 1 px of the square.
  d <- trace_contours(extract_lines(img, smooth = 0))
  expect_length(d$contours, 1L)
  expect_true(lx_closed(d$contours[[1]]))
  v <- lx_vertices(d)
  expect_lte(max(lx_square_distance(v[, 1], v[, 2], 48.2, 47.7, 25, angle)), 1)
  # The default Kuwahara pre-filter rounds the corners slightly; away from
  # the corners the sides stay within 1 px.
  d <- trace_contours(extract_lines(img))
  expect_length(d$contours, 1L)
  expect_true(lx_closed(d$contours[[1]]))
  v <- lx_vertices(d)
  dist <- lx_square_distance(v[, 1], v[, 2], 48.2, 47.7, 25, angle)
  u <- (v[, 1] - 48.2) * cos(angle) + (v[, 2] - 47.7) * sin(angle)
  w <- -(v[, 1] - 48.2) * sin(angle) + (v[, 2] - 47.7) * cos(angle)
  corner <- abs(abs(u) - 25) < 4 & abs(abs(w) - 25) < 4
  expect_lte(max(dist[!corner]), 1)
 }
 # Random squares (half-side 15 to 30 px, any angle, sub-pixel centre): the
 # documented corner bound is about 1.2 px with the pre-filter, 1 px without.
 set.seed(2026)
 for (i in 1:20) {
  half <- runif(1, 15, 30); angle <- runif(1, 0, pi / 2)
  cx <- 48 + runif(1, -0.5, 0.5); cy <- 48 + runif(1, -0.5, 0.5)
  img <- square(cx, cy, half, angle)
  v <- lx_vertices(trace_contours(extract_lines(img)))
  expect_lte(max(lx_square_distance(v[, 1], v[, 2], cx, cy, half, angle)), 1.25)
  v <- lx_vertices(trace_contours(extract_lines(img, smooth = 0)))
  expect_lte(max(lx_square_distance(v[, 1], v[, 2], cx, cy, half, angle)), 1)
 }
})

test_that("step and narrow ramp edges are traced within 1 px of the edge", {
 img <- lx_render(64, 64, function(x, y) x * cos(0.3) + y * sin(0.3) < 30.3)
 d <- trace_contours(extract_lines(img))
 expect_length(d$contours, 1L)
 v <- lx_vertices(d)
 expect_lte(max(abs(v[, 1] * cos(0.3) + v[, 2] * sin(0.3) - 30.3)), 1)
 # The line spans the image: its end points are within 2 px of the borders.
 expect_gt(diff(range(v[, 2])), 60)
 for (width in c(3, 5)) {
  ramp <- outer(1:64, 1:64, function(y, x) 0.2 + 0.7 * pmin(1, pmax(0, (x - 32.3) / width + 0.5)))
  v <- lx_vertices(trace_contours(extract_lines(ramp)))
  expect_lte(max(abs(v[, 1] - 32.3)), 1)
 }
})

test_that("the pre-filter suppresses brush texture that plain flow-DoG draws", {
 step <- lx_render(96, 96, function(x, y) x * cos(0.25) + y * sin(0.25) < 45.5)
 far_ink <- function(m) {
  ij <- which(m, arr.ind = TRUE)
  sum(abs(ij[, 2] * cos(0.25) + ij[, 1] * sin(0.25) - 45.5) > 3)
 }
 # Twenty textures. Measured over seeds 1 to 30 (fdog-xdog-1): the
 # pre-filtered mask has less off-edge ink than the plain one in every seed
 # and 8.4 times less in total. Over seeds 1 to 20: medians 21.5 px
 # (pre-filtered) and 151.5 px (plain), 8 seeds with none; thresholds below
 # leave about a factor of two of margin.
 res <- t(vapply(1:20, function(seed) {
  img <- pmin(pmax(step + lx_brush(96, 96, 500, seed = seed), 0), 1)
  smoothed <- extract_lines(img)
  d <- trace_contours(smoothed); z <- d$contours[[which.max(lx_lengths(d))]]
  c(smoothed = far_ink(smoothed), plain = far_ink(extract_lines(img, smooth = 0)),
    deviation = max(abs(z[, 1] * cos(0.25) + z[, 2] * sin(0.25) - 45.5)), longest = max(lx_lengths(d)))
 }, numeric(4)))
 expect_true(all(res[, "smoothed"] < res[, "plain"]))
 expect_gte(sum(res[, "plain"]), 4 * sum(res[, "smoothed"]))
 expect_lte(median(res[, "smoothed"]), 40)
 expect_gte(median(res[, "plain"]), 75)
 # The longest traced contour follows the step across the whole image.
 expect_true(all(res[, "deviation"] <= 1))
 expect_true(all(res[, "longest"] > 95))
})

test_that("integration along the flow joins a dashed stroke", {
 # Dashes 6 px long, 3 px apart, 2 px thick, from x = 9 to 88.
 img <- outer(1:48, 1:96, function(y, x)
  ifelse(y %in% 24:25 & ((x - 1) %% 9) < 6 & x > 8 & x < 89, 0.5, 0.85))
 joined <- trace_contours(extract_lines(img, sigma_m = 3, smooth = 0, min_component = 1))
 expect_length(joined$contours, 1L)
 v <- lx_vertices(joined)
 expect_gte(diff(range(v[, 1])), 75)
 expect_true(all(v[, 2] >= 23 & v[, 2] <= 26))
 broken <- trace_contours(extract_lines(img, sigma_m = 0, smooth = 0, min_component = 1))
 expect_gte(length(broken$contours), 8L)
 expect_lt(max(lx_lengths(broken)), 10)
})

test_that("thin dark strokes give a centre line and wide bars give two edges", {
 for (width in 1:4) {
  img <- matrix(0.85, 48, 64); img[24:(23 + width), 8:56] <- 0.3
  d <- trace_contours(extract_lines(img))
  expect_length(d$contours, 1L)
  z <- d$contours[[1]]; mid <- z[z[, 1] > 14 & z[, 1] < 50, 2]
  expect_true(all(abs(mid - (23.5 + width / 2)) <= 1))
 }
 img <- matrix(0.85, 48, 64); img[20:29, 8:56] <- 0.3
 v <- lx_vertices(trace_contours(extract_lines(img)))
 mid <- v[v[, 1] > 14 & v[, 1] < 50, 2]
 expect_setequal(unique(abs(mid - 24.5) > 3), TRUE)
 expect_true(all(abs(abs(mid - 24.5) - 5) <= 1))
})

# An anti-aliased stroke of the given width, 46 px long, through a sub-pixel
# centre, at an angle.
lx_stroke <- function(width, angle, fg, bg, off = c(0.3, 0.2)) {
 cx <- 36.5 + off[1]; cy <- 36.5 + off[2]
 img <- lx_render(72, 72, function(x, y) {
  a <- (x - cx) * cos(angle) + (y - cy) * sin(angle)
  b <- -(x - cx) * sin(angle) + (y - cy) * cos(angle)
  abs(a) < 23 & abs(b) < width / 2
 }, fg = fg, bg = bg, ss = 6)
 list(img = img, cx = cx, cy = cy)
}
lx_stroke_fit <- function(s, angle, ...) {
 d <- trace_contours(extract_lines(s$img, ...))
 if (!length(d$contours)) return(c(n = 0, offset = NA, span = 0))
 v <- lx_vertices(d)
 a <- (v[, 1] - s$cx) * cos(angle) + (v[, 2] - s$cy) * sin(angle)
 b <- -(v[, 1] - s$cx) * sin(angle) + (v[, 2] - s$cy) * cos(angle)
 c(n = length(d$contours), offset = max(abs(b[abs(a) < 18])), span = diff(range(a)))
}

test_that("strokes at any angle are drawn once, along their middle", {
 angles <- c(0, 0.2, 0.45, pi / 4, 1.1, 1.35, pi / 2)
 offsets <- list(c(0.3, 0.2), c(-0.45, 0.1))
 for (angle in angles) for (off in offsets) {
  # Dark strokes 2 and 3 px wide on light ground, and light strokes 1 to 3 px
  # wide on dark ground: one contour, within 1 px of the stroke's middle,
  # covering most of the 46 px stroke (measured: at least 36 px).
  for (case in list(c(2, 0.3, 0.85), c(3, 0.3, 0.85), c(1, 0.85, 0.3), c(2, 0.85, 0.3), c(3, 0.85, 0.3))) {
   fit <- lx_stroke_fit(lx_stroke(case[1], angle, case[2], case[3], off), angle)
   expect_identical(unname(fit["n"]), 1, info = paste(case, angle, collapse = " "))
   expect_lte(fit["offset"], 1)
   expect_gt(fit["span"], 34)
  }
  # A 1 px dark stroke is kept without the pre-filter at any angle.
  fit <- lx_stroke_fit(lx_stroke(1, angle, 0.3, 0.85, off), angle, smooth = 0)
  expect_identical(unname(fit["n"]), 1)
  expect_lte(fit["offset"], 1)
 }
 # Documented limit: the default pre-filter removes an oblique 1 px dark
 # stroke, but keeps one along a pixel row.
 expect_identical(unname(lx_stroke_fit(lx_stroke(1, 0.45, 0.3, 0.85), 0.45)["n"]), 0)
 expect_identical(unname(lx_stroke_fit(lx_stroke(1, 0, 0.3, 0.85), 0)["n"]), 1)
})

test_that("parallel thin strokes close together are drawn as outer outlines, as documented", {
 hatch <- function(gap) {
  img <- matrix(0.85, 48, 64); img[20:21, 8:56] <- 0.3; img[(22 + gap):(23 + gap), 8:56] <- 0.3
  d <- trace_contours(extract_lines(img))
  sort(unique(unlist(lapply(d$contours, function(z) z[z[, 1] > 14 & z[, 1] < 50, 2]))))
 }
 # 2 to 4 px apart: lines only outside the pair, none between the strokes.
 for (gap in 2:4) {
  rows <- hatch(gap)
  expect_true(length(rows) > 0 && all(rows < 20 | rows > 23 + gap))
 }
 # 6 px apart: one centre line per stroke.
 rows <- hatch(6)
 expect_length(rows, 2L)
 expect_true(all(abs(rows - c(20.5, 28.5)) <= 1))
})

test_that("quarter turns, reflections and translations transform the mask exactly", {
 set.seed(42)
 base <- lx_render(80, 72, function(x, y) (x - 30.4)^2 + (y - 33.1)^2 < 17^2 |
  (abs((x - 52) * cos(0.5) + (y - 50) * sin(0.5)) < 11 & abs(-(x - 52) * sin(0.5) + (y - 50) * cos(0.5)) < 11))
 grey <- pmin(pmax(base + matrix(runif(80 * 72, -0.08, 0.08), 80, 72), 0), 1)
 rgb <- array(c(grey, pmin(1, grey * 0.8 + 0.15), rev(grey) * 0.5 + grey * 0.5), c(80, 72, 3))
 cases <- list(list(grey, list()), list(grey, list(smooth = 0, sigma_m = 0)),
               list(rgb, list(colour = TRUE)))
 for (case in cases) {
  img <- case[[1]]; args <- case[[2]]
  run <- function(x) lx_plain(do.call(extract_lines, c(list(x), args)))
  m <- run(img)
  expect_gt(sum(m), 100)
  flip_x <- if (length(dim(img)) == 3L) img[, dim(img)[2]:1, , drop = FALSE] else img[, ncol(img):1]
  flip_y <- if (length(dim(img)) == 3L) img[dim(img)[1]:1, , , drop = FALSE] else img[nrow(img):1, ]
  tr <- if (length(dim(img)) == 3L) aperm(img, c(2, 1, 3)) else t(img)
  expect_identical(run(lx_rot90(img)), lx_rot90(m))
  expect_identical(run(flip_x), m[, ncol(m):1])
  expect_identical(run(flip_y), m[nrow(m):1, ])
  expect_identical(run(tr), t(m))
 }
 # Translation: the same pattern on a constant canvas at two offsets.
 pattern <- lx_render(40, 40, function(x, y) (x - 20.3)^2 + (y - 19.6)^2 < 12^2)
 canvas <- function(r, c) { out <- matrix(0.9, 90, 90); out[r + 1:40, c + 1:40] <- pattern; out }
 a <- lx_plain(extract_lines(canvas(10, 12)))
 b <- lx_plain(extract_lines(canvas(37, 41)))
 expect_gt(sum(a), 50)
 expect_identical(b[37 + 1:40, 41 + 1:40], a[10 + 1:40, 12 + 1:40])
 expect_identical(sum(a), sum(b))
})

test_that("the decision depends on lightness ratios, as documented", {
 core <- function(L, ...) {
  p <- mlvtoolbox:::.line_parameters(list(...))
  lx_plain(mlvtoolbox:::.line_core(list(L = L), nrow(L), ncol(L), p))
 }
 # A lightness step: drawn when the step is 40 percent of the darker side,
 # not when it is 10 percent (documented threshold about 20 percent).
 step <- function(dark, rel) { L <- matrix(dark * (1 + rel), 48, 48); L[, 1:24] <- dark; L }
 expect_gt(sum(core(step(0.4, 0.4), smooth = 0)), 0)
 expect_identical(sum(core(step(0.4, 0.1), smooth = 0)), 0L)
 # Without the pre-filter, scaling lightness leaves the mask unchanged.
 L <- lx_render(64, 64, function(x, y) (x - 30.2)^2 + (y - 33.7)^2 < 20^2, fg = 0.3, bg = 0.7)
 for (scale in c(0.5, 1.3)) expect_identical(core(L * scale, smooth = 0), core(L, smooth = 0))
})

test_that("extraction is deterministic and records its definition", {
 img <- lx_render(48, 48, function(x, y) (x - 24.5)^2 + (y - 24.5)^2 < 14^2)
 a <- extract_lines(img, sigma = 1.5, smooth = 3, min_component = 5)
 b <- extract_lines(img, sigma = 1.5, smooth = 3, min_component = 5)
 expect_identical(a, b)
 expect_type(a, "logical")
 expect_identical(dim(a), dim(img))
 e <- attr(a, "extraction")
 expect_identical(e$method, "fdog-xdog-1")
 expect_identical(e$parameters, list(sigma = 1.5, k = 1.6, tau = 0.98, eps = 0, phi = Inf, sigma_c = 2,
                                     sigma_m = 3, smooth = 3, colour = FALSE, min_component = 5L))
 # One set of defaults: the formals of extract_lines() are the defaults the
 # other entry points and feature specifications use.
 defaults <- lapply(formals(extract_lines)[-(1:2)], eval)
 expect_identical(defaults, mlvtoolbox:::.line_defaults)
 expect_identical(attr(extract_lines(img), "extraction")$parameters, mlvtoolbox:::.line_defaults)
 expect_identical(attr(xdog_lines(img), "extraction")$parameters, mlvtoolbox:::.line_defaults)
 expect_identical(line_extraction_stages(img)$extraction$parameters, mlvtoolbox:::.line_defaults)
 expect_identical(xdog_lines(img, sigma = 1.5, smooth = 3, min_component = 5), a)
})

test_that("greyscale, RGB, RGBA and file inputs are accepted", {
 grey <- lx_render(40, 40, function(x, y) abs(x - 20.2) < 9 & abs(y - 19.7) < 9)
 rgb <- array(grey, c(40, 40, 3))
 expect_identical(extract_lines(grey), extract_lines(rgb))
 expect_identical(extract_lines(array(grey, c(40, 40, 1))), extract_lines(grey))
 expect_identical(extract_lines(array(c(rgb, rep(1, 1600)), c(40, 40, 4))), extract_lines(rgb))
 transparent <- array(c(rgb, rep(0, 1600)), c(40, 40, 4))
 expect_identical(sum(extract_lines(transparent)), 0L)
 skip_if_not_installed("png")
 file <- tempfile(fileext = ".png"); png::writePNG(rgb, file)
 expect_identical(extract_lines(file), extract_lines(png::readPNG(file)))
})

test_that("constant, tiny and degenerate images give an empty mask", {
 for (v in c(0, 0.5, 1)) expect_identical(sum(extract_lines(matrix(v, 30, 30))), 0L)
 for (d in list(c(1, 1), c(2, 3), c(3, 1), c(5, 5))) {
  set.seed(prod(d))
  m <- extract_lines(matrix(runif(prod(d)), d[1], d[2]))
  expect_identical(dim(m), as.integer(d))
  expect_identical(sum(m), 0L)
 }
 expect_identical(sum(extract_lines(array(0.3, c(20, 20, 3)), colour = TRUE)), 0L)
})

test_that("an isoluminant colour boundary is not drawn, as documented", {
 lab <- function(rgb) mlvtoolbox:::.srgb_lab(rgb[1], rgb[2], rgb[3], TRUE)
 red <- c(0.8, 0.3, 0.3)
 target <- lab(red)$L
 g <- uniroot(function(g) lab(c(0.3, g, 0.3))$L - target, c(0, 1), tol = 1e-12)$root
 img <- array(0, c(40, 40, 3))
 for (ch in 1:3) { img[, 1:20, ch] <- red[ch]; img[, 21:40, ch] <- c(0.3, g, 0.3)[ch] }
 expect_lt(abs(lab(c(0.3, g, 0.3))$L - target), 1e-9)
 expect_gt(abs(lab(c(0.3, g, 0.3))$a - lab(red)$a), 0.5)
 expect_identical(sum(extract_lines(img, colour = TRUE)), 0L)
 expect_identical(sum(extract_lines(img)), 0L)
})

test_that("invalid arguments fail with clear messages", {
 img <- matrix(0.5, 10, 10)
 expect_error(extract_lines(img, sigma = 0), "sigma")
 expect_error(extract_lines(img, sigma = NA), "sigma")
 expect_error(extract_lines(img, sigma = c(1, 2)), "sigma")
 expect_error(extract_lines(img, k = 1), "k must")
 expect_error(extract_lines(img, tau = 1.5), "tau")
 expect_error(extract_lines(img, eps = Inf), "eps")
 expect_error(extract_lines(img, phi = 0), "phi")
 expect_error(extract_lines(img, sigma_c = -1), "sigma_c")
 expect_error(extract_lines(img, sigma_m = -1), "sigma_m")
 expect_error(extract_lines(img, smooth = -1), "smooth")
 expect_error(extract_lines(img, smooth = 100), "smooth")
 expect_error(extract_lines(img, colour = NA), "colour")
 expect_error(extract_lines(img, min_component = 0.5), "min_component")
 expect_error(extract_lines(img, method = "canny"), "must be \"fdog\"")
 expect_error(extract_lines(img * 3), "\\[0, 1\\]")
 expect_error(extract_lines(replace(img, 1, NA)), "\\[0, 1\\]")
 expect_error(extract_lines(array(0.5, c(4, 4, 2))), "greyscale")
 expect_error(extract_lines(list(1)), "greyscale")
 expect_error(extract_lines(matrix("a", 2, 2)), "greyscale")
 expect_error(extract_lines(c("a.png", "b.png")), "single")
 expect_error(line_extraction_stages(img, sharpen = 2), "unknown line extraction parameter")
 # The three entry points handle arguments the same way: named in full,
 # method exactly "fdog".
 expect_identical(xdog_lines(img, method = "fdog"), extract_lines(img))
 expect_identical(line_extraction_stages(img, method = "fdog")$mask, extract_lines(img))
 for (f in list(xdog_lines, line_extraction_stages)) {
  expect_error(f(img, 2), "must be named")
  expect_error(f(img, sig = 2), "unknown line extraction parameter")
  expect_error(f(img, method = "fd"), "must be \"fdog\"")
  expect_error(f(img, method = NA), "must be \"fdog\"")
 }
 expect_error(extract_lines(img, method = "fd"), "must be \"fdog\"")
 # Lower bounds keep every kernel well defined.
 expect_error(extract_lines(img, sigma = 0.05), "sigma must be")
 expect_error(extract_lines(img, sigma = 1e-200), "sigma must be")
 expect_error(extract_lines(img, sigma_c = 0.01), "sigma_c must be")
 expect_error(extract_lines(img, sigma_m = 0.05), "sigma_m must be a single number equal to 0 or")
 expect_error(extract_lines(img, smooth = 0.2), "smooth must be a single number equal to 0 or")
 expect_silent(extract_lines(img, sigma = 0.1, sigma_c = 0.1, sigma_m = 0.1, smooth = 0.5))
 expect_silent(extract_lines(img, sigma_m = 0, smooth = 0))
})

test_that("images with a single row or column are accepted", {
 for (d in list(c(1, 5, 3), c(5, 1, 3), c(1, 1, 3), c(1, 5, 4), c(1, 5, 1))) {
  m <- extract_lines(array(0.5, d))
  expect_identical(dim(m), as.integer(d[1:2]))
  expect_identical(sum(m), 0L)
 }
 expect_identical(dim(extract_lines(array(0.5, c(1, 5, 3)), colour = TRUE)), c(1L, 5L))
 expect_identical(dim(line_extraction_stages(array(0.5, c(5, 1, 3)))$lightness), c(5L, 1L))
 d <- trace_rgb(array(0.5, c(1, 6, 3)), xdog_lines)
 expect_length(d$contours, 0L)
})

test_that("a threshold above the whole response gives an empty mask, as documented", {
 img <- lx_render(48, 48, function(x, y) (x - 24.3)^2 + (y - 23.8)^2 < 14^2)
 s <- line_extraction_stages(img, eps = 1)
 expect_true(all(s$band))
 expect_identical(sum(s$mask), 0L)
})

test_that("line_extraction_stages returns every stage and the same mask", {
 img <- lx_render(40, 48, function(x, y) (x - 24.4)^2 + (y - 20.3)^2 < 12^2)
 s <- line_extraction_stages(img, min_component = 5)
 expect_named(s, c("lightness", "smoothed", "flow", "dog", "response", "xdog", "band", "mask", "extraction"))
 for (k in c("lightness", "smoothed", "dog", "response", "xdog", "band")) expect_identical(dim(s[[k]]), c(40L, 48L))
 expect_named(s$flow, c("x", "y", "anisotropy"))
 expect_true(all(abs(s$flow$x^2 + s$flow$y^2 - 1) < 1e-12 | (s$flow$x == 0 & s$flow$y == 0)))
 expect_true(all(s$flow$anisotropy >= 0 & s$flow$anisotropy <= 1))
 expect_identical(s$mask, extract_lines(img, min_component = 5))
 expect_identical(s$extraction, attr(s$mask, "extraction"))
 expect_true(all(s$xdog %in% c(0, 1)))
 expect_identical(s$band, s$response < 0)
 soft <- line_extraction_stages(img, phi = 50)
 expect_true(all(soft$xdog >= 0 & soft$xdog <= 1) && any(soft$xdog > 0 & soft$xdog < 1))
 # Lightness is CIE L*/100: white is 1, black is 0, sRGB mid-grey is 0.5338.
 expect_equal(line_extraction_stages(matrix(c(0, 1, 0.5, 0.5), 2))$lightness[, 1], c(0, 1))
 expect_equal(line_extraction_stages(matrix(0.5, 2, 2))$lightness[1, 1], 0.53389, tolerance = 1e-4)
})

test_that("trace_rgb records the line extractor's provider and parameters", {
 rgb <- array(0.9, c(40, 40, 3)); rgb[10:30, 12:28, ] <- 0.25
 d <- trace_rgb(rgb, backend = xdog_lines, min_component = 10)
 expect_identical(d$source$extraction_provider, "xdog_lines")
 expect_identical(d$source$provider_version, "fdog-xdog-1")
 expect_identical(d$source$extraction_parameters$min_component, 10L)
 expect_length(d$contours, 1L)
 # A backend without the attribute keeps the previous provenance.
 plain <- trace_rgb(rgb, function(image) image[, , 1] < 0.5)
 expect_null(plain$source$provider_version)
 expect_null(plain$source$extraction_parameters)
 bad <- function(image) structure(image[, , 1] < 0.5, extraction = list(method = 3))
 expect_error(trace_rgb(rgb, bad), "single method string")
})

test_that("feature specifications carry and validate the line extraction definition", {
 s <- feature_spec(extraction = list(lines = list(sigma = 1.5, min_component = 10)))
 expect_identical(unname(s$definitions["lines"]), "fdog-xdog-1")
 expect_identical(s$settings$extraction$provider, "xdog_lines")
 expect_identical(s$settings$extraction$provider_version, "fdog-xdog-1")
 expect_identical(s$settings$extraction$lines$sigma, 1.5)
 expect_identical(s$settings$extraction$lines$smooth, 5)
 expect_invisible(validate_feature_spec(s))
 expect_identical(feature_spec(extraction = list(provider = "xdog_lines", provider_version = "fdog-xdog-1"))$settings$extraction$lines,
                  mlvtoolbox:::.line_parameters(list()))
 # Specifications without line settings are unchanged.
 expect_false("lines" %in% names(feature_spec()$settings$extraction))
 expect_false("lines" %in% names(feature_spec()$definitions))
 # Stale or unknown definitions are rejected.
 expect_error(feature_spec(extraction = list(provider = "xdog_lines", provider_version = "fdog-xdog-0")), "unsupported line extraction definition")
 expect_error(feature_spec(extraction = list(provider = "xdog_lines", provider_version = "2")), "requires provider_version")
 expect_error(feature_spec(extraction = list(lines = list(sigma = 1.5), provider_version = "fdog-xdog-0", provider = "x")), "unsupported")
 expect_error(feature_spec(extraction = list(lines = list(blur = 2))), "unknown line extraction parameter")
 expect_error(feature_spec(extraction = list(lines = list(sigma = -1))), "sigma")
 expect_error(feature_spec(extraction = list(lines = 2)), "named list")
 expect_error(feature_spec(extraction = list(lines = list(), provider = "mine", provider_version = "v1")),
              "requires provider_version 'fdog-xdog-1'")
 stale <- s; stale$definitions["lines"] <- "fdog-xdog-0"
 expect_error(validate_feature_spec(stale), "noncanonical")
 edited <- s; edited$settings$extraction$lines$sigma <- 3
 expect_invisible(validate_feature_spec(edited))
 partial <- s; partial$settings$extraction$lines$sigma <- NULL
 expect_error(validate_feature_spec(partial), "noncanonical")
 # trace_rgb uses the specification's parameters and checks the backend's output.
 rgb <- array(0.9, c(40, 40, 3)); rgb[10:30, 12:28, ] <- 0.25
 d <- trace_rgb(rgb, xdog_lines, spec = s)
 expect_identical(d$source$provider_version, "fdog-xdog-1")
 expect_identical(d$source$extraction_parameters, s$settings$extraction$lines)
 expect_identical(d$feature_spec, s)
 expect_error(trace_rgb(rgb, xdog_lines, spec = s, sigma = 3), "spec or line extraction arguments")
 expect_error(trace_rgb(rgb, function(image, ...) image[, , 1] < 0.5, spec = s), "does not match")
 # A specification naming another provider version cannot record this
 # backend's output as its own.
 other <- feature_spec(extraction = list(provider = "mine", provider_version = "v1"))
 expect_error(trace_rgb(rgb, xdog_lines, spec = other), "reports extraction 'fdog-xdog-1'")
})
