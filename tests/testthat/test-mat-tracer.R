# Default branch tracer (branch_rule = "invariant", definition
# "direction-symmetric-tiemax-graph-2" onwards): order-free spur pruning that never
# deletes a skeleton component, merging through degree-two pixels.

poly <- function(xy) { xy <- rbind(xy, xy[1, ]); cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE]) }
fake_mat <- function(sk) structure(list(skeleton = sk, distance_map = matrix(2, nrow(sk), ncol(sk)),
  AOF = matrix(0, nrow(sk), ncol(sk)), parameters = list(branch_rule = "invariant")), class = "mlv_mat")
paths_of <- function(b) lapply(b, function(z) paste(z$row, z$column))
# Branch set up to direction (and rotation of closed paths), as sorted strings.
canon <- function(keys) sort(vapply(keys, function(key) {
  if (length(key) > 3 && key[1] == key[length(key)]) { u <- key[-length(key)]
    min(unlist(lapply(seq_along(u), function(s) { v <- u[c(s:length(u), seq_len(s - 1L))]; v <- c(v, v[1])
      c(paste(v, collapse = ";"), paste(rev(v), collapse = ";")) })))
  } else min(paste(key, collapse = ";"), paste(rev(key), collapse = ";")) }, ""))
components <- function(s) { lab <- matrix(0L, nrow(s), ncol(s)); k <- 0L
  for (i in which(s)) if (!lab[i]) { k <- k + 1L; q <- i; lab[i] <- k
   while (length(q)) { p <- q[1]; q <- q[-1]; r <- (p - 1) %% nrow(s) + 1; c <- (p - 1) %/% nrow(s) + 1
    for (dr in -1:1) for (dc in -1:1) { rr <- r + dr; cc <- c + dc
     if (rr < 1 || cc < 1 || rr > nrow(s) || cc > ncol(s)) next; j <- rr + nrow(s) * (cc - 1)
     if (s[j] && !lab[j]) { lab[j] <- k; q <- c(q, j) } } } }; lab }
small_shapes <- function() {
  o <- 8; sh <- function(xy) mlv_drawing(list(poly(xy)), ceiling(max(xy[, 1])) + o, ceiling(max(xy[, 2])) + o)
  sq <- function(s) sh(cbind(o + c(0, s, s, 0), o + c(0, 0, s, s)))
  rect <- function(a, b) sh(cbind(o + c(0, a, a, 0), o + c(0, 0, b, b)))
  circ <- function(d) { t <- seq(0, 2 * pi, length = 41)[-1]; sh(cbind(o + d / 2 * (1 + cos(t)), o + d / 2 * (1 + sin(t)))) }
  ring <- function(d, w) { t <- seq(0, 2 * pi, length = 41)[-1]
    mlv_drawing(list(poly(cbind(o + d / 2 * (1 + cos(t)), o + d / 2 * (1 + sin(t)))),
                     poly(cbind(o + d / 2 + (d / 2 - w) * cos(t), o + d / 2 + (d / 2 - w) * sin(t)))), d + 2 * o, d + 2 * o) }
  list(square_10 = sq(10), square_14 = sq(14), square_20 = sq(20), rect_20x10 = rect(20, 10), rect_20x14 = rect(20, 14),
       rect_20x18 = rect(20, 18), tri_14 = sh(cbind(o + c(0, 14, 7), o + c(14 * sqrt(3) / 2, 14 * sqrt(3) / 2, 0))),
       tri_20 = sh(cbind(o + c(0, 20, 10), o + c(10 * sqrt(3), 10 * sqrt(3), 0))), rtri_14 = sh(cbind(o + c(0, 14, 0), o + c(0, 14, 14))),
       circle_10 = circ(10), circle_16 = circ(16), circle_30 = circ(30), ring_20_4 = ring(20, 4), ring_40_10 = ring(40, 10),
       L_20_8 = sh(cbind(o + c(0, 8, 8, 20, 20, 0), o + c(0, 0, 12, 12, 20, 20))),
       T_20_8 = sh(cbind(o + c(0, 20, 20, 14, 14, 6, 6, 0), o + c(0, 0, 8, 8, 20, 20, 8, 8))),
       H_20_8 = sh(cbind(o + c(0, 8, 8, 12, 12, 20, 20, 12, 12, 8, 8, 0), o + c(0, 0, 6, 6, 0, 0, 20, 20, 14, 14, 20, 20))))
}

test_that("a 20 px square keeps its four tied diagonal arms and maps every contour pixel", {
  square <- mlv_drawing(list(rbind(c(10, 10, 30, 10), c(30, 10, 30, 30), c(30, 30, 10, 30), c(10, 30, 10, 10))), 40, 40)
  a <- analyse_mat(square, properties = c("parallelism", "separation", "mirror", "taper"))
  medial <- a$mat$skeleton & a$mat$distance_map > 0
  # Oracle: the skeleton is an X of four arms meeting at one junction; every arm
  # is shorter than min_branch_points = 10, and removing all four tied arms would
  # delete the component, so the junction is protected and all four are kept.
  expect_length(a$mat$branches, 4L)
  expect_true(all(vapply(a$mat$branches, nrow, 1L) < 10L))
  shared <- Reduce(intersect, paths_of(a$mat$branches))
  expect_length(shared, 1L)                                       # the four arms share only the junction
  kept <- matrix(FALSE, 40, 40); for (b in a$mat$branches) kept[cbind(b$row, b$column)] <- TRUE
  expect_identical(kept, medial)                                  # the whole skeleton is traced
  ink <- rasterize(square); expect_identical(sum(ink), 80L)
  for (p in names(a$drawing$mat_images)) expect_true(all(is.finite(a$drawing$mat_images[[p]][ink])))
  # the published chain is untouched by the default rule's change
  expect_identical(compute_mat(square, branch_rule = "mlv2023")$parameters$branch_rule_definition, "mlv2023-edgelink-lastwriter-1")
})

test_that("small closed shapes keep a nonempty skeleton in every component", {
  for (nm in names(sh <- small_shapes())) {
    m <- compute_mat(sh[[nm]]); medial <- m$skeleton & m$distance_map > 0
    lab <- components(medial); size <- tabulate(lab[lab > 0])
    expect_gt(length(m$branches), 0L, label = nm)
    covered <- unique(unlist(lapply(m$branches, function(b) lab[cbind(b$row, b$column)])))
    expect_setequal(covered, which(size >= 2L))                    # no component with a path is lost
    y <- map_mat_properties(sh[[nm]], compute_mat_properties(m))
    expect_gt(mean(is.finite(y$mat_images$parallelism[rasterize(sh[[nm]])])), .5, label = nm)
  }
})

test_that("short spurs are pruned shortest first, ties together, and long structure is kept", {
  tr <- function(sk, k = 10L) paths_of(trace_skeleton(fake_mat(sk), k))
  # a 41 px line with a 4 px perpendicular spur at its middle: the spur goes and
  # the line is merged through the former junction into one branch
  sk <- matrix(FALSE, 20, 50); sk[10, 5:45] <- TRUE; sk[11:13, 25] <- TRUE
  p <- tr(sk); expect_length(p, 1L); expect_identical(length(p[[1]]), 41L)
  expect_false(any(c("11 25", "12 25", "13 25") %in% unlist(p)))
  # a spur of at least min_branch_points pixels is not pruned
  sk2 <- sk; sk2[11:19, 25] <- TRUE; expect_length(tr(sk2), 3L)
  # the same spur is removed when it is shorter than a larger minimum
  expect_length(tr(sk2, 11L), 1L)
  # an untied fork at a line end: the shorter arm is pruned, the longer one is
  # merged into the line (sequential shortest-first pruning)
  f <- matrix(FALSE, 20, 50); f[10, 5:30] <- TRUE
  f[cbind(9:8, 31:32)] <- TRUE                                   # arm of 3 pixels with the junction
  f[cbind(11:14, 31:34)] <- TRUE                                 # arm of 5 pixels with the junction
  p <- tr(f); expect_length(p, 1L); expect_identical(length(p[[1]]), 30L)
  expect_true("14 34" %in% p[[1]]); expect_false("8 32" %in% p[[1]])
  # a tied fork (two arms of 4) would leave its junction an end: both arms kept
  g <- matrix(FALSE, 20, 50); g[10, 5:30] <- TRUE; g[cbind(9:7, 31:33)] <- TRUE; g[cbind(11:13, 31:33)] <- TRUE
  p <- tr(g); expect_length(p, 3L); expect_identical(sort(lengths(p)), c(4L, 4L, 26L))
  # a star of unequal short arms keeps the path through its two longest arms
  s <- matrix(FALSE, 30, 30); s[15, 10:19] <- TRUE; s[16:17, 15] <- TRUE   # arms of 6, 5 and 3 pixels at (15, 15)
  p <- tr(s); expect_length(p, 1L); expect_identical(length(p[[1]]), 10L)
  # a star of tied short arms is kept whole; a single short isolated path too
  x <- matrix(FALSE, 30, 30); x[15, 11:19] <- TRUE; x[11:19, 15] <- TRUE
  expect_identical(sort(lengths(tr(x))), rep(5L, 4))
  i <- matrix(FALSE, 10, 10); i[5, 3:6] <- TRUE; expect_identical(lengths(tr(i)), 4L)
  # interior junction-to-junction paths are not spurs and are never pruned
  h <- matrix(FALSE, 30, 40); h[15, 5:35] <- TRUE; h[3:27, 18] <- TRUE; h[3:27, 20] <- TRUE
  p <- tr(h); expect_length(p, 7L)
  expect_true(any(vapply(p, function(z) identical(z, c("15 18", "15 19", "15 20")) || identical(z, c("15 20", "15 19", "15 18")), TRUE)))
})

test_that("default tracing is idempotent, keeps every component and is lattice invariant on random skeletons", {
  ns <- asNamespace("mlvtoolbox")
  tf <- list(function(m) t(m)[, nrow(m):1], function(m) m[nrow(m):1, ncol(m):1], function(m) t(m)[ncol(m):1, ],
             function(m) m[nrow(m):1, ], function(m) m[, ncol(m):1], function(m) t(m), function(m) t(m[nrow(m):1, ncol(m):1]))
  keys <- function(paths, s, g = NULL) { id <- matrix(seq_along(s), nrow(s)); gid <- if (is.null(g)) id else g(id)
    lapply(paths, function(p) { rc <- arrayInd(gid[p], dim(s)); paste(rc[, 1], rc[, 2]) }) }
  set.seed(23)
  for (it in 1:120) {
    h <- sample(10:36, 1); w <- sample(10:36, 1)
    s <- skeletonize(matrix(runif(h * w) < runif(1, .2, .7), h, w)); k <- sample(3:12, 1)
    P <- ns$cpp_trace_graph_pruned(s, k); base <- canon(keys(P, s))
    for (g in tf) expect_identical(canon(keys(ns$cpp_trace_graph_pruned(g(s), k), s, g)), base)
    s2 <- s; s2[] <- FALSE; s2[unlist(P)] <- TRUE
    expect_identical(canon(keys(ns$cpp_trace_graph_pruned(s2, k), s)), base)          # pruning is at a fixed point
    lab <- components(s); size <- tabulate(lab[lab > 0])
    expect_setequal(unique(lab[unlist(P)]), which(size >= 2L))
    # junction pixels are shared, never repeated inside an open path
    for (p in P) if (p[1] != p[length(p)]) expect_false(anyDuplicated(p) > 0)
  }
})

test_that("default scores on small shapes are invariant to reversal, renumbering and lattice transforms", {
  n <- 41
  tf <- list(rot90 = function(x, y) cbind(n + 1 - y, x), rot180 = function(x, y) cbind(n + 1 - x, n + 1 - y),
             rot270 = function(x, y) cbind(y, n + 1 - x), flipx = function(x, y) cbind(n + 1 - x, y),
             flipy = function(x, y) cbind(x, n + 1 - y), transpose = function(x, y) cbind(y, x),
             anti = function(x, y) cbind(n + 1 - y, n + 1 - x))
  tfd <- function(d, g) mlv_drawing(lapply(d$contours, function(z) cbind(g(z[, 1], z[, 2]), g(z[, 3], z[, 4]))), n, n)
  tfi <- function(M, g) { out <- M; rc <- arrayInd(seq_along(M), dim(M)); p <- g(rc[, 2], rc[, 1]); out[cbind(p[, 2], p[, 1])] <- M; out }
  sq <- function(a, b, o) poly(cbind(o + c(0, a, a, 0), o + c(0, 0, b, b)))
  drawings <- list(square20 = mlv_drawing(list(sq(20, 20, 10)), n, n), rect = mlv_drawing(list(sq(20, 14, 10)), n, n),
    tri = mlv_drawing(list(poly(cbind(c(10, 30, 20), c(28, 28, 10)))), n, n),
    tee = mlv_drawing(list(poly(cbind(10 + c(0, 20, 20, 14, 14, 6, 6, 0), 10 + c(0, 0, 8, 8, 20, 20, 8, 8)))), n, n),
    ring = mlv_drawing(list(sq(24, 24, 8), sq(14, 14, 13)), n, n))
  P <- c("parallelism", "separation", "mirror", "taper")
  close_to <- function(a, b) for (p in P) { expect_identical(is.na(a[[p]]), is.na(b[[p]]))
    if (any(is.finite(a[[p]]))) expect_lt(max(abs(a[[p]] - b[[p]]), na.rm = TRUE), 1e-12) }
  for (d in drawings) {
    m <- compute_mat(d); expect_gt(length(m$branches), 0L)
    A <- map_mat_properties(d, compute_mat_properties(m, P))$mat_images
    rv <- m; rv$branches <- lapply(m$branches, function(b) b[nrow(b):1, ])
    close_to(map_mat_properties(d, compute_mat_properties(rv, P))$mat_images, A)
    set.seed(5); pm <- m; pm$branches <- m$branches[sample(length(m$branches))]
    close_to(map_mat_properties(d, compute_mat_properties(pm, P))$mat_images, A)
    m0 <- compute_mat(d, trace_branches = FALSE)
    A0 <- map_mat_properties(d, compute_mat_properties(m0, P))$mat_images
    close_to(A0, A)
    for (g in tf) {
      mt <- m0; for (f in c("distance_map", "AOF", "skeleton", "mask")) mt[[f]] <- tfi(m0[[f]], g)
      B <- map_mat_properties(tfd(d, g), compute_mat_properties(mt, P))$mat_images
      close_to(lapply(A0, tfi, g), B)
    }
  }
})

test_that("a loop left by pruning its spurs is an isolated loop with no start", {
  # A square ring of skeleton with four tied 4 px spurs pointing outwards. Each
  # junction keeps its two ring paths when its single spur goes, so all spurs are
  # pruned and the ring merges into one closed path. Its pixels include the
  # former junctions, which still touch pruned skeleton pixels; the loop is
  # recognised by its start pixel joining no other branch, not by skeleton degree.
  sk <- matrix(FALSE, 40, 40); sk[10, 10:30] <- sk[30, 10:30] <- TRUE; sk[10:30, 10] <- sk[10:30, 30] <- TRUE
  sk[15, 31:33] <- TRUE; sk[25, 7:9] <- TRUE; sk[31:33, 20] <- TRUE; sk[7:9, 22] <- TRUE
  m <- fake_mat(sk); m$distance_map[] <- 3 + (row(sk) %% 3) / 2 + (col(sk) %% 5) / 7
  d <- mlv_drawing(list(matrix(c(3, 3, 37, 3), 1)), 40, 40); m$mask <- rasterize(d)
  m$branches <- trace_skeleton(m, 10L)
  expect_length(m$branches, 1L)
  b <- m$branches[[1]]; expect_identical(nrow(b), 81L)                 # 80 ring pixels, closed
  expect_true(mlvtoolbox:::.is_isolated_loop(b, m$branches))
  P <- c("parallelism", "separation", "mirror", "taper")
  A <- compute_mat_properties(m, P)$score_images
  u <- b[-nrow(b), ]; k <- nrow(u)
  # start at a former junction pixel and elsewhere, both directions
  j <- which(paste(u$row, u$column) %in% c("15 30", "25 10", "30 20", "10 22"))
  for (s in c(j, 17L)) for (rv in c(FALSE, TRUE)) {
    ord <- c(s:k, seq_len(s - 1L)); if (rv) ord <- rev(ord)
    m2 <- m; m2$branches <- list(u[c(ord, ord[1]), ])
    B <- compute_mat_properties(m2, P)$score_images
    for (p in P) expect_lt(max(abs(A[[p]] - B[[p]]), na.rm = TRUE), 1e-12)
  }
})
