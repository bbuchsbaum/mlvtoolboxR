# The `arms` column of detect_junctions(). Expected arms are derived by hand
# from the construction. Drawing convention: counterclockwise degrees, image y
# down, so a segment drawn up the screen has orientation 90.
jd <- function(...) mlv_drawing(list(...), 100, 100)
sg <- function(...) matrix(c(...), ncol = 4, byrow = TRUE)
# Gaps from arms as documented: angles[k] = arms[k+1] - arms[k], and the last
# gap wraps: arms[1] + 360 - arms[n].
gaps_of <- function(a) diff(c(a, a[1] + 360))

test_that("default-rule arms are the directions in which contours leave the junction", {
  # T: bar leaves right (0) and left (180); the stem leaves down (270).
  j <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 10, 10, 20)))
  expect_equal(j$arms[[1]], c(0, 180, 270))
  expect_equal(j$angles[[1]], c(180, 90, 90))
  # X: two crossing lines leave in all four axis directions.
  j <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 0, 10, 20)))
  expect_equal(j$arms[[1]], c(0, 90, 180, 270))
  expect_equal(j$angles[[1]], rep(90, 4))
  # Diagonal T: the bar (10,10)-(30,30) leaves up-left (135) and down-right
  # (315); the stem (20,20)-(30,10) leaves up-right (45).
  j <- detect_junctions(jd(sg(10, 10, 30, 30), sg(20, 20, 30, 10)))
  expect_equal(j$arms[[1]], c(45, 135, 315), tolerance = 1e-12)
  expect_equal(j$angles[[1]], c(90, 180, 90), tolerance = 1e-12)
})

test_that("mlv2023 arms are MLV's own directions, including the M007 convention", {
  # Horizontal bar through the junction: MLV takes the direction from each
  # bar end towards the junction with y up: from (0,10) atan2(0, 10) = 0,
  # from (20,10) atan2(0, -10) = 180. The stem starts at the junction and
  # contributes its orientation, 270. Same as the default rule here.
  j <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 10, 10, 20)), rule = "mlv2023")
  expect_equal(j$arms[[1]], c(0, 180, 270))
  expect_equal(j$angles[[1]], c(180, 90, 90))
  # X: the vertical line gives atan2(10, 0) = 90 from (10,0) and
  # atan2(-10, 0) = 270 from (10,20).
  j <- detect_junctions(jd(sg(0, 10, 20, 10), sg(10, 0, 10, 20)), rule = "mlv2023")
  expect_equal(j$arms[[1]], c(0, 90, 180, 270))
  # Diagonal T (MLV-M007): from (10,10), atan2(20 - 10, 20 - 10) = 45; from
  # (30,30), atan2(-10, -10) = -135, i.e. 225. The stem contributes its
  # orientation 45. Arms 45, 45, 225; gaps 0, 180, 180.
  j <- detect_junctions(jd(sg(10, 10, 30, 30), sg(20, 20, 30, 10)), rule = "mlv2023")
  expect_equal(j$arms[[1]], c(45, 45, 225), tolerance = 1e-12)
  expect_equal(j$angles[[1]], c(0, 180, 180), tolerance = 1e-12)
  expect_equal(j$min_angle, 0, tolerance = 1e-12)
})

random_drawing <- function(seed) {
  set.seed(seed)
  mlv_drawing(lapply(seq_len(sample(5:30, 1)), function(k) {
    p <- runif(2, 0, 30); m <- sample(1:3, 1)
    xy <- rbind(p, p + apply(matrix(rnorm(2 * m, 0, 4), m, 2), 2, cumsum))
    cbind(xy[-nrow(xy), , drop = FALSE], xy[-1, , drop = FALSE])
  }), 40, 40)
}

test_that("arms reproduce angles on random drawings under both rules", {
  seen <- 0L
  for (seed in 4100 + seq_len(40)) {
    d <- random_drawing(seed)
    for (rule in c("exit", "mlv2023")) {
      j <- detect_junctions(d, rule = rule)
      expect_identical(lengths(j$arms), lengths(j$angles))
      for (k in seq_len(nrow(j))) {
        a <- j$arms[[k]]
        seen <- seen + 1L
        expect_false(is.unsorted(a))
        expect_true(all(a >= 0 & a <= 360))
        if (rule == "exit") {
          expect_equal(j$angles[[k]], gaps_of(a), tolerance = 1e-9)
          expect_equal(sum(j$angles[[k]]), 360, tolerance = 1e-9)
        } else {
          expect_equal(j$angles[[k]], gaps_of(a) %% 360, tolerance = 1e-9)
        }
        expect_equal(j$min_angle[k], min(j$angles[[k]]))
        expect_equal(j$max_angle[k], max(j$angles[[k]]))
      }
    }
  }
  expect_gt(seen, 50L)
})

test_that("default-rule arms do not depend on contour direction or numbering", {
  key <- function(j) order(round(j$x, 6), round(j$y, 6))
  compared <- 0L
  for (seed in 4200 + seq_len(30)) {
    d <- random_drawing(seed)
    set.seed(seed + 1L)
    perm <- sample(length(d$contours))
    flipped <- lapply(d$contours[perm], function(z) z[nrow(z):1, c(3, 4, 1, 2), drop = FALSE])
    e <- mlv_drawing(flipped, d$width, d$height)
    a <- detect_junctions(d); b <- detect_junctions(e)
    expect_equal(nrow(a), nrow(b))
    if (!nrow(a)) next
    a <- a[key(a), ]; b <- b[key(b), ]
    expect_equal(a$x, b$x, tolerance = 1e-9)
    expect_equal(a$y, b$y, tolerance = 1e-9)
    for (k in seq_len(nrow(a))) {
      expect_equal(a$arms[[k]], b$arms[[k]], tolerance = 1e-9)
      compared <- compared + 1L
    }
  }
  expect_gt(compared, 20L)
})

test_that("empty junction tables carry the arms column under both rules", {
  for (rule in c("exit", "mlv2023")) {
    empty <- detect_junctions(jd(sg(0, 0, 10, 10)), rule = rule)
    expect_equal(nrow(empty), 0L)
    expect_true(all(c("arms", "angles", "members") %in% names(empty)))
  }
  none <- detect_junctions(jd(sg(0, 0, 10, 0), sg(0, 20, 10, 20)))
  expect_equal(nrow(none), 0L)
  expect_true("arms" %in% names(none))
})
