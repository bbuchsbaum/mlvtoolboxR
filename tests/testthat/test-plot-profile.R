# plot_branch_profile(): branch quantities against image x, under the drawing.

profile_fixture <- function() {
  d <- mlv_drawing(list(matrix(c(5, 8, 75, 8), 1), matrix(c(5, 30, 75, 22), 1)), 80, 36)
  compute_mat_properties(compute_mat(d), c("parallelism", "separation"))
}

with_png <- function(code, width = 600, height = 500) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = width, height = height)
  on.exit({ grDevices::dev.off(); unlink(file) })
  force(code)
}

test_that("the returned data are the branch values that were plotted", {
  m <- profile_fixture()
  with_png(out <- plot_branch_profile(m, c("radius", "parallelism")))
  expect_s3_class(out, "data.frame")
  expect_named(out, c("column", "branch", "what", "x", "y", "value"))
  expect_setequal(unique(out$what), c("radius", "parallelism"))
  for (i in seq_along(m$branches)) {
    b <- m$branches[[i]]
    got <- out[out$branch == i & out$what == "parallelism", ]
    expect_equal(got$x, b$x)
    expect_equal(got$value, b$parallelism)
  }
})

test_that("branch selection, multiple columns and the panel hook work", {
  m <- profile_fixture()
  longest <- which.max(vapply(m$branches, nrow, integer(1)))
  calls <- character()
  hook <- function(column, row) calls <<- c(calls, paste(column, row))
  with_png(out <- plot_branch_profile(list(a = m, b = m), "separation",
    branches = "longest", discs = 5, col = c("red", "blue"), band = c(20, 60),
    ylim = list(separation = c(0, 0.5)), panel = hook))
  expect_identical(unique(out$branch), longest)
  expect_identical(sort(unique(out$column)), 1:2)
  expect_identical(calls, c("1 drawing", "1 separation", "2 drawing", "2 separation"))
  with_png(out <- plot_branch_profile(m, "parallelism", branches = list(1L)))
  expect_identical(unique(out$branch), 1L)
})

test_that("profiles share the image's x coordinates", {
  m <- profile_fixture()
  xr <- list()
  with_png(plot_branch_profile(m, "parallelism", panel = function(column, row)
    xr[[row]] <<- par("usr")[1:2]))
  expect_equal(xr$parallelism, xr$drawing)
  expect_lte(xr$drawing[1], 0.5); expect_gte(xr$drawing[2], ncol(m$mask) + 0.5)
})

test_that("graphics state is restored and invalid input fails clearly", {
  m <- profile_fixture()
  with_png({
    par(mfrow = c(1, 2), mar = c(2, 3, 4, 1))
    before <- par(c("mfrow", "mar", "cex", "las"))
    plot_branch_profile(m, "parallelism")
    expect_identical(par(c("mfrow", "mar", "cex", "las")), before)
    expect_error(plot_branch_profile(list(), "parallelism"), "mlv_mat")
    expect_error(plot_branch_profile(m, "taper"), "no column taper")
    expect_error(plot_branch_profile(m, NA_character_), "what")
    expect_error(plot_branch_profile(m, branches = 99), "branches")
    expect_error(plot_branch_profile(list(m, m), branches = list(1L)), "one entry per column")
    expect_error(plot_branch_profile(m, discs = -1), "discs")
    expect_error(plot_branch_profile(m, band = 1), "band")
    expect_error(plot_branch_profile(m, panel = 1), "panel")
  })
})

test_that("an unscored MAT and a MAT without branches are handled", {
  single <- compute_mat(mlv_drawing(list(matrix(c(10, 20, 50, 20), 1)), 60, 40))
  expect_length(single$branches, 0L)
  with_png(out <- plot_branch_profile(single, "radius"))
  expect_identical(nrow(out), 0L)
  unscored <- compute_mat(mlv_drawing(list(matrix(c(5, 8, 75, 8), 1),
                                           matrix(c(5, 30, 75, 22), 1)), 80, 36))
  expect_error(with_png(plot_branch_profile(unscored, "parallelism")), "compute_mat_properties")
})

test_that("blank rows, context curves and equal drawing heights work", {
  m <- profile_fixture()
  calls <- character()
  with_png(out <- plot_branch_profile(list(a = m, b = m), c("notes", "parallelism"),
    blank = "notes", context = TRUE,
    panel = function(column, row) calls <<- c(calls, paste(column, row))))
  expect_identical(calls, c("1 drawing", "1 notes", "1 parallelism",
                            "2 drawing", "2 notes", "2 parallelism"))
  expect_false("notes" %in% out$what)
  with_png({
    expect_error(plot_branch_profile(m, "notes"), "no column notes")
    expect_error(plot_branch_profile(m, context = NA), "context")
    expect_error(plot_branch_profile(m, blank = 1), "blank")
  })
  # Drawings of different canvas heights get equal panel heights.
  tall <- compute_mat_properties(compute_mat(mlv_drawing(list(
    matrix(c(5, 4, 75, 4), 1), matrix(c(5, 40, 75, 30), 1)), 80, 46)), "parallelism")
  heights <- numeric()
  with_png(plot_branch_profile(list(m, tall), "parallelism", panel = function(column, row)
    if (row == "drawing") heights[column] <<- abs(diff(par("usr")[3:4]))), width = 700)
  expect_equal(heights[1], heights[2], tolerance = 1e-6)
})

test_that("the medial axis is drawn in its own colour, never as ink", {
  skip_if_not_installed("png")
  m <- profile_fixture()
  b <- m$branches[[which.max(vapply(m$branches, nrow, integer(1)))]]
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 800, height = 700)
  where <- list()
  plot_branch_profile(m, "parallelism", branches = "longest", col = "black",
    panel = function(column, row) if (row == "drawing") {
      k <- which.min(abs(b$x - 40))
      ink <- which(m$mask, arr.ind = TRUE)[1, ]
      where$axis <<- c(grconvertX(b$x[k], "user", "device"), grconvertY(b$y[k], "user", "device"))
      where$ink <<- c(grconvertX(ink[2], "user", "device"), grconvertY(ink[1], "user", "device"))
    })
  grDevices::dev.off(); on.exit(unlink(file))
  img <- png::readPNG(file)
  px <- function(p) img[round(p[2]), round(p[1]), 1:3]
  axis <- px(where$axis)
  expect_false(all(axis < 0.2))                       # not black
  expect_gt(axis[1], axis[2] + 0.3)                   # red
  expect_true(all(px(where$ink) < 0.1))               # ink is black
})

test_that("stacked columns, row-specific context and the axis helper work", {
  m <- profile_fixture()
  calls <- character()
  with_png(plot_branch_profile(list(a = m, b = m), "parallelism", stack = TRUE,
    context = "parallelism", panel = function(column, row) calls <<- c(calls, row)),
    height = 900)
  expect_identical(calls, c("drawing", "parallelism", "drawing", "parallelism"))
  with_png({
    expect_error(plot_branch_profile(m, stack = NA), "stack")
    expect_error(plot_branch_profile(m, context = 1), "context")
    plot(m, "ink", legend = FALSE)
    expect_invisible(draw_mat_axis(m))
    expect_no_error(draw_mat_axis(m, branches = 1, col = "grey40", width = 1))
    expect_error(draw_mat_axis(list()), "mlv_mat")
    expect_error(draw_mat_axis(m, width = 0), "width")
  })
})
