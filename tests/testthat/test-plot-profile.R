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
