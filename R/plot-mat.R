#' Plot a stage of a medial-axis representation
#'
#' Draws one stage of the pipeline stored in an `mlv_mat`, in image
#' coordinates: x increases to the right, y increases downward, pixel centres
#' lie at integers and the aspect ratio is one. Ink is black in every stage.
#' Quantities are colour-coded with a legend that names the quantity and its
#' unit:
#'
#' * `"ink"`: the rasterized ink mask (black on white).
#' * `"distance"`: the signed distance map in pixels (sequential palette), with the ink in black.
#'   Ink pixels have nonpositive values; only an isolated single ink pixel is
#'   negative.
#' * `"flux"`: the average outward flux, stored as the sum of `samples`
#'   dot products (diverging palette with limits symmetric about zero, zero
#'   drawn as a neutral near-white; ink in black). The skeleton threshold,
#'   `threshold * samples`, is marked on the colour bar and drawn as a dashed
#'   isoline in the image.
#' * `"skeleton"`: the pixels whose flux exceeds the threshold (light red) and
#'   the thinned skeleton (dark red), with the ink in black.
#' * `"branches"`: traced branches, one colour per branch, each labelled with
#'   its index in `x$branches`.
#' * a branch score name present in `x$score_images` (for example
#'   `"parallelism"`; see [compute_mat_properties()]): the score at each
#'   skeleton pixel on a 0-1 sequential palette. Skeleton pixels without a
#'   score (for example the ends of a branch, or a branch too short to score)
#'   are drawn in light grey, and the legend lists them as "no score".
#'
#' The function changes no persistent graphics parameter. Room for the legend
#' is reserved inside the plotting region, to the right of the image (the same
#' room [plot.mlv_drawing()] reserves, so a drawing and a MAT stage of one
#' canvas share coordinates), and later calls such as [graphics::points()] use
#' the same image coordinates.
#'
#' @param x An `mlv_mat` from [compute_mat()].
#' @param what Stage to draw; see Description.
#' @param drawing Optional `mlv_drawing` whose segments are drawn on top, for
#'   example the drawing the MAT was computed from. Its dimensions must match.
#' @param legend Draw a colour bar or key to the right of the image.
#' @param zlim Optional colour limits for `"distance"`, `"flux"` and scores.
#'   The default flux limits are symmetric about zero.
#' @param col Optional colours. For a continuous stage, a vector replacing the
#'   default palette; for `"branches"`, colours recycled over the branches in
#'   index order.
#' @param main Panel title; the default names the stage. Titles use
#'   `cex.main = 1` and `font.main = 1` unless set in `...`.
#' @param axes Draw light pixel axes.
#' @param overlay_col,overlay_lwd Colour and width of the `drawing` overlay.
#' @param discs Draw the maximal inscribed disc (radius = the branch radius)
#'   at every `discs`-th point of every branch, as faint circles, omitting any
#'   disc that would cross the edge of the region shown; `0` (the default)
#'   draws none. Each disc is centred on the axis and touches the
#'   nearest ink, so its size shows the radius that every branch score is
#'   computed from.
#' @param levels Isolines drawn over the image: `NULL` (the default) draws
#'   the threshold isoline for `"flux"` and none otherwise; a numeric vector
#'   gives the levels of the stage's quantity (distance, flux or score) to
#'   draw; `numeric(0)` draws none.
#' @param xlim,ylim Optional image region to show, in pixel coordinates (for
#'   example `xlim = c(10, 30)` to zoom in). `ylim` may be given in either
#'   order; y always increases downward.
#' @param grid Draw a faint grid on the pixel boundaries, useful when zoomed
#'   in.
#' @param ... Further graphical parameters passed to [graphics::plot.window()]
#'   and [graphics::title()], such as `cex.main`.
#' @return `x`, invisibly.
#' @seealso [plot_branch_profile()] to plot a branch quantity against image
#'   position, under the drawing.
#' @examples
#' d <- mlv_drawing(list(matrix(c(5, 10, 45, 10), 1),
#'                       matrix(c(5, 22, 45, 22), 1)), 50, 32)
#' m <- compute_mat(d)
#' plot(m, "flux", drawing = d)
#' plot(m, "distance", levels = seq(1, 10))
#' plot(m, "branches", discs = 6)
#' plot(m, "ink", xlim = c(1, 12), ylim = c(5, 15), grid = TRUE, legend = FALSE)
#' @export
plot.mlv_mat <- function(x, what = "flux", drawing = NULL, legend = TRUE,
                         zlim = NULL, col = NULL, main = NULL, axes = FALSE,
                         overlay_col = "black", overlay_lwd = 0.6, discs = 0,
                         levels = NULL, xlim = NULL, ylim = NULL, grid = FALSE,
                         ...) {
  if (!inherits(x, "mlv_mat")) stop("x must be an mlv_mat", call. = FALSE)
  stages <- c("ink", "distance", "flux", "skeleton", "branches")
  scores <- names(x$score_images)
  if (!is.character(what) || length(what) != 1L || is.na(what) ||
      !what %in% c(stages, scores))
    stop("what must be one of ", paste(c(stages, scores), collapse = ", "),
         if (!length(scores)) " (branch scores need compute_mat_properties())",
         call. = FALSE)
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend))
    stop("legend must be TRUE or FALSE", call. = FALSE)
  if (!is.null(zlim) && (!is.numeric(zlim) || length(zlim) != 2L ||
      any(!is.finite(zlim)) || zlim[2] <= zlim[1]))
    stop("zlim must be two increasing finite numbers", call. = FALSE)
  if (!is.numeric(discs) || length(discs) != 1L || !is.finite(discs) || discs < 0)
    stop("discs must be a single nonnegative number", call. = FALSE)
  if (!is.null(levels) && (!is.numeric(levels) || any(!is.finite(levels))))
    stop("levels must be NULL or finite numbers", call. = FALSE)
  if (!is.logical(grid) || length(grid) != 1L || is.na(grid))
    stop("grid must be TRUE or FALSE", call. = FALSE)
  h <- nrow(x$mask); w <- ncol(x$mask)
  view <- .mat_view(xlim, ylim, w, h)
  if (!is.null(drawing)) {
    .check_drawing(drawing)
    if (!identical(c(as.integer(drawing$height), as.integer(drawing$width)), c(h, w)))
      stop("drawing and MAT image dimensions differ", call. = FALSE)
  }
  if (is.null(main)) main <- switch(what, ink = "Ink", distance = "Distance to ink",
    flux = "Outward flux", skeleton = "Thresholded flux and skeleton",
    branches = "Traced branches", paste("Branch score:", what))

  # Work out what the stage draws before opening the plot, so the legend's
  # width can be measured and reserved.
  bar <- NULL; key <- NULL; note <- NULL; z <- NULL; limits <- NULL; colours <- NULL
  medial <- x$skeleton & is.finite(x$distance_map) & x$distance_map > 0
  unscored <- NULL
  if (what %in% c("distance", "flux") || what %in% scores) {
    sc <- .mat_scale(x, what, zlim, col)
    z <- sc$z; limits <- sc$limits; colours <- sc$colours
    if (what %in% scores) unscored <- x$skeleton & !is.finite(z)
    if (is.null(limits)) note <- "no finite values" else {
      bar <- sc[c("colours", "limits", "label", "mark")]
      if (!is.null(unscored) && any(unscored)) bar$missing <- list(label = "no score",
        fill = .mat_unscored_colour)
    }
  } else if (what == "skeleton") {
    key <- list(labels = c("ink", "above\nthreshold", "skeleton"),
                fill = c("black", .mat_above_colour, .mat_skeleton_colour))
  } else if (what == "branches") {
    key <- .mat_branch_key(length(x$branches), col)
  } else {
    key <- list(labels = c("ink", if (!is.null(drawing)) "contour"),
                fill = c("black", if (!is.null(drawing)) NA),
                line = c(NA, if (!is.null(drawing)) overlay_col))
  }

  graphics::plot.new()
  pad <- 0
  win <- function(pad) graphics::plot.window(xlim = c(view$x[1], view$x[2] + pad),
    ylim = rev(view$y), asp = 1, xaxs = "i", yaxs = "i", ...)
  win(pad)
  if (legend) {
    # The same room as plot.mlv_drawing() reserves for its colour bar, so a
    # drawing and a MAT stage of the same canvas share coordinates. It scales
    # with the image, as a replayed (narrower) copy of the figure does.
    pad <- 0.3 * max(diff(view$x), diff(view$y))
    win(pad)
  }
  # rasterImage() places the first matrix row at ytop, which is the image top
  # because the y axis is reversed; image(useRaster = TRUE) would flip it.
  draw <- function(colours) graphics::rasterImage(grDevices::as.raster(colours),
    0.5, h + 0.5, w + 0.5, 0.5, interpolate = FALSE)
  layer <- function(m, colour) if (any(m))
    draw(matrix(ifelse(m, colour, "transparent"), h, w))
  graphics::clip(view$x[1], view$x[2], view$y[1], view$y[2])

  if (what == "ink") {
    layer(x$mask, "black")
  } else if (!is.null(z)) {
    if (!is.null(limits)) draw(.mat_colour_matrix(z, colours, limits))
    if (what %in% c("distance", "flux")) layer(x$mask, "black")
    if (what %in% scores) {
      layer(x$mask, "black"); layer(unscored, .mat_unscored_colour)
    }
  } else if (what == "skeleton") {
    above <- x$AOF > x$parameters$threshold * x$parameters$samples
    layer(above & !x$mask, .mat_above_colour); layer(x$mask, "black")
    layer(x$skeleton, .mat_skeleton_colour)
  } else {
    layer(x$mask, "black")
    layer(x$skeleton & !medial, "grey50")
    pal <- if (!is.null(col)) col else .mat_branch_colours
    n <- length(x$branches)
    if (!n) note <- "no branches" else {
      if (discs > 0) .mat_draw_discs(x$branches, discs, view)
      for (i in seq_len(n)) {
        b <- x$branches[[i]]
        graphics::rect(b$x - 0.5, b$y - 0.5, b$x + 0.5, b$y + 0.5,
                       col = pal[(i - 1L) %% length(pal) + 1L], border = NA)
      }
      .mat_branch_labels(x, seq_len(n), view)
    }
  }
  if (discs > 0 && what != "branches" && length(x$branches))
    .mat_draw_discs(x$branches, discs, view)
  if (grid) .mat_pixel_grid(view)
  if (!is.null(z) && !is.null(limits)) {
    lv <- if (is.null(levels)) (if (what == "flux") bar$mark) else levels
    lv <- lv[lv > min(z, na.rm = TRUE) & lv < max(z, na.rm = TRUE)]
    if (length(lv)) {
      zz <- z; zz[!is.finite(zz)] <- NA
      graphics::contour(seq_len(w), seq_len(h), t(zz), levels = lv, add = TRUE,
        drawlabels = FALSE, lwd = if (what == "flux" && is.null(levels)) 0.9 else 0.5,
        lty = if (what == "flux" && is.null(levels)) 2 else 1,
        col = if (what == "flux" && is.null(levels)) "black" else "grey35")
    }
  }
  if (!is.null(drawing) && length(drawing$contours)) {
    s <- do.call(rbind, drawing$contours)
    graphics::segments(s[, 1], s[, 2], s[, 3], s[, 4], col = overlay_col, lwd = overlay_lwd)
  }
  usr <- graphics::par("usr")
  graphics::clip(usr[1], usr[2], min(usr[3:4]), max(usr[3:4]))
  graphics::rect(view$x[1], view$y[2], view$x[2], view$y[1], border = "grey80", lwd = 0.6)
  if (!is.null(note)) graphics::text(mean(view$x), mean(view$y), note, col = "grey30")
  if (axes) {
    graphics::axis(1, col = "grey70", col.axis = "grey40", cex.axis = 0.7)
    graphics::axis(2, col = "grey70", col.axis = "grey40", cex.axis = 0.7, las = 1)
  }
  do.call(graphics::title, utils::modifyList(
    list(main = main, cex.main = 1, font.main = 1), list(...)))
  # Drawn through recordGraphics(), so that a replay on another device size
  # (knitr's phone copy) re-measures the text and re-places the legend.
  if (legend && !is.null(bar)) .replayable(.mat_colour_bar, bar, view)
  if (legend && !is.null(key)) .replayable(.mat_key, key, view)
  invisible(x)
}

# The values, limits, colours and legend label of a continuous stage, shared
# by plot.mlv_mat() and colour_key.mlv_mat() so that a panel and its key agree.
.mat_scale <- function(x, what, zlim = NULL, col = NULL) {
  scores <- names(x$score_images)
  if (!is.character(what) || length(what) != 1L || is.na(what) ||
      !what %in% c("distance", "flux", scores))
    stop("what must be one of ", paste(c("distance", "flux", scores), collapse = ", "),
         call. = FALSE)
  z <- switch(what, distance = x$distance_map, flux = x$AOF, x$score_images[[what]])
  limits <- if (!is.null(zlim)) zlim else if (what %in% scores) c(0, 1) else .mat_range(z)
  if (what == "flux" && is.null(zlim) && !is.null(limits)) limits <- c(-1, 1) * max(abs(limits))
  list(z = z, limits = limits, colours = if (!is.null(col)) col else .mat_palette(what),
       label = switch(what, distance = "distance\n(px)",
                      flux = paste0("flux, sum\nof ", x$parameters$samples),
                      paste0(what, "\n(0", .range_dash(), "1)")),
       mark = if (what == "flux" && !is.null(x$parameters$threshold))
         x$parameters$threshold * x$parameters$samples)
}

#' @rdname colour_key
#' @param what For a medial-axis representation: the continuous stage
#'   (`"distance"`, `"flux"` or a branch score), as in [plot.mlv_mat()]. The
#'   flux key marks the skeleton threshold.
#' @param col For a medial-axis representation: colours replacing the default
#'   palette, as in [plot.mlv_mat()].
#' @export
colour_key.mlv_mat <- function(x, what = "flux", zlim = NULL, col = NULL, ...) {
  if (!is.null(zlim) && (!is.numeric(zlim) || length(zlim) != 2L ||
      any(!is.finite(zlim)) || zlim[2] <= zlim[1]))
    stop("zlim must be two increasing finite numbers", call. = FALSE)
  sc <- .mat_scale(x, what, zlim, col)
  if (is.null(sc$limits)) stop("the stage has no finite values", call. = FALSE)
  colour_key.default(sc$colours, sc$limits, sub("\n", " ", sc$label, fixed = TRUE),
                     mark = sc$mark, mark_label = if (!is.null(sc$mark)) "threshold", ...)
}

# Colours used by more than one stage. Above-threshold and skeleton pixels
# take light and dark reds, the colours of positive flux in the flux stage.
.mat_above_colour <- "#F2C4BA"
.mat_skeleton_colour <- "#A1181D"
.mat_unscored_colour <- "grey75"

# Okabe-Ito colours (without black), recycled over branches in index order.
# Blue and orange come late: blue means negative flux in the stages before,
# and the vignettes reserve blue and orange for default and published values.
.mat_branch_colours <- c("#CC79A7", "#009E73", "#D55E00", "#56B4E9", "#E69F00",
                         "#0072B2", "#F0E442")

# Default 101-colour palettes. hcl.colors() can return NA for a colour whose
# HCL coordinates fall outside sRGB (in R 4.5.1 the neutral midpoint of
# "Blue-Red 3" is NA), and an NA colour is drawn transparent;
# .no_na_colours() (shared with the drawing keys) interpolates such entries.
.mat_palette <- function(what, n = 101L) {
  spec <- switch(what, distance = list("Purples 3", TRUE), flux = list("Blue-Red 3", FALSE),
                 list("Viridis", FALSE))
  .no_na_colours(grDevices::hcl.colors(n, spec[[1]], rev = spec[[2]]))
}

# Colour index of each value; values outside the limits take the end colours
# and missing values are transparent.
.mat_colour_matrix <- function(z, colours, limits) {
  k <- 1L + floor((pmin(pmax(z, limits[1]), limits[2]) - limits[1]) /
                    diff(limits) * (length(colours) - 1L) + 0.5)
  k[!is.finite(z)] <- NA
  matrix(ifelse(is.na(k), "transparent", colours[k]), nrow(z), ncol(z))
}

.mat_range <- function(z) {
  if (!any(is.finite(z))) return(NULL)
  r <- range(z[is.finite(z)])
  if (r[1] == r[2]) r <- r + c(-0.5, 0.5)
  r
}

# The image region shown: pixel-edge limits in x and y (y given top first).
.mat_view <- function(xlim, ylim, w, h) {
  check <- function(v, name) if (!is.null(v) && (!is.numeric(v) || length(v) != 2L ||
    any(!is.finite(v)) || v[1] == v[2])) stop(name, " must be two different finite numbers",
    call. = FALSE)
  check(xlim, "xlim"); check(ylim, "ylim")
  list(x = if (is.null(xlim)) c(0.5, w + 0.5) else sort(xlim),
       y = if (is.null(ylim)) c(0.5, h + 0.5) else sort(ylim))
}

# Legend for the branches stage: one swatch per branch in its own colour, up
# to six; beyond that the first five and a count of the rest.
.mat_branch_key <- function(n, col = NULL) {
  if (!n) return(list(labels = "ink", fill = "black"))
  pal <- if (is.null(col)) .mat_branch_colours else col
  shown <- if (n <= 6L) n else 5L
  fill <- pal[(seq_len(shown) - 1L) %% length(pal) + 1L]
  labels <- paste("branch", seq_len(shown))
  if (shown < n) {
    labels <- c(labels, sprintf("+ %d more", n - shown))
    fill <- c(fill, NA)
  }
  list(labels = c("ink", labels), fill = c("black", fill))
}

# Maximal inscribed discs (outlines) at every k-th point of each branch; with a
# view, only discs that lie wholly inside it, so none is cut by the frame.
.mat_draw_discs <- function(branches, every, view = NULL, col = "grey35") {
  every <- max(1L, as.integer(round(every)))
  t <- seq(0, 2 * pi, length.out = 73L)
  drawn <- vector("list", length(branches))
  for (k in seq_along(branches)) {
    b <- branches[[k]]
    n <- nrow(b); if (!n) next
    start <- ceiling(every / 2)
    i <- if (start > n) ceiling(n / 2) else seq(start, n, by = every)
    i <- i[i >= 1 & i <= n & is.finite(b$radius[i]) & b$radius[i] > 0]
    if (!is.null(view)) i <- i[b$x[i] - b$radius[i] >= view$x[1] &
      b$x[i] + b$radius[i] <= view$x[2] & b$y[i] - b$radius[i] >= view$y[1] &
      b$y[i] + b$radius[i] <= view$y[2]]
    for (j in i) graphics::polygon(b$x[j] + b$radius[j] * cos(t),
      b$y[j] + b$radius[j] * sin(t), border = grDevices::adjustcolor(col, 0.55),
      col = NA, lwd = 0.6)
    if (length(i)) graphics::points(b$x[i], b$y[i], pch = 16, cex = 0.25,
                                    col = grDevices::adjustcolor(col, 0.8))
    drawn[[k]] <- as.integer(i)
  }
  invisible(drawn)
}

# Branch index labels with a white halo, each placed at the candidate spot
# around its branch's midpoint that lies farthest from ink, skeleton and the
# labels already placed.
.mat_branch_labels <- function(x, index, view, cex = 0.75, col = "black") {
  occupied <- which(x$mask | x$skeleton, arr.ind = TRUE)
  occ <- cbind(occupied[, 2], occupied[, 1])
  placed <- matrix(numeric(), 0, 2)
  d <- 1.3 * abs(graphics::strheight("8", cex = cex))
  angles <- seq(0, 2 * pi, length.out = 13L)[-13L]
  for (i in index) {
    b <- x$branches[[i]]; if (!nrow(b)) next
    k <- ceiling(nrow(b) / 2)
    cand <- rbind(cbind(b$x[k] + d * cos(angles), b$y[k] + d * sin(angles)),
                  cbind(b$x[k] + 1.8 * d * cos(angles), b$y[k] + 1.8 * d * sin(angles)))
    inside <- cand[, 1] > view$x[1] + d & cand[, 1] < view$x[2] - d &
      cand[, 2] > view$y[1] + d & cand[, 2] < view$y[2] - d
    if (any(inside)) cand <- cand[inside, , drop = FALSE]
    dist <- function(q, p) sqrt((q[, 1] - p[1])^2 + (q[, 2] - p[2])^2)
    centre <- matrix(c(b$x[k], b$y[k]), 1)
    # Clearance from ink, skeleton and earlier labels; ties go to the nearer ring.
    clearance <- apply(cand, 1, function(p) min(dist(occ, p),
      if (nrow(placed)) dist(placed, p) - d else Inf) - 0.05 * dist(centre, p))
    p <- cand[which.max(clearance), ]
    .halo_text(p[1], p[2], i, cex = cex, font = 2,
               col = if (is.function(col)) col(i) else col)
    placed <- rbind(placed, p)
  }
}

.halo_text <- function(x, y, labels, cex = 0.75, col = "black", halo = "white", ...) {
  dx <- 0.06 * graphics::strwidth("M", cex = cex)
  dy <- 0.06 * graphics::strheight("M", cex = cex)
  for (a in seq(0, 2 * pi, length.out = 9L)[-9L])
    graphics::text(x + dx * cos(a), y + dy * sin(a), labels, cex = cex, col = halo, xpd = NA, ...)
  graphics::text(x, y, labels, cex = cex, col = col, xpd = NA, ...)
}

# Faint lines on the pixel boundaries inside the view.
.mat_pixel_grid <- function(view) {
  gx <- seq(ceiling(view$x[1] - 0.5), floor(view$x[2] - 0.5)) + 0.5
  gy <- seq(ceiling(view$y[1] - 0.5), floor(view$y[2] - 0.5)) + 0.5
  graphics::abline(v = gx, h = gy, col = grDevices::adjustcolor("grey60", 0.5), lwd = 0.4)
}

.mat_legend_cex <- 0.7

.mat_bar_labels <- function(bar) {
  at <- pretty(bar$limits, 4); at <- at[at >= bar$limits[1] & at <= bar$limits[2]]
  labels <- format(at, trim = TRUE)
  if (!is.null(bar$mark) && bar$mark >= bar$limits[1] && bar$mark <= bar$limits[2]) {
    keep <- abs(at - bar$mark) > 0.12 * diff(bar$limits)
    at <- c(at[keep], bar$mark)
    labels <- c(labels[keep], sprintf("%.1f", bar$mark))
  }
  list(at = at, labels = labels)
}

# Converts inches to user units on the current plot (x and y separately).
.inch_x <- function(v) v * diff(graphics::par("usr")[1:2]) / graphics::par("pin")[1]
.inch_y <- function(v) v * abs(diff(graphics::par("usr")[3:4])) / graphics::par("pin")[2]

# A vertical colour bar right of the image: the title (quantity and unit)
# above the bar, ticks and labels to its right, all within the reserved pad.
.mat_colour_bar <- function(bar, view) {
  cex <- .mat_legend_cex
  x0 <- view$x[2] + .inch_x(0.08); x1 <- x0 + .inch_x(0.1)
  title <- strsplit(bar$label, "\n")[[1]]
  line_h <- .inch_y(1.25 * graphics::strheight("M", "inches", cex = cex))
  y_top <- view$y[1] + 1.3 * length(title) * line_h + .inch_y(0.06)
  y_bottom <- view$y[2]
  if (!is.null(bar$missing)) y_bottom <- y_bottom - 1.6 * line_h
  if (y_bottom - y_top < .inch_y(0.3)) y_top <- y_bottom - .inch_y(0.3)
  graphics::text(x0, view$y[1], bar$label, adj = c(0, 1), cex = cex, xpd = NA)
  n <- length(bar$colours)
  graphics::rasterImage(grDevices::as.raster(matrix(rev(bar$colours), ncol = 1L)),
                        x0, y_bottom, x1, y_top, interpolate = FALSE, xpd = NA)
  graphics::rect(x0, y_bottom, x1, y_top, border = "grey40", lwd = 0.5, xpd = NA)
  pos <- function(v) y_bottom - (v - bar$limits[1]) / diff(bar$limits) * (y_bottom - y_top)
  tl <- .mat_bar_labels(bar)
  if (!is.null(bar$mark) && bar$mark >= bar$limits[1] && bar$mark <= bar$limits[2])
    graphics::segments(x0, pos(bar$mark), x1, pos(bar$mark), lwd = 1, lty = 2, xpd = NA)
  graphics::segments(x1, pos(tl$at), x1 + .inch_x(0.03), pos(tl$at), col = "grey40",
                     lwd = 0.5, xpd = NA)
  graphics::text(x1 + .inch_x(0.05), pos(tl$at), tl$labels, adj = c(0, 0.5), cex = cex, xpd = NA)
  if (!is.null(bar$missing)) {
    ym <- view$y[2] - 0.6 * line_h
    graphics::rect(x0, ym - .inch_y(0.065), x1, ym + .inch_y(0.065),
                   col = bar$missing$fill, border = NA, xpd = NA)
    graphics::text(x1 + .inch_x(0.05), ym, bar$missing$label, adj = c(0, 0.5), cex = cex,
                   xpd = NA)
  }
}

# A key of filled squares (or short lines) and labels right of the image.
.mat_key <- function(key, view) {
  cex <- .mat_legend_cex
  x0 <- view$x[2] + .inch_x(0.08); x1 <- x0 + .inch_x(0.1)
  line <- if (is.null(key$line)) rep(NA, length(key$labels)) else key$line
  step <- .inch_y(1.7 * graphics::strheight("M", "inches", cex = cex))
  n <- lengths(strsplit(key$labels, "\n"))
  y <- view$y[1] + step * (cumsum(n) - n / 2)
  for (i in seq_along(key$labels)) {
    if (!is.na(line[i])) graphics::segments(x0, y[i], x1, y[i], col = line[i], lwd = 1.5, xpd = NA)
    else if (!is.na(key$fill[i])) graphics::rect(x0, y[i] - .inch_y(0.05), x1,
      y[i] + .inch_y(0.05), col = key$fill[i], border = NA, xpd = NA)
    graphics::text(x1 + .inch_x(0.04), y[i], key$labels[i], adj = c(0, 0.5), cex = cex,
                   xpd = NA)
  }
}
