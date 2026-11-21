#' Plot branch quantities under the drawing they come from
#'
#' Draws each medial-axis representation as a column: the ink and the chosen
#' branches on top, and below it one panel per quantity, plotted against the
#' image x coordinate on the same horizontal scale as the image. A point on a
#' branch and its score therefore lie on the same vertical line, so a reader
#' can look straight down from a feature of the drawing to its effect on the
#' score. Positions on the y axis are easier to compare than colours, which
#' makes small differences between branches readable.
#'
#' The figure fills the current device page, like [graphics::pairs()]. Panel
#' heights follow the image's aspect ratio, and text size follows the device
#' width (smaller on narrow devices). The graphics state is restored on exit.
#'
#' @param x An `mlv_mat` scored by [compute_mat_properties()], or a list of
#'   them (one column each; list names become column titles).
#' @param what Branch quantities to plot, one row each: `"radius"` or any
#'   score present in the branches, for example `c("parallelism",
#'   "separation")`.
#' @param branches Which branches to show: `NULL` (all), `"longest"`, an
#'   integer vector of branch indices used for every column, or a list with
#'   one such entry per column.
#' @param discs Draw the maximal inscribed disc at every `discs`-th point of
#'   the shown branches (see [plot.mlv_mat()]); `0` draws none.
#' @param col Colours of the profile curves, recycled over the columns; the
#'   default is black. Where a column shows several branches, they are told
#'   apart by line type and by their index at the line ends. The axis in the
#'   drawing row is always drawn in the medial-axis colour (the skeleton red
#'   of `plot(x, "skeleton")`), half a pixel wide, whatever `col` is.
#' @param ylim Optional y limits: one numeric pair for every row, or a named
#'   list of pairs keyed by `what`. The default is `c(0, 1)` for scores and
#'   from zero to the largest value for `"radius"`.
#' @param band Optional x range (two numbers, in pixels) shaded in every
#'   panel, drawing included, for example the stretch over which a summary is
#'   computed.
#' @param context Draw the other columns' curves faintly behind each profile,
#'   so that every column can be compared with the rest: `TRUE` for every row,
#'   `FALSE` for none, or the names of the rows to give context.
#' @param stack Place the columns one below another instead of side by side,
#'   each at the full device width (useful on narrow devices or when a drawing
#'   needs labels).
#' @param blank Names in `what` that are not branch quantities but empty rows,
#'   drawn with only their row label on the shared x scale, for `panel` to
#'   fill (for example a diagram of fitted spans). Their y range is 0 to 1.
#' @param labels Label each line with its branch index at its right end when
#'   a column shows more than one branch.
#' @param main Column titles; the default is `names(x)`.
#' @param panel Optional function called after each panel is drawn, with the
#'   panel's coordinate system active, as `panel(column, row)`: `column` is
#'   the column index and `row` is `"drawing"` or one of `what`. Use it to
#'   add annotations.
#' @param cex Text size; the default scales with the device width.
#' @param lwd Width of the profile lines.
#' @return Invisibly, a data frame with one row per plotted point: `column`,
#'   `branch`, `what`, `x`, `y` (image coordinates of the branch point) and
#'   `value`.
#' @seealso [plot.mlv_mat()], [compute_mat_properties()].
#' @examples
#' d <- mlv_drawing(list(matrix(c(5, 8, 75, 8), 1),
#'                       matrix(c(5, 30, 75, 22), 1)), 80, 36)
#' m <- compute_mat_properties(compute_mat(d), c("parallelism", "separation"))
#' plot_branch_profile(m, c("radius", "parallelism"), branches = "longest",
#'                     discs = 8)
#' @export
plot_branch_profile <- function(x, what = "parallelism", branches = NULL, discs = 0,
                                col = NULL, ylim = NULL, band = NULL, labels = TRUE,
                                main = NULL, panel = NULL, cex = NULL, lwd = 1.6,
                                context = FALSE, blank = character(), stack = FALSE) {
  mats <- if (inherits(x, "mlv_mat")) list(x) else x
  if (!is.list(mats) || !length(mats) || !all(vapply(mats, inherits, logical(1), "mlv_mat")))
    stop("x must be an mlv_mat or a list of them", call. = FALSE)
  nc <- length(mats)
  if (!is.character(what) || !length(what) || anyNA(what))
    stop("what must name branch quantities", call. = FALSE)
  if (!is.character(blank)) stop("blank must be a character vector", call. = FALSE)
  for (m in mats) {
    have <- unique(unlist(lapply(m$branches, names)))
    miss <- setdiff(what, c(have, blank))
    if (length(m$branches) && length(miss))
      stop("branches have no column ", paste(miss, collapse = ", "),
           " (scores need compute_mat_properties())", call. = FALSE)
  }
  if (!is.numeric(discs) || length(discs) != 1L || !is.finite(discs) || discs < 0)
    stop("discs must be a single nonnegative number", call. = FALSE)
  if (!is.null(band) && (!is.numeric(band) || length(band) != 2L || any(!is.finite(band))))
    stop("band must be two finite numbers", call. = FALSE)
  if (!is.null(panel) && !is.function(panel)) stop("panel must be a function", call. = FALSE)
  if (is.logical(context) && length(context) == 1L && !is.na(context)) {
    context <- if (context) what else character()
  } else if (!is.character(context)) stop("context must be TRUE, FALSE or row names", call. = FALSE)
  pick <- if (is.list(branches) && !is.character(branches)) branches else rep(list(branches), nc)
  if (length(pick) != nc) stop("branches must have one entry per column", call. = FALSE)
  shown <- lapply(seq_len(nc), function(k) .profile_branches(mats[[k]], pick[[k]]))
  if (is.null(main)) main <- names(mats)
  if (is.null(main)) main <- rep("", nc)
  main <- rep_len(as.character(main), nc)

  # Data: one row per plotted point.
  rows <- list()
  for (k in seq_len(nc)) for (i in shown[[k]]) for (v in setdiff(what, blank)) {
    b <- mats[[k]]$branches[[i]]
    if (nrow(b)) rows[[length(rows) + 1L]] <- data.frame(column = k, branch = i, what = v,
      x = b$x, y = b$y, value = as.numeric(b[[v]]), stringsAsFactors = FALSE)
  }
  data <- if (length(rows)) do.call(rbind, rows) else data.frame(column = integer(),
    branch = integer(), what = character(), x = numeric(), y = numeric(), value = numeric())
  limits <- lapply(what, function(v) {
    if (is.list(ylim) && !is.null(ylim[[v]])) return(ylim[[v]])
    if (is.numeric(ylim) && length(ylim) == 2L) return(ylim)
    if (v == "radius") {
      r <- data$value[data$what == v & is.finite(data$value)]
      return(c(0, if (length(r)) max(pretty(c(0, max(r)))) else 1))
    }
    c(0, 1)
  })
  names(limits) <- what

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  din <- graphics::par("din")
  s <- if (is.null(cex)) min(1, max(0.72, din[1] / 7)) else cex
  csi <- 0.2 * s                          # one text line, in inches
  left <- 2.3 * csi; right <- 0.25 * csi
  label_h <- 1.25 * csi                   # room for the row label above each profile
  title_h <- if (any(nzchar(main))) 1.5 * csi else 0.1
  if (!is.logical(stack) || length(stack) != 1L || is.na(stack))
    stop("stack must be TRUE or FALSE", call. = FALSE)
  blocks <- if (stack) 1L else nc          # columns side by side
  page_h <- if (stack) din[2] / nc else din[2]
  pw <- din[1] / blocks - left - right
  # Each drawing is cropped vertically to its ink, its shown branches and their
  # discs; the full width is kept so that x matches the profiles below.
  # Every column gets the same height, centred on its content, so the drawings
  # of a row share their top and bottom.
  yview <- lapply(seq_len(nc), function(k) .profile_yview(mats[[k]], shown[[k]]))
  tall <- max(vapply(yview, diff, numeric(1)))
  yview <- lapply(seq_len(nc), function(k) {
    h <- nrow(mats[[k]]$mask); v <- mean(yview[[k]]) + c(-0.5, 0.5) * tall
    v <- v - min(0, v[1] - 0.5) - max(0, v[2] - h - 0.5)
    c(max(0.5, v[1]), min(h + 0.5, v[2]))
  })
  aspect <- max(vapply(seq_len(nc), function(k) diff(yview[[k]]) / ncol(mats[[k]]$mask),
                       numeric(1)))
  draw_h <- pw * aspect + title_h + 0.04
  bottom <- 2.6 * csi
  avail <- page_h - draw_h - bottom
  # Blank rows (annotations) get half the height of a data row.
  weight <- ifelse(what %in% blank, 0.5, 1)
  row_h <- avail / length(what)          # includes label_h
  # Too little room: shrink the drawing. The floor leaves room for each row's
  # margins also when knitr replays the figure on a narrower (shorter) device.
  if (row_h < 0.75) {
    row_h <- min(0.75, (page_h - bottom) / (length(what) + 1))
    draw_h <- page_h - bottom - row_h * length(what)
  }
  heights <- c(draw_h, row_h * length(what) * weight / sum(weight))
  heights[length(heights)] <- heights[length(heights)] + bottom
  # knitr's phone copy replays the figure at about half the height while the
  # margins keep their size in inches: give every row at least twice its
  # margins, taking the room from the drawing.
  need <- 2 * c(title_h + 0.04, rep(label_h + 0.04, length(what))) + 0.1
  need[length(need)] <- need[length(need)] + 2 * bottom
  rows <- pmax(heights[-1], need[-1])
  heights <- c(max(need[1], page_h - sum(rows)), rows)
  cells <- seq_len((length(what) + 1L) * nc)
  if (stack) graphics::layout(matrix(cells, ncol = 1L), heights = rep(heights, nc))
  else graphics::layout(matrix(cells, ncol = nc), heights = heights)
  graphics::par(cex = s, mgp = c(1.9, 0.45, 0), tcl = -0.25, las = 1, xpd = FALSE)
  pal <- function(k, i) if (is.null(col)) "black" else rep_len(col, nc)[k]
  lty_of <- function(k, i) if (length(shown[[k]]) > 1L) match(i, shown[[k]]) else 1L

  for (k in seq_len(nc)) {
    m <- mats[[k]]
    graphics::par(mai = c(0.04, left, title_h, right))
    plot.mlv_mat(m, "ink", legend = FALSE, main = main[k], cex.main = 1, ylim = yview[[k]])
    view <- .mat_view(NULL, yview[[k]], ncol(m$mask), nrow(m$mask))
    if (!is.null(band)) graphics::rect(min(band), view$y[1], max(band), view$y[2],
      col = grDevices::adjustcolor("grey50", 0.12), border = NA)
    if (discs > 0) .mat_draw_discs(m$branches[shown[[k]]], discs, view)
    .mat_draw_axis(m$branches[shown[[k]]])
    if (labels && length(shown[[k]]) > 1L)
      .mat_branch_labels(m, shown[[k]], view, col = .mat_axis_col)
    if (!is.null(panel)) panel(k, "drawing")
    xr <- graphics::par("usr")[1:2]
    for (j in seq_along(what)) {
      v <- what[j]; last <- j == length(what)
      graphics::par(mai = c(if (last) bottom + 0.04 else 0.04, left, label_h, right))
      graphics::plot.new()
      graphics::plot.window(xlim = xr, ylim = limits[[v]], xaxs = "i")
      usr <- graphics::par("usr")
      if (!is.null(band)) graphics::rect(min(band), usr[3], max(band), usr[4],
                                         col = "grey93", border = NA)
      if (k == 1L || stack) graphics::mtext(.profile_label(v), side = 3, line = 0.15, adj = 0,
                                   at = xr[1], cex = 0.8 * s, col = "grey20")
      if (v %in% blank) {
        if (!is.null(panel)) panel(k, v)
        next
      }
      graphics::axis(2, col = "grey60", col.axis = "grey25", lwd = 0.6, cex.axis = 0.8,
                     at = .profile_ticks(v, limits[[v]]), labels = FALSE)
      if (k == 1L || nc == 1L || stack) {
        # Score axes are labelled at their ends only; the middle tick marks
        # the midpoint, so a zoomed range (0 to 0.5) reads as one.
        at <- .profile_ticks(v, limits[[v]])
        lab <- if (v == "radius") at else at[c(1, length(at))]
        graphics::axis(2, at = lab, labels = format(lab, drop0trailing = TRUE),
                       lwd = 0, col.axis = "grey25", cex.axis = 0.8)
      }
      if (last) {
        # A range frame over the whole image width, so no data lie past its end.
        ix <- pretty(c(1, ncol(m$mask)), 4)
        ix <- unique(c(1, ix[ix > 1 & ix <= ncol(m$mask)]))
        graphics::axis(1, at = ix, col = "grey60", col.axis = "grey25", lwd = 0,
                       lwd.ticks = 0.6, cex.axis = 0.8)
        graphics::segments(1, usr[3], ncol(m$mask), usr[3], col = "grey60", lwd = 0.6, xpd = NA)
        graphics::mtext("x (px)", side = 1, line = 1.6, cex = 0.85 * s)
      }
      graphics::segments(xr[1], usr[3], xr[1], usr[4], col = "grey60", lwd = 0.6)
      if (v %in% context) for (kk in setdiff(seq_len(nc), k)) for (i in shown[[kk]]) {
        d <- data[data$column == kk & data$branch == i & data$what == v, ]
        graphics::lines(d$x, d$value, col = "grey80", lwd = 1)
      }
      for (i in shown[[k]]) {
        d <- data[data$column == k & data$branch == i & data$what == v, ]
        if (!nrow(d)) next
        graphics::lines(d$x, d$value, col = pal(k, i), lwd = lwd, lty = lty_of(k, i))
        # Points with no finite neighbour, or a branch whose finite values all
        # lie at one x, would draw no visible line: mark them as dots.
        fin <- is.finite(d$value)
        one <- fin & !(c(FALSE, utils::head(fin, -1)) | c(fin[-1], FALSE))
        if (any(fin) && diff(range(d$x[fin])) == 0) one <- fin
        if (any(one)) graphics::points(d$x[one], d$value[one], pch = 16, cex = 0.5,
                                       col = pal(k, i))
        if (labels && length(shown[[k]]) > 1L && any(is.finite(d$value))) {
          e <- which(is.finite(d$value)); e <- e[which.max(d$x[e])]
          graphics::text(d$x[e], d$value[e], i, pos = 4, offset = 0.25, cex = 0.75,
                         font = 2, col = pal(k, i), xpd = NA)
        }
      }
      if (!is.null(panel)) panel(k, v)
    }
  }
  invisible(data)
}

.profile_branches <- function(m, which) {
  n <- length(m$branches)
  if (!n) return(integer())
  if (is.null(which)) return(seq_len(n))
  if (identical(which, "longest")) return(which.max(vapply(m$branches, nrow, integer(1))))
  if (!is.numeric(which) || any(!is.finite(which)) || any(which < 1 | which > n))
    stop("branches must be NULL, \"longest\" or indices between 1 and ", n, call. = FALSE)
  as.integer(which)
}

# Vertical extent (pixel edges, top first) of the ink and the shown branches,
# with a small margin, within the image.
.profile_yview <- function(m, shown) {
  h <- nrow(m$mask)
  y <- which(rowSums(m$mask) > 0)
  for (i in shown) y <- c(y, m$branches[[i]]$y)
  if (!length(y)) return(c(0.5, h + 0.5))
  c(max(0.5, min(y) - 2.5), min(h + 0.5, max(y) + 2.5))
}

# Score axes are ticked at both ends and the middle, so a zoomed range (for
# example 0 to 0.5) reads as one.
.profile_ticks <- function(v, lim) if (v == "radius") pretty(lim, 3)[pretty(lim, 3) <= lim[2]] else
  c(lim[1], mean(lim), lim[2])

.profile_label <- function(v) if (v == "radius") "radius (px)" else v
