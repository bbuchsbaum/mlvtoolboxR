# Report T/X events in whole-photo coordinates for each model view.
library(mlvtoolbox)
args <- commandArgs(TRUE)
folder <- args[1L]
panel <- read.csv(file.path(folder, "panel.csv"), colClasses = c(id = "character"))
trace_segments <- function(segments, width, height) {
  coords <- as.matrix(segments[, c("x1", "y1", "x2", "y2")])
  lines <- lapply(seq_len(nrow(coords)), function(i) coords[i, , drop = FALSE])
  d <- mlv_drawing(lines, width, height)
  j <- detect_junctions(trace_contours(rasterize(d), simplify = 1))
  j[j$type %in% c("T", "X"), c("x", "y", "type")]
}

one_photo <- function(i) {
  photo <- panel[i, ]
  read_view <- function(suffix) read.csv(file.path(folder, paste0(photo$id, suffix)))
  tag <- function(j, method) {
    j$id <- rep(photo$id, nrow(j))
    j$method <- rep(method, nrow(j))
    j
  }
  detect_crop <- function(s) {
    left <- s$left[1L]
    top <- s$top[1L]
    width <- s$width[1L]
    height <- s$height[1L]
    j <- trace_segments(s, width, height)
    # Exclude a four-pixel margin only on artificial cuts through the photo.
    inside <- (left == 0 | j$x > 4) & (top == 0 | j$y > 4) &
      (left + width == photo$width | j$x <= width - 4) &
      (top + height == photo$height | j$y <= height - 4)
    j <- j[inside, ]
    j$x <- j$x + left
    j$y <- j$y + top
    tag(j, paste("crop", left, top, sep = "_"))
  }
  s <- read_view(".tiles.csv")
  views <- lapply(split(s, list(s$left, s$top), drop = TRUE), detect_crop)
  if ("--include-full" %in% args) {
    suffixes <- c(mlsd_mask = ".segments.csv", mlsd_clahe = ".clahe.segments.csv")
    for (method in names(suffixes)) {
      full <- trace_segments(read_view(suffixes[[method]]), photo$width, photo$height)
      views[[method]] <- tag(full, method)
    }
  }
  do.call(rbind, views)
}

empty <- data.frame(x = numeric(), y = numeric(), type = character(), id = character(),
                    method = character())
events <- do.call(rbind, c(list(empty), lapply(seq_len(nrow(panel)), one_photo)))
output <- if ("--include-full" %in% args) "all-views.csv" else "crop-junctions.csv"
write.csv(events, file.path(folder, output), row.names = FALSE)
