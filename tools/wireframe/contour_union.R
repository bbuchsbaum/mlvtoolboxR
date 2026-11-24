# Development-only follow-up. Keep complete observed curves within the ink
# budget of the package FDoG drawing for the same photograph.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Expected staged Toronto assets and output drawings")
pkgload::load_all(".", quiet = TRUE)
roster <- read.csv(file.path(args[1], "panel.csv"), colClasses = c(id = "character"))
rows <- roster[roster$panel == "development", ]
read_paths <- function(path) {
  obj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  lapply(obj$contours, function(x) matrix(unlist(x$segments), ncol = 4, byrow = TRUE))
}
for (i in seq_len(nrow(rows))) {
  id <- rows$id[i]
  w <- rows$width[i]
  h <- rows$height[i]
  prefix <- file.path(args[2], id)
  budget <- sum(png::readPNG(paste0(prefix, ".fdog.drawing.png")) > 0.5)
  ink <- matrix(FALSE, h, w)
  remaining <- budget
  selected <- list()
  for (source in c("mlsd", "fdog")) {
    paths <- read_paths(paste0(prefix, ".", source, ".contours.json"))
    if (source == "fdog") {
      control <- if (length(paths)) do.call(rbind, paths) else matrix(numeric(), 0L, 4L)
      colnames(control) <- c("x1", "y1", "x2", "y2")
      write.csv(control, paste0(prefix, ".fdog_retrace.segments.csv"), row.names = FALSE)
    }
    if (!length(paths)) next
    drawing <- mlv_drawing(paths, w, h)
    hits <- lapply(seq_along(paths), function(j) .raster_indices(drawing, j))
    # Favor long strokes covering previously undrawn regions. Nearness
    # affects ordering only. Each accepted curve retains every segment.
    near <- distance_transform(ink)$distance <= 2
    utility <- vapply(hits, function(p) sum(!near[p]), integer(1))
    for (j in order(-utility, seq_along(paths))) {
      new <- sum(!ink[hits[[j]]])
      if (new > 0L && new <= remaining) {
        ink[hits[[j]]] <- TRUE
        remaining <- remaining - new
        selected[[length(selected) + 1L]] <- paths[[j]]
      }
    }
  }
  stopifnot(sum(ink) <= budget)
  segments <- if (length(selected)) do.call(rbind, selected) else matrix(numeric(), 0L, 4L)
  colnames(segments) <- c("x1", "y1", "x2", "y2")
  write.csv(segments, paste0(prefix, ".union.segments.csv"), row.names = FALSE)
  groups <- rep(seq_along(selected), vapply(selected, nrow, integer(1)))
  write.csv(data.frame(group = groups, segments), paste0(prefix, ".union_vector.segments.csv"),
            row.names = FALSE)
}
plan <- tempfile(fileext = ".csv")
write.csv(expand.grid(id = rows$id, method = c("union", "union_vector", "fdog_retrace")), plan, row.names = FALSE)
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--vanilla", "tools/wireframe/contour_drawings.R",
                    shQuote(args[1]), shQuote(args[2]), shQuote(plan)))
unlink(plan)
if (status != 0L) stop("Tracing the development contour union failed")
