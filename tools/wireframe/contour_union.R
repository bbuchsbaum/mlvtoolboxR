# Development-only follow-up. Keep complete observed curves within the ink
# budget of the package FDoG drawing for the same photograph.
args <- commandArgs(TRUE)
if (!(length(args) %in% c(2L, 3L)))
  stop("Expected staged Toronto assets, output drawings, and optional roster")
pkgload::load_all(".", quiet = TRUE)
roster <- read.csv(file.path(args[1], "panel.csv"), colClasses = c(id = "character"))
optional_roster <- length(args) == 3L
if (optional_roster) {
  requested <- read.csv(args[3], colClasses = c(id = "character"))
  if (!("id" %in% names(requested)) || anyDuplicated(requested$id))
    stop("optional roster must have distinct id values")
  if (any(!requested$id %in% roster$id))
    stop("optional roster contains IDs absent from staged inputs")
  rows <- roster[match(requested$id, roster$id), , drop = FALSE]
  for (field in intersect(c("category", "panel", "width", "height", "pixel_sha256"), names(requested))) {
    expected <- as.character(requested[[field]])
    # The evaluation roster calls this frozen panel ``reserve-test`` while
    # staged photos and geometry retain Toronto's source label ``reserve``.
    if (field == "panel") expected[expected == "reserve-test"] <- "reserve"
    if (!identical(as.character(rows[[field]]), expected))
      stop("optional roster disagrees with staged inputs for ", field)
  }
} else {
  rows <- roster[roster$panel == "development", ]
}
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
    if (!optional_roster && source == "fdog") {
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
  if (!optional_roster) {
    write.csv(segments, paste0(prefix, ".union.segments.csv"), row.names = FALSE)
  }
  groups <- rep(seq_along(selected), vapply(selected, nrow, integer(1)))
  write.csv(data.frame(group = groups, segments), paste0(prefix, ".union_vector.segments.csv"),
            row.names = FALSE)
}
plan <- tempfile(fileext = ".csv")
methods <- if (optional_roster) "union_vector" else c("union", "union_vector", "fdog_retrace")
write.csv(expand.grid(id = rows$id, method = methods), plan, row.names = FALSE)
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--vanilla", "tools/wireframe/contour_drawings.R",
                    shQuote(args[1]), shQuote(args[2]), shQuote(plan)))
unlink(plan)
if (status != 0L) stop("Tracing the development contour union failed")
