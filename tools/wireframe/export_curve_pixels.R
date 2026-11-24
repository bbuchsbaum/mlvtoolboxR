# Export one zero-based, column-major pixel path for every source contour.
# Usage: Rscript export_curve_pixels.R DRAWINGS ROSTER OUTPUT

read_contour_paths <- function(path) {
  obj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (is.null(obj$contours)) stop("contours missing from ", path, call. = FALSE)
  lapply(obj$contours, function(x) matrix(unlist(x$segments), ncol = 4L, byrow = TRUE))
}

export_curve_pixels <- function(drawings, roster, output) {
  roster <- read.csv(roster, colClasses = c(id = "character"))
  if (!("id" %in% names(roster))) stop("ROSTER must contain id", call. = FALSE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(roster))) {
    id <- roster$id[i]
    png_paths <- file.path(drawings,
                           paste(id, c("fdog", "mlsd"), "drawing.png", sep = "."))
    dimensions <- lapply(png_paths, function(path) {
      if (!file.exists(path)) stop("source drawing PNG missing: ", path, call. = FALSE)
      dim(png::readPNG(path))[1:2]
    })
    if (!identical(dimensions[[1]], dimensions[[2]]))
      stop("FDoG and M-LSD source PNG dimensions disagree for ", id, call. = FALSE)
    height <- dimensions[[1]][1]
    width <- dimensions[[1]][2]
    result <- list()
    for (source in c("fdog", "mlsd")) {
      path <- file.path(drawings, paste(id, source, "contours.json", sep = "."))
      curves <- read_contour_paths(path)
      drawing <- mlv_drawing(curves, width, height)
      hits <- lapply(seq_along(curves), function(j) .raster_indices(drawing, j))
      # Prove the exported per-curve pixels reconstruct the complete source
      # raster.  R's matrix indexing is the required Fortran/column-major order.
      reconstructed <- matrix(FALSE, height, width)
      if (length(hits)) reconstructed[unlist(hits, use.names = FALSE)] <- TRUE
      source_mask <- png::readPNG(png_paths[[match(source, c("fdog", "mlsd"))]]) > 0.5
      if (!identical(reconstructed, source_mask))
        stop("per-curve pixels do not reconstruct the full ", source,
             " source mask for ", id, call. = FALSE)
      result[[source]] <- lapply(hits, function(x) as.integer(x - 1L))
    }
    jsonlite::write_json(result, file.path(output, paste0(id, ".json")),
                         auto_unbox = FALSE)
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 3L) stop("Expected DRAWINGS ROSTER OUTPUT", call. = FALSE)
  suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
  export_curve_pixels(args[1], args[2], args[3])
}
