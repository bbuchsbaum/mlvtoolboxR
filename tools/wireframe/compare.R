# Run from the repository root with the package built at the revision under test.
library(mlvtoolbox)
folder <- commandArgs(TRUE)[1L]
panel <- read.csv(file.path(folder, "panel.csv"), colClasses = c(id = "character"))
model_version <- "mlsd-512-large-453cafa"

compare_photo <- function(id) {
  path <- function(suffix) file.path(folder, paste0(id, suffix))
  rgb <- read_mlv_image(path(".small.png"))
  witness <- id %in% c("00079423", "00078318")
  record <- function(drawing, method) {
    junctions <- detect_junctions(drawing)[, c("x", "y", "type")]
    if (witness) {
      write.csv(as.data.frame(drawing)[, c("x1", "y1", "x2", "y2")],
                path(paste0(".", method, ".traces.csv")), row.names = FALSE)
    }
    junctions$id <- rep(id, nrow(junctions))
    junctions$method <- rep(method, nrow(junctions))
    junctions
  }

  results <- lapply(c(0, 6), function(gap) {
    spec <- feature_spec(extraction = list(lines = list(junction_gap = gap), simplify = 1))
    drawing <- trace_rgb(rgb, backend = xdog_lines, spec = spec)
    method <- paste0("fdog", gap)
    if (witness) {
      stages <- line_extraction_stages(rgb, junction_gap = gap)
      png::writePNG(stages$mask * 1, path(paste0(".", method, ".mask.png")))
      stopifnot(identical(drawing$contours, trace_contours(stages$mask, simplify = 1)$contours))
    }
    record(drawing, method)
  })

  read_drawing <- function(suffix) {
    segments <- as.matrix(read.csv(path(suffix)))
    lines <- lapply(seq_len(nrow(segments)), function(i) segments[i, , drop = FALSE])
    mlv_drawing(lines, dim(rgb)[2L], dim(rgb)[1L])
  }

  trace_model <- function(native, method, version, parameters) {
    # Each prediction file comes from prepare.py processing this exact RGB image.
    mask <- rasterize(native)
    attr(mask, "extraction") <- list(method = version, parameters = parameters)
    provider <- function(image) {
      stopifnot(identical(image, rgb))
      mask
    }
    spec <- feature_spec(extraction = list(provider = "mlsd", provider_version = version, simplify = 1))
    if (witness) png::writePNG(mask * 1, path(paste0(".", method, ".mask.png")))
    record(trace_rgb(rgb, backend = provider, spec = spec), method)
  }

  native <- read_drawing(".segments.csv")
  enhanced <- read_drawing(".clahe.segments.csv")
  annotated <- read_drawing(".annot.csv")
  model_parameters <- list(score = 0.10, distance = 20)
  clahe_parameters <- c(model_parameters, list(clip_limit = 2, tile_grid = c(8L, 8L)))
  do.call(rbind, c(results, list(
    trace_model(native, "mlsd_mask", model_version, model_parameters),
    record(native, "mlsd_vector"),
    trace_model(enhanced, "mlsd_clahe", paste0(model_version, "-clahe-lab"), clahe_parameters),
    record(trace_contours(rasterize(annotated), simplify = 1), "annotation_mask"),
    record(annotated, "annotation_vector")
  )))
}
write.csv(do.call(rbind, lapply(panel$id, compare_photo)), file.path(folder, "junctions.csv"),
          row.names = FALSE)
