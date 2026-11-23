# Called by evaluate_hawp.py with an explicit photo/threshold plan.
library(mlvtoolbox)
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1L]
plan <- read.csv(args[2L], colClasses = c(id = "character"))
detect_photo <- function(id) {
  rgb <- read_mlv_image(file.path(folder, paste0(id, ".small.png")))
  predictions <- read.csv(file.path(folder, paste0(id, ".hawp.csv")))
  at_cutoff <- function(cutoff) {
    coords <- as.matrix(predictions[predictions$score > cutoff, c("x1", "y1", "x2", "y2")])
    native <- mlv_drawing(lapply(seq_len(nrow(coords)), function(i) coords[i, , drop = FALSE]),
                          dim(rgb)[2L], dim(rgb)[1L])
    provider <- function(image) {
      stopifnot(identical(image, rgb))
      rasterize(native)
    }
    spec <- feature_spec(extraction = list(provider = "hawp", provider_version = "hawpv2-92ae446",
                                           simplify = 1))
    traced <- trace_rgb(rgb, backend = provider, spec = spec)
    events <- function(drawing, method) {
      j <- detect_junctions(drawing)[, c("x", "y", "type")]
      j$id <- rep(id, nrow(j))
      j$method <- rep(sprintf("hawp_%s_%.2f", method, cutoff), nrow(j))
      j
    }
    rbind(events(native, "vector"), events(traced, "mask"))
  }
  do.call(rbind, lapply(plan$cutoff[plan$id == id], at_cutoff))
}
write.csv(do.call(rbind, lapply(unique(plan$id), detect_photo)), args[3L], row.names = FALSE)
