# Trace the screened sources and export per-curve working-frame pixel pools.
#
# Usage, from the package root:
#   Rscript --vanilla tools/wireframe/source_pool/trace_sources.R \
#     .agent-work/toronto/staged .agent-work/contour-study/draw-toronto \
#     .agent-work/source-pool docs/evidence/wireframe/curve-rank/development.csv
#
# Inputs in the source-pool folder come from prepare_photos.py. For each
# development photo this writes pools/<id>.json with zero-based column-major
# raster indices per curve for fdog, mlsd (as in .agent-work/contour-proposal)
# and the two new sources, plus each new curve's ordering key, and
# draw/<id>.<source>.{contours.json,drawing.png} in the working frame.
args <- commandArgs(TRUE)
if (length(args) != 4L) stop("Expected staged, draw-toronto, source-pool and roster paths")
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
staged <- args[1]
base <- args[2]
out <- args[3]
roster <- read.csv(args[4], colClasses = c(id = "character"))
stopifnot(nrow(roster) == 60L, all(roster$panel == "development"))
panel <- read.csv(file.path(staged, "panel.csv"), colClasses = c(id = "character"))
dir.create(file.path(out, "pools"), showWarnings = FALSE)
dir.create(file.path(out, "draw"), showWarnings = FALSE)

W <- 256L
H <- 192L

# The staging transform for one-based pixel centres: (p - 0.5) * w / W + 0.5.
to_working <- function(z, from_w, from_h) {
  sx <- W / from_w
  sy <- H / from_h
  cbind((z[, 1] - 0.5) * sx + 0.5, (z[, 2] - 0.5) * sy + 0.5,
        (z[, 3] - 0.5) * sx + 0.5, (z[, 4] - 0.5) * sy + 0.5)
}
pool_of <- function(drawing) {
  lapply(seq_along(drawing$contours), function(j) .raster_indices(drawing, j) - 1L)
}
read_segments <- function(file) {
  obj <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  lapply(obj$contours, function(x) matrix(unlist(x$segments), ncol = 4, byrow = TRUE))
}
export <- function(drawing, prefix) {
  png::writePNG(1.0 * rasterize(drawing), paste0(prefix, ".drawing.png"))
  contours <- lapply(seq_along(drawing$contours), function(k) list(
    id = as.character(drawing$ids[k]),
    segments = unname(split(drawing$contours[[k]], row(drawing$contours[[k]])))
  ))
  jsonlite::write_json(list(schema_version = 1, contours = contours, sites = list()),
                       paste0(prefix, ".contours.json"), auto_unbox = TRUE, digits = 10)
}
polyline_length <- function(z) sum(sqrt((z[, 3] - z[, 1])^2 + (z[, 4] - z[, 2])^2))

log <- list()
for (i in seq_len(nrow(roster))) {
  id <- roster$id[i]
  size <- panel[match(id, panel$id), c("width", "height")]
  stopifnot(size$width == W, size$height == H)
  pools <- list()

  # Existing sources, exactly as .agent-work/contour-proposal/paths.R.
  for (method in c("fdog", "mlsd")) {
    seg <- read_segments(file.path(base, paste(id, method, "contours.json", sep = ".")))
    pools[[method]] <- pool_of(mlv_drawing(seg, W, H))
  }

  # Check that the FDoG call reproduces the existing source drawing at 256 px.
  small <- read_mlv_image(file.path(staged, paste0(id, ".small.png")))
  fdog256 <- trace_contours(xdog_lines(small), simplify = 1)
  old <- png::readPNG(file.path(base, paste0(id, ".fdog.drawing.png"))) > 0.5
  fdog256_same <- identical(rasterize(fdog256), old)

  # FDoG on the 512 x 384 LANCZOS resize, mapped to the working frame.
  big <- read_mlv_image(file.path(out, "fdog512", paste0(id, ".512.png")))
  stopifnot(dim(big)[1] == 384L, dim(big)[2] == 512L)
  d512 <- trace_contours(xdog_lines(big), simplify = 1)
  f512 <- mlv_drawing(lapply(d512$contours, to_working, 512, 384), W, H, ids = d512$ids)
  pools$fdog512 <- pool_of(f512)
  pools$fdog512_length <- vapply(f512$contours, polyline_length, numeric(1))
  export(f512, file.path(out, "draw", paste(id, "fdog512", sep = ".")))

  # Structured edges: thresholded NMS mask at 800 x 600, traced, then mapped.
  mask <- png::readPNG(file.path(out, "sed", paste0(id, ".mask.png")))
  if (length(dim(mask)) == 3L) mask <- mask[, , 1]
  mask <- mask > 0.5
  stopifnot(nrow(mask) == 600L, ncol(mask) == 800L)
  strength <- matrix(readBin(file.path(out, "sed", paste0(id, ".strength.f32")), "double",
                             n = 800L * 600L, size = 4L, endian = "little"),
                     nrow = 600L, ncol = 800L, byrow = TRUE)
  stopifnot(identical(mask, strength >= 0.1))
  raw <- trace_contours(mask)
  sed <- simplify_contours(raw, 1)
  direct <- trace_contours(mask, simplify = 1)
  stopifnot(identical(sed$contours, direct$contours), identical(sed$ids, direct$ids))
  # Strength: mean NMS edge strength over the traced pixel chain of each curve.
  chain_strength <- vapply(raw$contours, function(z) {
    xy <- rbind(z[, 1:2, drop = FALSE], z[nrow(z), 3:4])
    mean(strength[unique((xy[, 1] - 1) * 600 + xy[, 2])])
  }, numeric(1))
  names(chain_strength) <- raw$ids
  sedw <- mlv_drawing(lapply(direct$contours, to_working, 800, 600), W, H, ids = direct$ids)
  stopifnot(identical(sedw$ids, direct$ids))
  pools$sed <- pool_of(sedw)
  pools$sed_strength <- unname(chain_strength[sedw$ids])
  export(sedw, file.path(out, "draw", paste(id, "sed", sep = ".")))

  jsonlite::write_json(pools, file.path(out, "pools", paste0(id, ".json")),
                       auto_unbox = FALSE, digits = NA)
  log[[i]] <- data.frame(
    id = id, fdog256_reproduced = fdog256_same,
    fdog_curves = length(pools$fdog), mlsd_curves = length(pools$mlsd),
    fdog512_curves = length(pools$fdog512), sed_curves = length(pools$sed),
    sed_mask_pixels = sum(mask), sed_min_chain_strength = min(chain_strength)
  )
  cat(id, "fdog256", fdog256_same, "fdog512", length(pools$fdog512), "sed",
      length(pools$sed), "\n")
}
log <- do.call(rbind, log)
write.csv(log, file.path(out, "trace-log.csv"), row.names = FALSE)
cat("FDoG at 256 reproduced on", sum(log$fdog256_reproduced), "of", nrow(log), "photos\n")
