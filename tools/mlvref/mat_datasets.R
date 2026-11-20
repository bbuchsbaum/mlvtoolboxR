# Track B: end-to-end per-pixel MAT scores, R vs stored MATLAB, on the exported datasets.
# Usage (repository root):
#   MLVREF_LIB=<lib> Rscript tools/mlvref/mat_datasets.R <per_dataset> <cores> [seed]
# Reads .agent-work/mlvref/<dataset>/{segments,drawings,mat_samples}.csv.gz (local natural data)
# and writes per-drawing aggregate rows to .agent-work/track-b/e2e/<dataset>.rds. Only
# aggregates leave the machine (see mat_datasets_summary.R).
source("tools/mlvref/mat_common.R")
args <- commandArgs(TRUE)
per <- if (length(args) >= 1) as.integer(args[1]) else 10L
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
seed <- if (length(args) >= 3) as.integer(args[3]) else 20261119L
src <- ".agent-work/mlvref"; out <- Sys.getenv("E2E_OUT", ".agent-work/track-b/e2e"); dir.create(out, FALSE, TRUE)
THR_ANGLE <- asin(0.3 * pi / 2) * 180 / pi   # computeAllMATfromVecLD uses the omitted-threshold 0.3 (MLV-M004)
all_configs <- list(
  compat = list(morphology = "hmt8", aof = "mlv2023", branch_rule = "mlv2023", emulate = FALSE),
  default = list(morphology = "minimal", aof = "centered", branch_rule = "invariant", emulate = FALSE),
  default_hmt8 = list(morphology = "hmt8", aof = "centered", branch_rule = "invariant", emulate = FALSE),
  pratt = list(morphology = "pratt", aof = "mlv2023", branch_rule = "invariant", emulate = FALSE),
  hmt8 = list(morphology = "hmt8", aof = "mlv2023", branch_rule = "invariant", emulate = FALSE),
  hmt8_mlvmap = list(morphology = "hmt8", aof = "mlv2023", branch_rule = "invariant", emulate = TRUE),
  hmt8_centered = list(morphology = "hmt8", aof = "centered", branch_rule = "invariant", emulate = FALSE))
pick <- Sys.getenv("E2E_CONFIGS"); configs <- if (nzchar(pick)) all_configs[strsplit(pick, ",")[[1]]] else all_configs
props <- c("parallelism", "separation", "mirror")

one_drawing <- function(k, seg, dims, samp) {
  s <- seg[seg$drawing == k, ]; s <- s[(s$x2 - s$x1)^2 + (s$y2 - s$y1)^2 > 0, ]
  cont <- lapply(split(s, s$contour), function(z) as.matrix(z[order(z$segment), c("x1", "y1", "x2", "y2")]))
  dm <- dims[dims$drawing == k, ]
  x <- mlv_drawing(unname(cont), dm$width, dm$height)
  sm <- samp[samp$drawing == k, ]
  rows <- list(); contour_rows <- list()
  for (cn in names(configs)) { cf <- configs[[cn]]; t0 <- proc.time()[["elapsed"]]
    m <- compute_mat(x, threshold_angle = THR_ANGLE, morphology = cf$morphology, aof = cf$aof, branch_rule = cf$branch_rule)
    m <- tryCatch(compute_mat_properties(m, props), error = function(e) e)
    if (inherits(m, "error")) { rows[[length(rows) + 1]] <- data.frame(drawing = k, config = cn, property = props, n = NA,
      covered = NA, within_1e4 = NA, within_1e2 = NA, sum_abs = NA, q50 = NA, q90 = NA, mean_matlab = NA, mean_r_same = NA,
      mean_r_all = NA, n_r_all = NA, skel = NA, branches = NA, seconds = NA, error = conditionMessage(m)); next }
    if (cf$emulate) { rc <- unique(sm[, c("y", "x")]); v <- map_at(m, as.matrix(rc), matlab_image = TRUE)
      imgs <- lapply(props, function(p) { im <- matrix(NA_real_, x$height, x$width); im[as.matrix(rc)] <- v[[p]]; im }); names(imgs) <- props
    } else { y <- map_mat_properties(x, m); imgs <- y$mat_images }
    el <- proc.time()[["elapsed"]] - t0
    for (p in props) { st <- sm[sm$property == p, ]; if (!nrow(st)) next
      # Per-contour means: MATLAB's contour mean is the mean nonzero mapped score over the
      # contour's insertShape raster, reproduced by the package rasterizer (MLV-R010).
      mimg <- matrix(NA_real_, x$height, x$width); mimg[cbind(st$y, st$x)] <- st$score
      cm <- t(vapply(seq_along(x$contours), function(i) { ix <- ns$.raster_indices(x, i)
        a <- mimg[ix]; b <- imgs[[p]][ix]; c(if (any(is.finite(a))) mean(a, na.rm = TRUE) else NA, if (any(is.finite(b))) mean(b, na.rm = TRUE) else NA) }, double(2)))
      contour_rows[[length(contour_rows) + 1]] <- data.frame(drawing = k, config = cn, property = p, matlab = cm[, 1], r = cm[, 2])
      v <- imgs[[p]][cbind(st$y, st$x)]; e <- abs(v - st$score)
      rows[[length(rows) + 1]] <- data.frame(drawing = k, config = cn, property = p, n = nrow(st),
        covered = sum(is.finite(v)), within_1e4 = sum(e <= 1e-4, na.rm = TRUE), within_1e2 = sum(e <= 1e-2, na.rm = TRUE),
        sum_abs = sum(e, na.rm = TRUE), q50 = median(e, na.rm = TRUE), q90 = unname(quantile(e, .9, na.rm = TRUE)),
        mean_matlab = mean(st$score), mean_r_same = mean(v, na.rm = TRUE),
        mean_r_all = mean(imgs[[p]][is.finite(imgs[[p]])]), n_r_all = sum(is.finite(imgs[[p]])),
        skel = sum(m$skeleton), branches = length(m$branches), seconds = el, error = NA_character_) } }
  out <- do.call(rbind, rows); attr(out, "contours") <- do.call(rbind, contour_rows); out
}
sets <- sort(basename(dirname(Sys.glob(file.path(src, "*", "COMPLETE")))))
sets <- if (nzchar(Sys.getenv("E2E_SETS"))) strsplit(Sys.getenv("E2E_SETS"), ",")[[1]] else sets; sets <- sets[file.exists(file.path(src, sets, "mat_samples.csv.gz"))]
for (ds in sets) { fo <- file.path(out, paste0(ds, ".rds")); if (file.exists(fo)) { cat(ds, "done\n"); next }
  seg <- read.csv(file.path(src, ds, "segments.csv.gz"))[, c("drawing", "contour", "segment", "x1", "y1", "x2", "y2")]
  dims <- read.csv(file.path(src, ds, "drawings.csv.gz"))
  samp <- read.csv(file.path(src, ds, "mat_samples.csv.gz"))
  set.seed(seed); ids <- sort(sample(unique(samp$drawing), min(per, length(unique(samp$drawing)))))
  only <- Sys.getenv("E2E_DRAWINGS"); if (nzchar(only)) ids <- as.integer(strsplit(only, ",")[[1]])
  t0 <- Sys.time()
  res <- parallel::mclapply(ids, function(k) tryCatch(one_drawing(k, seg, dims, samp),
    error = function(e) { message(ds, " drawing ", k, ": ", conditionMessage(e)); NULL }), mc.cores = cores)
  bad <- ids[vapply(res, function(r) is.null(r) || !is.data.frame(r), TRUE)]
  ok <- Filter(is.data.frame, res)
  res <- do.call(rbind, ok); attr(res, "failed_drawings") <- bad
  attr(res, "contours") <- do.call(rbind, lapply(ok, attr, "contours"))
  saveRDS(res, fo); cat(ds, length(ids), "drawings", format(Sys.time() - t0), "\n")
}
