# Track B: aggregate end-to-end receipts (per dataset x config x property) from mat_datasets.R output.
indir <- Sys.getenv("E2E_OUT", ".agent-work/track-b/e2e2"); files <- Sys.glob(file.path(indir, "*.rds"))
r <- do.call(rbind, lapply(files, function(f) { z <- readRDS(f); z$dataset <- sub("\\.rds$", "", basename(f))
  if (is.null(z$error)) z$error <- NA_character_; z }))
err <- r[!is.na(r$error), ]; r <- r[is.na(r$error), ]
if (nrow(err)) { cat("config failures (drawings):\n"); print(table(unique(err[, c("dataset", "drawing", "config")])$config)) }
r$family <- sub("_[0-9].*$", "", sub("_vecLD$", "", r$dataset))
agg <- function(g) {
  s <- split(r, g)
  do.call(rbind, lapply(s, function(z) data.frame(z[1, intersect(c("family", "dataset", "config", "property"), names(z)), drop = FALSE],
    drawings = length(unique(paste(z$dataset, z$drawing))), pixels = sum(z$n), covered = sum(z$covered) / sum(z$n),
    within_1e4 = sum(z$within_1e4) / sum(z$n), within_1e2 = sum(z$within_1e2) / sum(z$n),
    mean_abs = { f <- is.finite(z$sum_abs); sum(z$sum_abs[f]) / sum(z$covered[f]) },  # one IAPS drawing stores 11 separation = Inf (radius 0) median_drawing_q50 = median(z$q50), median_drawing_q90 = median(z$q90),
    drawing_mean_cor = { f <- is.finite(z$mean_matlab) & is.finite(z$mean_r_same); cor(z$mean_matlab[f], z$mean_r_same[f]) },
    drawing_mean_cor_own_pixels = { f <- is.finite(z$mean_matlab) & is.finite(z$mean_r_all); cor(z$mean_matlab[f], z$mean_r_all[f]) },
    sec_per_drawing = mean(z$seconds), row.names = NULL)))
}
by_family <- agg(list(r$family, r$config, r$property))
overall <- agg(list(r$config, r$property)); overall$family <- "ALL"
out <- rbind(by_family[, names(overall)], overall)
out <- out[order(out$family, out$property, out$config), ]
write.csv(out, "docs/evidence/mlvref/mat/e2e_summary.csv", row.names = FALSE)
# Per-contour means (MATLAB: mean nonzero mapped score over the contour raster).
cc <- do.call(rbind, lapply(files, function(f) { z <- attr(readRDS(f), "contours"); if (is.null(z)) return(NULL)
  z$dataset <- sub("\\.rds$", "", basename(f)); z }))
if (!is.null(cc)) {
  cc$family <- sub("_[0-9].*$", "", sub("_vecLD$", "", cc$dataset))
  csum <- function(z) { z$matlab[!is.finite(z$matlab)] <- NA; ok <- is.finite(z$matlab) & is.finite(z$r); e <- abs(z$matlab[ok] - z$r[ok])
    data.frame(contours = nrow(z), both_defined = sum(ok), matlab_only = sum(is.finite(z$matlab) & !is.finite(z$r)),
      r_only = sum(!is.finite(z$matlab) & is.finite(z$r)), cor = cor(z$matlab[ok], z$r[ok]),
      median_abs = median(e), q90_abs = unname(quantile(e, .9)), within_0.01 = mean(e <= .01), within_0.05 = mean(e <= .05),
      mean_bias = mean(z$r[ok] - z$matlab[ok])) }
  g <- split(cc, list(cc$config, cc$property), drop = TRUE)
  ctab <- do.call(rbind, lapply(names(g), function(n) cbind(family = "ALL", config = g[[n]]$config[1], property = g[[n]]$property[1], csum(g[[n]]))))
  g2 <- split(cc, list(cc$family, cc$config, cc$property), drop = TRUE)
  ctab <- rbind(ctab, do.call(rbind, lapply(names(g2), function(n) cbind(family = g2[[n]]$family[1], config = g2[[n]]$config[1], property = g2[[n]]$property[1], csum(g2[[n]])))))
  write.csv(ctab, "docs/evidence/mlvref/mat/e2e_contour_summary.csv", row.names = FALSE)
  print(ctab[ctab$family == "ALL", ], digits = 3, row.names = FALSE)
}
print(overall[, c("config", "property", "drawings", "covered", "within_1e4", "within_1e2", "mean_abs", "drawing_mean_cor", "sec_per_drawing")], digits = 3, row.names = FALSE)
cat("failed drawings:", sum(vapply(files, function(f) length(attr(readRDS(f), "failed_drawings")), 1L)), "\n")
