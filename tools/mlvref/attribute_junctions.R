# Attribute R-vs-MATLAB junction differences to individual mechanisms by
# walking from a faithful emulation of the pinned MLV pipeline (V0) to the
# mlvtoolbox definition one switch at a time:
#   V0 MLV faithful (validated against stored outputs)
#   V1 + interior-arm orientation fixed           (MLV-M007)
#   V2 + transitive merging                       (MLV-M006)
#   V3 + negative-coordinate intersections kept   (MLV-M008)
#   V4 + no 1e-4 endpoint snapping
#   V5 + no short-contour exclusion
#   V6 + zero-length segments removed             (MLV-M015)
#   V7 + neighbourhood-exit arms                  (MLV-R005/R006 definition)
#   V8 + Euclidean absolute gap                   (MLV-M018; package default)
#   R  mlvtoolbox::detect_junctions defaults (rule = "exit", gap_metric = "euclidean")
# V8>R checks that the emulated definition explains R completely; the V6>V7
# differences are classified by arm-rule features (arm_features).
# Side metric V6>V6euclid: effect of a Euclidean absolute gap (MLV-M018).
# Usage: Rscript tools/mlvref/attribute_junctions.R [lib] stored [cores] [tag]
args <- commandArgs(trailingOnly = TRUE)
src <- readLines("tools/mlvref/compare_junctions.R"); n <- grep("^sets <- complete_sets", src)
eval(parse(text = src[1:(n - 1)]))
source("tools/mlvref/mlv_junction_emulator.R")
tag <- if (length(args) >= 4) args[4] else ""
detail_root <- paste0(".agent-work/track-a/detail/attribution", tag)
dir.create(detail_root, recursive = TRUE, showWarnings = FALSE)

variants <- list(
  V0 = list(),
  V1 = list(angles = "fixed"),
  V2 = list(angles = "fixed", merge = "transitive"),
  V3 = list(angles = "fixed", merge = "transitive", negative = FALSE),
  V4 = list(angles = "fixed", merge = "transitive", negative = FALSE, snap = FALSE),
  V5 = list(angles = "fixed", merge = "transitive", negative = FALSE, snap = FALSE, short = FALSE))

exact_count <- function(A, B, k) {
  cmp <- compare_one(A, B, k); m <- cmp$matched
  ex <- if (is.null(m)) logical() else m$pos_ok & m$type_a == m$type_b & m$contours_ok & m$angles_ok & m$minmax_ok
  list(nA = nrow(A), nB = nrow(B), exact = sum(ex), cmp = cmp, ex = ex)
}

remap <- function(E, b) {
  if (!nrow(E)) return(E)
  cs <- split_ids(E$contours); ss <- split_ids(E$segments)
  E$contours <- vapply(cs, function(v) paste(b$contour_map[v], collapse = ";"), "")
  E$segments <- vapply(seq_along(cs), function(i) paste(mapply(function(ci, si) b$map[[ci]][si],
    cs[[i]], ss[[i]]), collapse = ";"), "")
  E
}

# Features of an R junction that make the mlvtoolbox arm rule differ from
# the (M007-corrected) MLV rule; r = endpoint radius.
arm_features <- function(d, p, members, r = 2) {
  walk <- offset <- both_near <- FALSE
  for (i in seq_len(nrow(members))) {
    z <- d$contours[[members$contour_index[i]]]; si <- members$segment_index[i]; n <- nrow(z)
    dd <- function(v) sum((v - p)^2)
    da <- dd(z[si, 1:2]); db <- dd(z[si, 3:4])
    if (da < r^2 && db < r^2) both_near <- TRUE
    if (da >= r^2 && db >= r^2) {
      u <- z[si, 3:4] - z[si, 1:2]; w <- p - z[si, 1:2]
      if (abs(u[1] * w[2] - u[2] * w[1]) / sqrt(sum(u^2)) > 1e-9) offset <- TRUE
    } else {
      # arm segments examined by MLV: the incident one and its neighbour across
      # the near vertex; R walks further whenever their far vertex is inside r.
      if (da < db) {
        if (si < n && db <= r^2) walk <- TRUE
        if (si > 2 && dd(z[si - 1, 1:2]) <= r^2) walk <- TRUE
      } else {
        if (si > 1 && da <= r^2) walk <- TRUE
        if (si < n - 1 && dd(z[si + 1, 3:4]) <= r^2) walk <- TRUE
      }
    }
  }
  c(walk = walk, interior_offset = offset, both_near = both_near)
}

one_drawing <- function(s, w, h, A, k) {
  s <- s[order(s$contour, s$segment), ]
  cont <- unname(lapply(split(s, s$contour), function(p) as.matrix(p[, c("x1", "y1", "x2", "y2")])))
  V <- lapply(variants, function(v) do.call(mlv_emulate, c(list(contours = cont), v)))
  b <- build_drawing(s, w, h)
  dc <- unname(b$drawing$contours)
  V$V6 <- remap(do.call(mlv_emulate, c(list(contours = dc), variants$V5)), b)
  V$V7 <- remap(do.call(mlv_emulate, c(list(contours = dc), variants$V5, list(arms = "exit"))), b)
  V$V8 <- remap(do.call(mlv_emulate, c(list(contours = dc), variants$V5, list(arms = "exit", gap = "euclidean"))), b)
  V6e <- remap(do.call(mlv_emulate, c(list(contours = dc), variants$V5, list(gap = "euclidean"))), b)
  R <- with_minmax(r_junctions(b)); raw <- attr(R, "raw")
  chain <- c(list(S = with_minmax(A)), V, list(R = R))
  steps <- lapply(seq_len(length(chain) - 1L), function(i) exact_count(chain[[i]], chain[[i + 1L]], k))
  names(steps) <- paste(names(chain)[-length(chain)], names(chain)[-1L], sep = ">")
  sr <- exact_count(chain$S, R, k); ge <- exact_count(V$V6, V6e, k)
  last <- steps[["V6>V7"]]; m <- last$cmp$matched
  feats <- NULL
  r_v7 <- V$V7
  if (!is.null(m) && any(!last$ex)) {
    bad <- which(!last$ex)
    feats <- do.call(rbind, lapply(bad, function(t) {
      jj <- m$cand[t]
      rj <- match_points(cbind(r_v7$x[jj], r_v7$y[jj]), cbind(raw$x, raw$y))$j
      if (!length(rj)) return(NULL); jj <- rj
      f <- arm_features(b$drawing, c(raw$x[jj], raw$y[jj]), raw$members[[jj]])
      data.frame(drawing = k, t(f), arms_differ = m$n_arms_a[t] != m$n_arms_b[t],
                 type_differ = m$type_a[t] != m$type_b[t], pos_ok = m$pos_ok[t], contours_ok = m$contours_ok[t])
    }))
  }
  list(counts = data.frame(drawing = k, step = c(names(steps), "S>R", "V6>V6euclid"),
         nA = c(vapply(steps, `[[`, 0, "nA"), sr$nA, ge$nA), nB = c(vapply(steps, `[[`, 0, "nB"), sr$nB, ge$nB),
         exact = c(vapply(steps, `[[`, 0, "exact"), sr$exact, ge$exact)),
       residual = list(V6 = V$V6, V7 = V$V7, R = R, cmp = last$cmp, ex = last$ex, features = feats,
                       ref_only = NROW(last$cmp$ref_only), cand_only = NROW(last$cmp$cand_only)))
}

sets <- complete_sets(export_root)
all_counts <- list()
for (ds in sets) {
  t0 <- Sys.time()
  segs <- read_gz(file.path(export_root, ds, "segments.csv"))
  draws <- read_gz(file.path(export_root, ds, "drawings.csv"))
  sj <- read_gz(file.path(export_root, ds, "junctions.csv"))
  out <- parallel::mclapply(draws$drawing, function(k) {
    i <- which(draws$drawing == k)
    one_drawing(segs[segs$drawing == k, ], draws$width[i], draws$height[i], sj[sj$drawing == k, ], k)
  }, mc.cores = cores)
  bad <- vapply(out, inherits, NA, "try-error")
  if (any(bad)) stop(ds, ": ", as.character(out[bad][[1]]))
  saveRDS(out, file.path(detail_root, paste0(ds, ".rds")))
  cnt <- do.call(rbind, lapply(out, `[[`, "counts")); cnt$dataset <- ds
  all_counts[[ds]] <- cnt
  agg <- aggregate(cbind(nA, nB, exact) ~ step, cnt, sum)
  message(ds, " (", round(as.numeric(Sys.time() - t0, units = "secs")), "s): ",
          paste(sprintf("%s %d/%d", agg$step, agg$exact, agg$nA), collapse = ", "))
}
cnt <- do.call(rbind, all_counts)
agg <- aggregate(cbind(nA, nB, exact) ~ dataset + step, cnt, sum)
utils::write.csv(agg, file.path(receipt_root, sprintf("attribution-by-dataset%s.csv", tag)), row.names = FALSE)
tot <- aggregate(cbind(nA, nB, exact) ~ step, cnt, sum)
utils::write.csv(tot, file.path(receipt_root, sprintf("attribution-total%s.csv", tag)), row.names = FALSE)
print(tot)
feat <- do.call(rbind, lapply(sets, function(ds) {
  out <- readRDS(file.path(detail_root, paste0(ds, ".rds")))
  f <- do.call(rbind, lapply(out, function(o) o$residual$features))
  if (is.null(f)) return(NULL); f$dataset <- ds; f
}))
if (!is.null(feat)) {
  fc <- aggregate(drawing ~ walk + interior_offset + both_near + arms_differ + type_differ + pos_ok + contours_ok,
                  feat, length)
  names(fc)[names(fc) == "drawing"] <- "junctions"
  utils::write.csv(fc, file.path(receipt_root, sprintf("residual-features%s.csv", tag)), row.names = FALSE)
  print(fc)
}
