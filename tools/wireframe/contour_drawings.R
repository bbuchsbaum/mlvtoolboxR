# Trace predicted segments and score the junctions formed by the R drawing.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Expected the photo folder, output folder, and roster")
pkgload::load_all(".", quiet = TRUE)
rows <- read.csv(args[3], colClasses = c(id = "character"))
photos <- read.csv(file.path(args[1], "panel.csv"), colClasses = c(id = "character"))
for (i in seq_len(nrow(rows))) {
  id <- rows$id[i]
  method <- rows$method[i]
  prefix <- file.path(args[2], paste(id, method, sep = "."))
  size <- photos[match(id, photos$id), c("width", "height")]
  if (method == "artist") {
    refs <- jsonlite::fromJSON(file.path(args[1], paste0(id, ".json")), simplifyVector = FALSE)
    contours <- lapply(refs$polylines, function(path) {
      xy <- matrix(unlist(path$points), ncol = 2, byrow = TRUE)
      cbind(head(xy, -1), tail(xy, -1))
    })
    drawing <- mlv_drawing(contours, size$width, size$height, zero_length = "drop")
  } else if (method %in% c("union_vector", "curve_rank")) {
    seg <- read.csv(paste0(prefix, ".segments.csv"))
    indices <- split(seq_len(nrow(seg)), seg$group)
    paths <- lapply(indices, function(ix) as.matrix(seg[ix, c("x1","y1","x2","y2"), drop = FALSE]))
    drawing <- mlv_drawing(unname(paths), size$width, size$height, zero_length = "drop")
  } else if (method == "fdog") {
    rgb <- read_mlv_image(file.path(args[1], paste0(id, ".small.png")))
    drawing <- trace_contours(xdog_lines(rgb), simplify = 1)
  } else {
    seg <- as.matrix(read.csv(paste0(prefix, ".segments.csv")))
    paths <- lapply(seq_len(nrow(seg)), function(k) seg[k, , drop = FALSE])
    native <- mlv_drawing(paths, size$width, size$height, zero_length = "drop")
    drawing <- trace_contours(rasterize(native), simplify = 1)
  }
  png::writePNG(1.0 * rasterize(drawing), paste0(prefix, ".drawing.png"))
  if (method != "artist") {
    # Export only predicted T/X geometry. Artist contours already live in the
    # staged source JSON and are not duplicated here.
    sites <- detect_junctions(drawing)
    sites <- sites[sites$type %in% c("T", "X"), ]
    export_contours <- lapply(seq_along(drawing$contours), function(k) list(
      id=as.character(drawing$ids[k]),
      segments=unname(split(as.matrix(drawing$contours[[k]]), row(as.matrix(drawing$contours[[k]]))))
    ))
    export_sites <- lapply(seq_len(nrow(sites)), function(k) list(
      junction_id=sites$junction_id[k], x=sites$x[k], y=sites$y[k],
      type=as.character(sites$type[k]), arms=unname(sites$arms[[k]]),
      members=lapply(seq_len(nrow(sites$members[[k]])), function(m) list(
        contour_id=as.character(sites$members[[k]]$contour_id[m]),
        contour_index=sites$members[[k]]$contour_index[m],
        segment_index=sites$members[[k]]$segment_index[m]
      ))
    ))
    jsonlite::write_json(list(schema_version=1, contours=export_contours,
                              sites=export_sites), paste0(prefix, ".contours.json"),
                         auto_unbox=TRUE, digits=10)
    table <- sites[, c("x", "y", "type")]
    # R uses mathematical orientation, the benchmark vectors have image y down.
    table$directions <- vapply(sites$arms, function(a) paste(-a, collapse = " "), "")
    write.csv(table, paste0(prefix, ".sites.csv"), row.names = FALSE)
  }
}
