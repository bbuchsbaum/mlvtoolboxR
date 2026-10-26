# Exact output-preservation receipt for separately timed installed packages.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
before <- readRDS(args[1])
after <- readRDS(args[2])
stopifnot(identical(names(before$outputs), names(after$outputs)))
rows <- data.frame(
  workload = names(before$outputs),
  identical = vapply(names(before$outputs), function(name)
    identical(before$outputs[[name]], after$outputs[[name]]), logical(1))
)
receipt <- list(
  baseline = args[1], candidate = args[2],
  inputs_identical = identical(before$inputs, after$inputs),
  outputs = rows
)
jsonlite::write_json(receipt, args[3], pretty = TRUE, auto_unbox = TRUE)
print(rows)
stopifnot(receipt$inputs_identical, all(rows$identical))
