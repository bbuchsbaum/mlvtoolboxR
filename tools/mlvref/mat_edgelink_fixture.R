# Track B: synthetic skeletons for the edgelink/cleanedgelist port. Writes masks to
# inst/fixtures/edgelink/<name>.mask.csv; tools/mlvref/mat_edgelink_fixture.m then
# runs the pinned MLV traceSkeleton/GetConSeg (edgelink, minlength 10) under Octave
# and writes <name>.octave.branches.csv (branch, row, col). Synthetic data only.
source("tools/mlvref/mat_common.R")
out <- "inst/fixtures/edgelink"
d1 <- mlv_drawing(list(matrix(c(8,8,52,8, 52,8,52,52, 52,52,8,52, 8,52,8,8),4,byrow=TRUE),
                       matrix(c(30,8,30,52),1), matrix(c(8,30,52,30),1)), 60, 60)
d2 <- mlv_drawing(list(matrix(c(5,10,55,14),1), matrix(c(5,30,55,24),1), matrix(c(10,50,50,40, 50,40,58,55),2,byrow=TRUE),
                       matrix(c(20,2,25,58),1)), 60, 60)
d3 <- mlv_drawing(list(matrix(c(30,5,55,30, 55,30,30,55, 30,55,5,30, 5,30,30,5),4,byrow=TRUE), matrix(c(30,5,30,55),1),
                       matrix(c(12,12,20,20),1)), 60, 60)
for (nm in c("grid", "fan", "kite")) {
  d <- switch(nm, grid = d1, fan = d2, kite = d3)
  sk <- compute_mat(d, aof = "mlv2023", trace_branches = FALSE, min_component = 3)$skeleton
  write.table(sk * 1L, file.path(out, paste0(nm, ".mask.csv")), sep = ",", row.names = FALSE, col.names = FALSE)
  cat(nm, sum(sk), "skeleton pixels\n")
}
