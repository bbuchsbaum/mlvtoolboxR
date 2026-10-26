test_that("public designer reconstructs verified conditions and preserves identity", {
 z<-designer_v1_fixtures()$set_knapsack
 result<-do.call(design_subsets,c(list(drawings=z$drawings),z$arguments))
 expect_equal(result$status,"optimal");expect_equal(result$diagnostics$objective,270,tolerance=1e-8)
 expect_equal(result$solver$absolute_gap,0)
 expect_identical(as.data.frame(result),result$selection)
 expect_output(print(result),"optimal")
 for(i in seq_along(z$drawings)) {
  original<-z$drawings[[i]]
  ids<-unlist(lapply(c(result$conditions,list(unassigned=result$unassigned)),function(s)s[[i]]$ids),use.names=FALSE)
  expect_setequal(ids,original$ids);expect_equal(anyDuplicated(ids),0L)
  for(set in result$conditions) {
   x<-set[[i]]
   expect_null(x$mat_properties);expect_null(x$junctions)
   expect_identical(x$contours,original$contours[match(x$ids,original$ids)])
   expect_identical(x$feature_spec,z$arguments$spec)
  }
 }
 a<-match(result$selection$condition,c("unassigned","A","B"))-1L
 expect_true(reference_evaluate(z$drawings,z$arguments,a)$feasible)
 expect_equal(length(result$diagnostics$per_drawing),2)
 expect_identical(result$input_geometry,lapply(z$drawings,mlvtoolbox:::.geometry_snapshot))
 repeated<-do.call(design_subsets,c(list(drawings=z$drawings),z$arguments))
 expect_identical(result$selection,repeated$selection)
 file<-tempfile();on.exit(unlink(file));saveRDS(result,file)
 expect_identical(readRDS(file),result)
})

test_that("no-solution and time-limited public outputs are unambiguous", {
 cases<-designer_v1_fixtures()
 for(name in c("per_drawing_infeasible","exhausted")) {
  z<-cases[[name]];out<-do.call(design_subsets,c(list(drawings=z$drawings),z$arguments))
  expect_identical(out$status,z$expected$status)
  expect_null(out$conditions);expect_null(out$diagnostics)
  expect_true(all(is.na(out$selection$condition)))
  expect_output(print(out),out$status)
 }
 z<-cases$ties;z$arguments$time_limit<-0
 out<-do.call(design_subsets,c(list(drawings=z$drawings),z$arguments))
 expect_equal(out$status,"budget_exhausted")
 expect_equal(out$solver$nodes,0)
})

test_that("requested relations are fresh diagnostics and never selection constraints", {
 z<-designer_v1_fixtures()$ties
 out<-do.call(design_subsets,c(list(drawings=z$drawings,diagnostics="junctions"),z$arguments))
 expect_equal(nrow(out$relational_diagnostics$A$tie$junctions),0)
 expect_equal(out$solver$objective,0)
 z$arguments$spec$settings$mat$min_branch_points<-2L
 out<-do.call(design_subsets,c(list(drawings=z$drawings,diagnostics="separation"),z$arguments))
 expect_s3_class(out$relational_diagnostics$A$tie$mat$mat,"mlv_mat")
 expect_null(out$conditions$A$tie$mat_properties)
 expect_error(do.call(design_subsets,c(list(drawings=z$drawings,diagnostics="coverage"),z$arguments)),"diagnostics")
})

test_that("physical inputs retain distinct calibrations through subset reconstruction", {
 z<-designer_v1_fixtures()$small_set;a<-z$arguments
 a$spec<-feature_spec(c("length","curvature","orientation"),normalize=FALSE,
   units="visual_degrees",pixels_per_degree=1,estimator="arclength",curvature_scale=.1)
 z$drawings[[2]]$contours<-lapply(z$drawings[[2]]$contours,function(x)x*2)
 out<-do.call(design_subsets,c(list(drawings=z$drawings,pixels_per_degree=c(1,2)),a))
 expect_equal(out$status,"optimal")
 expect_equal(out$diagnostics$objective,90,tolerance=1e-8)
 expect_equal(out$input_specifications$b$settings$pixels_per_degree,2)
 expect_equal(out$conditions$B$b$feature_spec$settings$pixels_per_degree,2)
})
