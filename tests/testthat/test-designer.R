test_that("bounded exact selection matches every frozen tiny optimum", {
 cases<-designer_v1_fixtures()
 for(name in names(cases)) {
  z<-cases[[name]];a<-z$arguments
  p<-do.call(mlvtoolbox:::.subset_prepare,c(list(drawings=z$drawings),a[setdiff(names(a),c("time_limit","max_nodes"))]))
  result<-mlvtoolbox:::.subset_solve(p,a$time_limit,a$max_nodes)
  if(identical(z$expected$status,"budget_exhausted")) {
   expect_equal(result$status,"budget_exhausted",info=name)
   expect_null(result$assignment);next
  }
  ref<-reference_exhaustive(z$drawings,a)
  if(ref$feasible==0) {
   expect_equal(result$status,"infeasible",info=name);expect_null(result$assignment);next
  }
  expect_equal(result$status,"optimal",info=name)
  expect_equal(result$objective,ref$optimum,tolerance=1e-8,info=name)
  check<-reference_evaluate(z$drawings,a,result$assignment)
  expect_true(check$feasible,info=name)
  expect_equal(result$objective,check$objective,tolerance=1e-8)
  expect_equal(result$absolute_gap,0)
  expect_identical(mlvtoolbox:::.subset_solve(p)$assignment,result$assignment)
 }
})

test_that("budgeted search reports honest bounds, with independent feasibility", {
 z<-designer_v1_fixtures()$knapsack;a<-z$arguments
 p<-do.call(mlvtoolbox:::.subset_prepare,c(list(drawings=z$drawings),a[setdiff(names(a),c("time_limit","max_nodes"))]))
 zero<-mlvtoolbox:::.subset_solve(p,time_limit=0)
 expect_equal(zero$status,"budget_exhausted");expect_equal(zero$nodes,0)
 # This cap is intentionally far below complete search; either interrupted
 # status is valid, but never an optimality/infeasibility certificate.
 limited<-mlvtoolbox:::.subset_solve(p,max_nodes=100)
 expect_true(limited$status%in%c("feasible","budget_exhausted"))
 expect_gte(limited$upper_bound,270-1e-8)
 if(!is.null(limited$assignment)) {
  expect_true(reference_evaluate(z$drawings,a,limited$assignment)$feasible)
  expect_gte(limited$upper_bound,limited$objective)
 }
 expect_error(mlvtoolbox:::.subset_solve(p,time_limit=-1),"time_limit")
 expect_error(mlvtoolbox:::.subset_solve(p,max_nodes=0),"max_nodes")
})

test_that("designer rejects unsupported or ambiguous measurement constraints", {
 z<-designer_v1_fixtures()$ties;a<-z$arguments
 prepare<-function(...) {
  args<-c(list(drawings=z$drawings),a[setdiff(names(a),c("time_limit","max_nodes"))])
  args[names(list(...))]<-list(...);do.call(mlvtoolbox:::.subset_prepare,args)
 }
 expect_error(prepare(drawings=list()),"nonempty")
 expect_error(prepare(spec=feature_spec()),"normalize=FALSE")
 expect_error(prepare(spec=feature_spec("junction_type",normalize=FALSE)),"additive")
 expect_error(prepare(nuisance=c(spatial_coverage=1)),"nuisance")
 expect_error(prepare(nuisance=c(length=-1)),"nuisance")
 expect_error(prepare(nuisance=1),"nuisance")
 expect_error(prepare(nuisance=c(horizontal_vertical=1)),"included")
 expect_error(prepare(count_bounds=c(1,1.5)),"count_bounds")
 expect_error(prepare(length_bounds=c(5,4)),"length_bounds")
 expect_error(prepare(count_imbalance=-1),"count_imbalance")
 expect_error(prepare(pixels_per_degree=c(10,20)),"calibration")
 expect_error(prepare(drawings=rep(z$drawings,17)),"unique")
 lots<-rep(unname(z$drawings),17)
 expect_error(prepare(drawings=lots),"64 eligible")
})

test_that("signed objectives and feasible per-drawing constraints are exact", {
 z<-designer_v1_fixtures()$ties;a<-z$arguments
 drawing<-z$drawings[[1]]
 drawing$contours[[2]]<-matrix(c(10,10,10,15),1)
 drawings<-list(x=drawing)
 a$target<-"orientation";a$nuisance<-numeric()
 p<-do.call(mlvtoolbox:::.subset_prepare,c(list(drawings=drawings),a[setdiff(names(a),c("time_limit","max_nodes"))]))
 out<-mlvtoolbox:::.subset_solve(p)
 expect_equal(out$objective,10,tolerance=1e-8)
 expect_equal(out$objective,reference_exhaustive(drawings,a)$optimum)
 expect_equal(out$assignment[2],2L)
 a$target<-"length";a$mode<-"per_drawing";drawings<-list(x=drawing,y=drawing)
 p<-do.call(mlvtoolbox:::.subset_prepare,c(list(drawings=drawings),a[setdiff(names(a),c("time_limit","max_nodes"))]))
 out<-mlvtoolbox:::.subset_solve(p)
 expect_equal(out$status,"optimal");expect_equal(out$objective,0,tolerance=1e-8)
 expect_true(reference_evaluate(drawings,a,out$assignment)$feasible)
 empty<-mlv_drawing(list(),10,10)
 p<-do.call(mlvtoolbox:::.subset_prepare,c(list(drawings=list(empty)),a[setdiff(names(a),c("time_limit","max_nodes"))]))
 expect_equal(mlvtoolbox:::.subset_solve(p)$status,"infeasible")
})
