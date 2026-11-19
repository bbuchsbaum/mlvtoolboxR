test_that("fused statistics preserve qualified 0.2.0 public values and metadata", {
 witness<-readRDS(system.file("fixtures","statistics-0.2.0.rds",package="mlvtoolbox"))
 for(i in seq_along(witness$cases)) {
  case<-witness$cases[[i]]
  expect_witness(do.call(getExportedValue("mlvtoolbox",case$fun),case$args),
                 case$expected,info=paste("baseline case",i,case$fun))
 }
})

test_that("native grouped sums preserve ordered R sums on adversarial inputs", {
 set.seed(9501)
 for(i in seq_len(50)) {
  n<-if(i==1)0 else sample(10:500,1);ng<-sample(1:9,1);nb<-sample(1:11,1)
  groups<-sample(c(seq_len(ng),NA_integer_),n,TRUE)
  bins<-sample(c(seq_len(nb),NA_integer_),n,TRUE)
  weights<-sample(c(-1e16,1e16,-1,0,1,1e-100,NA_real_,NaN),n,TRUE)*runif(n)
  expected<-matrix(0,ng,nb)
  for(g in seq_len(ng))for(b in seq_len(nb)) {
   keep<-which(!is.na(groups)&!is.na(bins)&groups==g&bins==b)
   expected[g,b]<-sum(weights[keep],na.rm=TRUE)
  }
  expect_identical(mlvtoolbox:::cpp_grouped_sums(weights,groups,bins,ng,nb),expected)
 }
 expect_identical(mlvtoolbox:::cpp_grouped_sums(numeric(),integer(),integer(),0L,3L),matrix(0,0,3))
 expect_identical(mlvtoolbox:::cpp_grouped_sums(rep(.Machine$double.xmax,2),c(1L,1L),c(1L,1L),1L,1L),matrix(Inf,1,1))
 expect_identical(mlvtoolbox:::cpp_grouped_sums(rep(-.Machine$double.xmax,2),c(1L,1L),c(1L,1L),1L,1L),matrix(-Inf,1,1))
 expect_identical(mlvtoolbox:::cpp_grouped_sums(c(Inf,-Inf),c(1L,1L),c(1L,1L),1L,1L),matrix(NaN,1,1))
})

test_that("native contour lengths retain separate squared terms and R summation", {
 set.seed(9502)
 for(scale in c(1e-140,1e-50,1,1e50,1e140)) {
  contours<-lapply(c(1,2,35),function(n)matrix(rnorm(n*4)*scale,n,4))
  expected<-vapply(contours,function(z)sum(sqrt(rowSums((z[,3:4,drop=FALSE]-z[,1:2,drop=FALSE])^2))),double(1))
  expect_identical(mlvtoolbox:::cpp_contour_lengths(contours),expected)
 }
 expect_identical(mlvtoolbox:::cpp_contour_lengths(list()),numeric())
 expect_error(mlvtoolbox:::cpp_contour_lengths(list(matrix(1,2,3))),"four contour columns")
})

test_that("compiled histogram entry points reject invalid dimensions and indices", {
 expect_error(mlvtoolbox:::cpp_grouped_sums(1,1L,1L,50000L,50000L),"index capacity")
 expect_error(mlvtoolbox:::cpp_grouped_sums(1,integer(),1L,1L,1L),"inconsistent")
 expect_error(mlvtoolbox:::cpp_grouped_sums(1,0L,1L,1L,1L),"index")
 expect_error(mlvtoolbox:::cpp_grouped_sums(1,1L,2L,1L,1L),"index")
 geo<-mlvtoolbox:::cpp_geometry_histograms
 expect_error(geo(list(),1L,list(),TRUE,1),"inconsistent")
 expect_error(geo(list(),0L,list(c(0,1)),TRUE,1),"property")
 expect_error(geo(list(),1L,list(0),TRUE,1),"edges")
 expect_error(geo(list(),1L,list(c(1,0)),TRUE,1),"edges")
 expect_error(geo(list(),1L,list(c(0,NA)),TRUE,1),"edges")
 expect_error(geo(list(),integer(),list(),TRUE,0),"scale")
 expect_error(geo(list(matrix(1,2,3)),integer(),list(),TRUE,1),"four contour columns")
})

test_that("fused summaries use current geometry and conserve physical length", {
 d<-mlv_drawing(list(matrix(c(0,0,3,4),1)),10,10)
 enriched<-compute_contour_properties(d,properties="length")
 enriched$contours[[1]][1,3:4]<-c(0,10)
 f<-feature_table(enriched,"length",normalize=FALSE)
 expect_equal(sum(as.numeric(f[1,-1])),10)
 s<-feature_spec(c("orientation","length","curvature","horizontal_vertical"),
   units="visual_degrees",pixels_per_degree=5,normalize=FALSE,bins=13)
 f<-feature_table(d,spec=s)
 for(p in s$settings$properties)expect_equal(sum(as.numeric(f[1,grep(paste0("^",p,"_"),names(f))])),1)
 empty<-contour_statistics(mlv_drawing(list(),1,1),"length",normalize=TRUE)
 expect_equal(empty$value,rep(0,8))
})


test_that("fused remainder handling preserves axes and subnormal directions", {
 tiny<-.Machine$double.xmin*.Machine$double.eps
 d<-mlv_drawing(lapply(c(tiny,-tiny,1e-300,-1e-300,0,1,-1),function(y)
  matrix(c(0,0,1,y),1)),10,10)
 st<-as.data.frame(d)
 for(bins in c(1L,5L,13L,31L)) {
  expected<-matrix(0,length(d$ids),bins)
  index<-(floor((st$orientation%%180)*bins/180+.5)%%bins)+1L
  expected[cbind(seq_along(d$ids),index)]<-st$length
  actual<-feature_table(d,"orientation",bins=bins,normalize=FALSE)
  expect_identical(unname(as.matrix(actual[,-1,drop=FALSE])),expected)
 }
 expect_identical(mlvtoolbox:::cpp_contour_lengths(list(matrix(c(0L,0L,3L,4L),1))),5)
 expect_error(mlvtoolbox:::cpp_contour_lengths(list(matrix("x",1,4))),"numeric contour")
 expect_error(mlvtoolbox:::cpp_geometry_histograms(list(matrix("x",1,4)),integer(),list(),TRUE,1),"numeric contour")
})


test_that("invalid edited geometry is rejected before native bin indexing", {
 d<-mlv_drawing(list(matrix(c(0,0,1,1),1)),10,10)
 for(endpoint in list(c(0,0),c(NA_real_,1),c(Inf,1))) {
  invalid<-d;invalid$contours[[1]][1,3:4]<-endpoint
  expect_error(feature_table(invalid,"orientation"),"finite positive segment lengths")
 }
})
