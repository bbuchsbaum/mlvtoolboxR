parity_read<-function(name)unname(as.matrix(utils::read.csv(
 system.file("fixtures","parity",name,package="mlvtoolbox"),header=FALSE)))
parity_drawing<-function(name) {
 z<-parity_read(paste0(name,".input.csv"))
 mlv_drawing(lapply(unique(z[,1]),function(i)z[z[,1]==i,2:5,drop=FALSE]),256,192)
}
test_that("diverse contour measurements agree with executed source references", {
 for(name in c("integer0","integer1","integer2","subpixel","fractional")) {
  d<-parity_drawing(name);z<-as.data.frame(d)
  actual<-unname(as.matrix(z[,c("length","orientation","curvature")]))
  expect_equal(actual,parity_read(paste0(name,".octave.geometry.csv")),tolerance=1e-12)
  expect_equal(unname(do.call(rbind,rotate_drawing(d,37)$contours)),
               parity_read(paste0(name,".octave.rotate37.csv")),tolerance=1e-12)
  if(grepl("^integer",name))
   expect_lt(max(abs(actual-parity_read(paste0(name,".python.geometry.csv")))),3e-8)
  for(bins in c(5,8,13))for(prop in c("length","orientation","curvature"))for(norm in c(FALSE,TRUE)) {
   h<-contour_statistics(d,prop,bins,by="contour",normalize=norm)
   actual<-matrix(h$value,ncol=bins,byrow=TRUE)
   if(!(bins==13&&prop=="orientation")) {
    field<-paste0(if(norm)"norm_" else "",prop,bins,".csv")
    expect_equal(actual,parity_read(paste0(name,".octave.",field)),tolerance=1e-12)
    if(grepl("^integer",name))
     expect_equal(actual,parity_read(paste0(name,".python.",if(prop=="length")"adapted_" else "",field)),tolerance=1e-12)
   }
   target<-if(norm)rep(10000,length(d$ids)) else vapply(d$contours,function(s)sum(sqrt(rowSums((s[,3:4,drop=FALSE]-s[,1:2,drop=FALSE])^2))),double(1))
   expect_equal(rowSums(actual),target,tolerance=1e-12)
  }
 }
})
test_that("odd orientation bins obey axial reversal and a defined upper-edge rule", {
 d<-mlv_drawing(list(matrix(c(5,5,5,15),1),matrix(c(5,15,5,5),1)),20,20)
 h<-feature_table(d,"orientation",bins=13,normalize=FALSE)
 expect_identical(as.numeric(h[1,-1]),as.numeric(h[2,-1]))
 expect_equal(which(h[1,-1]>0),8L)
 expect_equal(h$orientation_8,c(10,10))
 for(bins in c(1,3,5,7,13,31)) {
  h<-feature_table(d,"orientation",bins=bins)
  expect_identical(as.numeric(h[1,-1]),as.numeric(h[2,-1]))
  expect_equal(sum(h[1,-1]),10000)
 }
})
test_that("grouped statistics preserve the checked non-boundary output exactly", {
 w<-readRDS(system.file("fixtures","statistics-before-optimization.rds",package="mlvtoolbox"))
 for(case in w$cases) {
  actual<-do.call(contour_statistics,c(list(w$drawing),case$args))
  expect_s3_class(attr(actual,"measurement_metadata")$schema,"mlv_feature_schema")
  attr(actual,"measurement_metadata")<-NULL
  expect_witness(actual,case$expected)
 }
 features<-feature_table(w$drawing,c("length","orientation","curvature","parallelism"),13)
 expect_s3_class(feature_schema(features),"mlv_feature_schema")
 attr(features,"feature_schema")<-NULL
 attr(features,"measurement_metadata")<-NULL
 expect_witness(features,w$features)
 expect_identical(feature_table(w$drawing,c("length","length","curvature")),
                  feature_table(w$drawing,c("length","curvature")))
 expect_equal(nrow(contour_statistics(mlv_drawing(list(),10,10),"length")),8)
 expect_equal(nrow(contour_statistics(mlv_drawing(list(),10,10),"length",by="contour")),0)
 expect_error(contour_statistics(w$drawing,bins=2^31),"positive integer")
})
test_that("diverse MAT stages agree when the nearest-boundary field is shared", {
 for(name in c("rectangle","ring","cross","diagonal","sparse","isolated")) {
  ink<-parity_read(paste0(name,".input.csv"))>0
  m<-compute_mat(ink,trace_branches=FALSE)
  D<-parity_read(paste0(name,".octave.distance.csv"))
  I<-parity_read(paste0(name,".octave.index.csv"))
  A<-parity_read(paste0(name,".octave.aof.csv"))
  S<-parity_read(paste0(name,".octave.skeleton.csv"))>0
  expect_lt(max(abs(m$distance_map-D)),2e-6)
  expect_lt(max(abs(average_outward_flux(D,I)-A)),1e-10)
  expect_identical(skeletonize(A>(2/pi*sin(28*pi/180))*60),S)
  h<-nrow(D);r<-row(D);c<-col(D);a<-m$nearest_boundary-1;b<-I-1
  expect_equal(unname((a%%h+1-r)^2+(a%/%h+1-c)^2),
                unname((b%%h+1-r)^2+(b%/%h+1-c)^2),tolerance=0)
 }
 ink<-parity_read("isolated.input.csv")>0
 expect_equal(compute_mat(ink)$distance_map[ink],c(-1,-1))
})
test_that("branch comparisons preserve correct contracts and expose the terminal fit defect", {
 for(name in c("constant","varying","curved","terminal_bend","short")) {
  z<-parity_read(paste0(name,".branch.csv"));b<-data.frame(x=z[,1],y=z[,2],radius=z[,3])
  ref<-parity_read(paste0(name,".octave.scores.csv"))
  expect_lt(max(abs(score_branch(b,"parallelism")-ref[,1])),1e-10)
  expect_lt(max(abs(score_branch(b,"separation")-ref[,2])),1e-10)
  if(name %in% c("constant","varying","short"))
   expect_lt(max(abs(score_branch(b,"mirror")-ref[,3])),1e-10)
 }
 xy<-parity_read("terminal_bend.branch.csv")[,1:2]
 fit<-fit_line_segments(xy)
 for(i in seq_len(nrow(fit$spans))) {
  z<-xy[fit$spans[i,1]:fit$spans[i,2],,drop=FALSE]
  if(nrow(z)>2) {
   chord<-z[nrow(z),]-z[1,];unit<-chord/sqrt(sum(chord^2))
   offset<-sweep(z[-c(1,nrow(z)),,drop=FALSE],2,z[1,],"-")
   perpendicular<-offset-outer(as.numeric(offset%*%unit),unit)
   expect_lte(mean(sqrt(rowSums(perpendicular^2))),1)
  }
 }
 expect_gt(nrow(fit$spans),1)
})
test_that("spatially filtered MAT mapping agrees with exhaustive nearest samples", {
 d<-mlv_drawing(list(matrix(c(2,4,35,4),1),matrix(c(2,16,35,16),1)),40,20)
 branches<-list(
  data.frame(row=rep(7,8),column=4:11,radius=rep(3,8),parallelism=c(0,.2,NA,.4,.6,.8,1,1)),
  data.frame(row=rep(12,8),column=25:32,radius=rep(4,8),parallelism=rep(.5,8)))
 mat<-structure(list(mask=rasterize(d),branches=branches,
  score_images=list(parallelism=matrix(NA_real_,20,40)),parameters=list()),class="mlv_mat")
 ink<-which(mat$mask);rc<-arrayInd(ink,dim(mat$mask))
 for(radius in c(.1,1,4,100)) {
  expected<-matrix(NA_real_,20,40)
  for(b in branches) {
   fp<-mlvtoolbox:::.tangent_points(b)
   for(q in seq_along(ink)) {
    ds<-rowSums(sweep(fp$points,2,rc[q,],"-")^2);at<-which.min(ds)
    value<-b$parallelism[fp$source[at]]
    if(sqrt(ds[at])<radius&&is.finite(value))
     expected[ink[q]]<-max(c(expected[ink[q]],value),na.rm=TRUE)
   }
  }
  actual<-map_mat_properties(d,mat,radius)
  expect_equal(actual$mat_images$parallelism,expected,tolerance=0)
  expect_identical(actual$mat_mapping_masks$parallelism,is.finite(expected))
 }
})
test_that("EDT and AOF agree with reference fields at multiple sampling radii", {
 for(name in c("features_random","features_ties")) {
  mask<-parity_read(paste0(name,".input.csv"))>0
  dt<-distance_transform(mask)
  D<-parity_read(paste0(name,".octave.distance.csv"))
  I<-parity_read(paste0(name,".octave.index.csv"))
  expect_lt(max(abs(dt$distance-D)),2e-6)
  expect_lt(max(abs(dt$distance-parity_read(paste0(name,".python.distance.csv")))),1e-12)
  for(eps in c(1,2)) {
   ref<-parity_read(paste0(name,".octave.aof",eps,".csv"))
   expect_lt(max(abs(average_outward_flux(D,I,epsilon=eps)-ref)),1e-10)
   expect_gt(max(abs(parity_read(paste0(name,".python.aof",eps,".csv"))-ref)),1)
  }
 }
})

test_that("tangent hoisting preserves every saved floating-point sample exactly", {
 cases<-readRDS(system.file("fixtures","tangents-before-optimization.rds",package="mlvtoolbox"))
 for(case in cases)expect_identical(mlvtoolbox:::.tangent_points(case$branch),case$expected)
})
