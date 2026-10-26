# Save independently installed, qualified 0.2.0 public outputs before optimization.
args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
stopifnot(as.character(packageVersion("mlvtoolbox"))=="0.2.0")
set.seed(64102)
mk<-function(n) {
 p<-cbind(cumsum(c(20,rnorm(n,0,2))),cumsum(c(20,rnorm(n,0,2))))
 cbind(p[-nrow(p),,drop=FALSE],p[-1,,drop=FALSE])
}
angles<-c(0,90,180,270,360,180/13/2,90-1e-12,90+1e-12,270-1e-12,270+1e-12,
          seq(0,360,length.out=33))
axis<-mlv_drawing(lapply(angles,function(a)matrix(c(0,0,5*cos(a*pi/180),5*sin(a*pi/180)),1)),100,100)
random<-mlv_drawing(lapply(c(1,2,8,21),mk),100,100)
w<-readRDS(system.file("fixtures","statistics-before-optimization.rds",package="mlvtoolbox"))
drawings<-list(random=random,axis=axis,mixed=w$drawing,empty=mlv_drawing(list(),100,100))
cases<-list()
record<-function(fun,args)cases[[length(cases)+1L]]<<-list(fun=fun,args=args,expected=do.call(get(fun),args))
for(name in names(drawings)) {
 d<-drawings[[name]]
 for(bins in c(1,5,13,31))for(norm in c(FALSE,TRUE))for(by in c("drawing","contour"))
  record("contour_statistics",list(x=d,properties=c("orientation","length","curvature","horizontal_vertical"),bins=bins,by=by,normalize=norm))
 for(props in list("length",c("orientation","orientation"),character()))
  record("feature_table",list(x=d,properties=props,bins=13))
}
for(estimator in c("mlv2023","turning","arclength"))for(unit in c("pixels","visual_degrees")) {
 s<-feature_spec(c("orientation","length","curvature","horizontal_vertical"),
   bins=7,units=unit,pixels_per_degree=20,estimator=estimator,curvature_scale=2)
 for(d in list(random,drawings$empty)) {
  record("feature_table",list(x=d,spec=s))
  record("contour_statistics",list(x=d,spec=s,by="drawing"))
 }
 record("feature_table",list(x=random,spec=s,pixels_per_degree=40))
}
for(case in w$cases)record("contour_statistics",c(list(x=w$drawing),case$args))
saveRDS(list(package="0.2.0",source_manifest="docs/evidence/release-a/package-source-manifest.json",
 seed=64102,cases=cases),args[2],compress="xz")
cat("Saved",length(cases),"baseline cases\n")
