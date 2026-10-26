# Installed-package benchmark. Setup and validation are outside warmed timings.
args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
make_drawing<-function(n) {
 contours<-lapply(seq_len(n),function(i) {
  x<-seq(2,20,length.out=21)+(i-1)%%32*30
  y<-10+3*sin(seq(0,2*pi,length.out=21))+(i-1)%/%32*30
  cbind(head(x,-1),head(y,-1),tail(x,-1),tail(y,-1))
 })
 mlv_drawing(contours,1000,ceiling(n/32)*30+20)
}
many<-make_drawing(1000L)
grid<-mlv_drawing(c(lapply(seq(10,190,10),function(x)matrix(c(x,5,x,195),1)),
 lapply(seq(10,190,10),function(y)matrix(c(5,y,195,y),1))),200,200)
theta<-seq(0,2*pi,length.out=121)
rings<-mlv_drawing(lapply(c(25,45,70,95),function(radius) {
 x<-128+radius*cos(theta);y<-128+.8*radius*sin(theta)
 cbind(head(x,-1),head(y,-1),tail(x,-1),tail(y,-1))
}),256,256)
properties<-c("orientation","length","curvature","horizontal_vertical")
jobs<-list(features_1000x20=function()feature_table(many,properties,bins=16),
 statistics_1000x20=function()contour_statistics(many,properties,bins=16,by="contour"),
 junction_grid_19x19=function()detect_junctions(grid),
 mat_rings_256=function()analyse_mat(rings,min_branch_points=3),
 raster_rings_256=function()rasterize(rings,line_width=3))
if(length(args)>=3)jobs<-jobs[grepl(args[3],names(jobs))]
rows<-list();outputs<-list()
for(name in names(jobs)) {
 f<-jobs[[name]]
 first<-unname(system.time(value<-f())[3])
 number<-1L
 repeat {elapsed<-unname(system.time(for(i in seq_len(number))f())[3]);if(elapsed>=.15||number>=128)break;number<-number*2L}
 times<-replicate(5,unname(system.time(for(i in seq_len(number))f())[3])/number)
 rows[[name]]<-data.frame(workload=name,first_seconds=first,median_seconds=median(times),
  min_seconds=min(times),max_seconds=max(times),iterations=number,output_bytes=as.numeric(object.size(value)))
 outputs[[name]]<-value
 Rprof(file.path(args[2],paste0(name,".Rprof")),interval=.001)
 for(k in 1:3)f()
 Rprof(NULL)
 capture.output(summaryRprof(file.path(args[2],paste0(name,".Rprof")))$by.total,
                file=file.path(args[2],paste0(name,"-profile.txt")))
 cat(name,median(times),"seconds\n")
}
out<-do.call(rbind,rows);out$source_commit<-system2("git",c("rev-parse","HEAD"),stdout=TRUE)
out$R<-R.version.string
write.csv(out,file.path(args[2],"timings.csv"),row.names=FALSE)
saveRDS(list(inputs=list(many=many,grid=grid,rings=rings),outputs=outputs),
 file.path(args[2],"outputs.rds"))
