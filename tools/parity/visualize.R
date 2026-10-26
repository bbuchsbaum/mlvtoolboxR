args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()));library(mlvtoolbox)
base<-"inst/fixtures/parity";dest<-args[2];dir.create(dest,recursive=TRUE,showWarnings=FALSE)
read<-function(name)unname(as.matrix(read.csv(file.path(base,name),header=FALSE)))
names<-c("rectangle","ring","cross","diagonal","sparse","isolated");results<-list();rows<-list()
for(name in names) {
 mask<-read(paste0(name,".input.csv"))>0;m<-compute_mat(mask,trace_branches=FALSE)
 a<-read(paste0(name,".octave.aof.csv"));s<-read(paste0(name,".octave.skeleton.csv"))>0
 common<-m$skeleton&s;union<-m$skeleton|s
 rows[[name]]<-data.frame(case=name,R_skeleton_pixels=sum(m$skeleton),
  octave_skeleton_pixels=sum(s),different_skeleton_pixels=sum(m$skeleton!=s),
  intersection_over_union=if(any(union))sum(common)/sum(union) else 1,
  nearest_label_disagreements=sum(m$nearest_boundary!=read(paste0(name,".octave.index.csv"))),
  max_AOF_difference=max(abs(m$AOF-a)),
  max_distance_difference=max(abs(m$distance_map-read(paste0(name,".octave.distance.csv")))))
 results[[name]]<-list(mask=mask,mat=m,aof=a,skeleton=s)
}
write.csv(do.call(rbind,rows),file.path(dest,"mat-topologies.csv"),row.names=FALSE)
colors<-grDevices::colorRampPalette(c("#2166ac","#f7f7f7","#b2182b"))(129)
layer<-function(z,palette,limits,title)graphics::image(t(z[nrow(z):1,,drop=FALSE]),
 col=palette,zlim=limits,axes=FALSE,xlab="",ylab="",main=title,useRaster=TRUE,asp=1)
png(file.path(dest,"mat-parity.png"),width=1650,height=1150,res=140)
par(mfrow=c(4,5),mar=c(.5,.5,2,.5),oma=c(3,0,3,0))
for(name in c("rectangle","ring","cross","sparse")) {
 r<-results[[name]];range<-max(abs(c(r$mat$AOF,r$aof)));diff<-r$mat$AOF-r$aof
 layer(r$mask,c("white","black"),c(0,1),paste(name,"input"))
 layer(r$mat$AOF,colors,c(-range,range),"R AOF")
 layer(r$aof,colors,c(-range,range),"Octave AOF")
 layer(diff,colors,c(-range,range),"R minus Octave")
 overlay<-ifelse(r$mat$skeleton&r$skeleton,1,ifelse(r$mat$skeleton,2,ifelse(r$skeleton,3,0)))
 layer(overlay,c("white","black","#d73027","#2166ac"),c(0,3),"Skeleton overlap")
}
mtext("Equal distances can still give different MATs at nearest-pixel ties",outer=TRUE,line=1,font=2,cex=1.15)
mtext("Skeletons: black = shared; red = R only; blue = Octave only. AOF panels share a scale within each row.",
 outer=TRUE,side=1,line=1,cex=.8)
dev.off()
