library(mlvtoolbox)
source("tools/voronoi/fixtures.R")
cases<-voronoi_fixtures();destination<-"docs/evidence/voronoi-spike"
png(file.path(destination,"parallel-gap-comparison.png"),width=1500,height=1050,res=150)
par(mfrow=c(2,2),mar=c(4.5,4.8,3,1),mgp=c(2.8,.8,0))
for(name in c("parallel","narrow"))for(method in c("Boost finite linear probe","Raster AOF")) {
 x<-cases[[name]]$drawing;z<-as.data.frame(x)
 ylim<-if(name=="parallel")c(9,15) else c(9.8,10.6)
 plot(NA,xlim=c(4,26),ylim=ylim,xlab="x (pixels)",ylab="y (pixels)",
      main=paste(name,method,sep=" — "))
 segments(z$x1,z$y1,z$x2,z$y2,col="gray45",lwd=2)
 if(method=="Boost finite linear probe") {
  edge<-read.csv(file.path(destination,paste0(name,"-selected.csv")))
  segments(edge$x1,edge$y1,edge$x2,edge$y2,col="#0072B2",lwd=2)
 } else {
  p<-read.csv(file.path(destination,paste0(name,"-raster-skeleton.csv")))
  points(p$column,p$row,pch=16,cex=.6,col="#D55E00")
 }
}
dev.off()
