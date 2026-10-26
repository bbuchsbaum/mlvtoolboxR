# Fixed before prototype execution. All coordinates use pixels; canvas is 40x40.
voronoi_fixtures <- function() {
 draw<-function(z)mlv_drawing(lapply(seq_len(nrow(z)),function(i)z[i,,drop=FALSE]),40,40)
 parallel<-rbind(c(5,10,25,10),c(5,14,25,14))
 list(
  parallel=list(drawing=draw(parallel),axis_y=12,axis_length=20),
  shifted_subquantum=list(drawing=draw(parallel+.0004),axis_y=12.0004,axis_length=20),
  shifted_quantum=list(drawing=draw(parallel+.0006),axis_y=12.0006,axis_length=20),
  shifted_tenth=list(drawing=draw(parallel+.1),axis_y=12.1,axis_length=20),
  narrow=list(drawing=draw(rbind(c(5,10,25,10),c(5,10.4,25,10.4))),axis_y=10.2,axis_length=20),
  rectangle=list(drawing=draw(rbind(c(5,5,25,5),c(25,5,25,15),c(25,15,5,15),c(5,15,5,5)))),
  crossing=list(drawing=draw(rbind(c(5,5,25,25),c(5,25,25,5)))),
  curved=list(drawing=draw(rbind(c(5,10,25,10),c(10,15,20,15),c(30,5,30,25)))),
  collapse=list(drawing=draw(rbind(c(5,10,25,10),c(5,10.0004,25,10.0004))),error="overlap|topology"),
  near_miss=list(drawing=draw(rbind(c(5,10,25,10),c(15,10.0004,15,20))),error="topology"),
  overlap=list(drawing=draw(rbind(c(5,10,25,10),c(15,10,30,10))),error="overlap"),
  overflow=list(drawing=draw(rbind(c(2200000,0,2200001,0))),error="int32"))
}
