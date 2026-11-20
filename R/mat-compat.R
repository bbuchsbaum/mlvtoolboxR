# Published-compatible ("mlv2023") branch rule: MLV traceSkeleton via Kovesi
# edgelink + cleanedgelist, MLV fitLineSegments, and mapMATtoContour order.
# These reproduce published MLV values, including recorded defects
# (MLV-M016; MLV-M028 direction dependence; MLV-M029 last-writer mapping and
# enumeration-order ties; MLV-M030 duplicated merge points).

# Line-by-line port of Kovesi's cleanedgelist (as vendored with MLV), kept only
# as a small-input test reference for cpp_cleanedgelist: it holds dense (2E)^2
# node and E^2 edge matrices, so it is unsuitable for production-size edge lists.
.mlv_cleanedgelist <- function(edgelist,minlength) {
 st<-new.env()
 st$E<-edgelist;Nedges<-length(edgelist);Nnodes<-2L*Nedges
 node<-matrix(0,Nnodes,2)
 for(n in seq_len(Nedges)){node[2*n-1,]<-edgelist[[n]][1,];node[2*n,]<-edgelist[[n]][nrow(edgelist[[n]]),]}
 key<-node[,1]*1e6+node[,2]
 st$A<-outer(key,key,"==")*1
 diag(st$A)<-0
 eof<-function(k)(k+1L)%/%2L
 st$B<-matrix(0,Nedges,Nedges)
 w<-which(st$A==1,arr.ind=TRUE)
 if(nrow(w))st$B[cbind(eof(w[,1]),eof(w[,2]))]<-1
 st$Nc<-colSums(st$A);st$Ec<-colSums(st$B)
 empty<-function(n)is.null(st$E[[n]])
 len<-function(e)sum(sqrt(rowSums((e[-nrow(e),,drop=FALSE]-e[-1,,drop=FALSE])^2)))
 removeedge<-function(n){
  st$E[n]<-list(NULL)
  st$Ec<-st$Ec-st$B[n,];st$Ec[n]<-0
  st$B[n,]<-0;st$B[,n]<-0
  d<-c(2*n-1,2*n)
  st$Nc<-st$Nc-st$A[d[1],];st$Nc<-st$Nc-st$A[d[2],]
  st$A[d,]<-0;st$A[,d]<-0
 }
 mergenodes<-function(n1,n2){
  e1<-eof(n1);e2<-eof(n2);s1<-2*e1-1;f1<-2*e1;s2<-2*e2-1;f2<-2*e2
  if(e1==e2)return(invisible())
  if(!st$A[n1,n2])stop("Attempt to merge nodes that are not connected",call.=FALSE)
  flip1<-n1%%2==1;flip2<-n2%%2==0
  a<-st$E[[e1]];b<-st$E[[e2]]
  if(!flip1&&!flip2){st$E[[e1]]<-rbind(a,b);st$A[f1,]<-st$A[f2,];st$A[,f1]<-st$A[,f2];st$Nc[f1]<-st$Nc[f2]}
  else if(!flip1&&flip2){st$E[[e1]]<-rbind(a,b[nrow(b):1,,drop=FALSE]);st$A[f1,]<-st$A[s2,];st$A[,f1]<-st$A[,s2];st$Nc[f1]<-st$Nc[s2]}
  else if(flip1&&!flip2){st$E[[e1]]<-rbind(a[nrow(a):1,,drop=FALSE],b)
   st$A[s1,]<-st$A[f1,];st$A[,s1]<-st$A[,f1];st$A[f1,]<-st$A[f2,];st$A[,f1]<-st$A[,f2]
   st$Nc[s1]<-st$Nc[f1];st$Nc[f1]<-st$Nc[f2]}
  else {st$E[[e1]]<-rbind(a[nrow(a):1,,drop=FALSE],b[nrow(b):1,,drop=FALSE])
   st$A[s1,]<-st$A[f1,];st$A[,s1]<-st$A[,f1];st$A[f1,]<-st$A[s2,];st$A[,f1]<-st$A[,s2]
   st$Nc[s1]<-st$Nc[f1];st$Nc[f1]<-st$Nc[s2]}
  st$B[e1,]<-pmax(st$B[e1,],st$B[e2,]);st$B[,e1]<-pmax(st$B[,e1],st$B[,e2]);st$B[e1,e1]<-0
  st$Ec<-colSums(st$B);st$Nc<-colSums(st$A)
  removeedge(e2)
 }
 info<-function(n){
  if(empty(n))return(list(spurdegree=0,spurnode=0,startnode=0,sconns=0,endnode=0,econns=0))
  sn<-2*n-1;en<-2*n;sc<-st$Nc[sn];ec<-st$Nc[en]
  if(sc==0&&ec>=1){sd<-ec;sp<-en} else if(sc>=1&&ec==0){sd<-sc;sp<-sn} else {sd<-0;sp<-0}
  list(spurdegree=sd,spurnode=sp,startnode=sn,sconns=sc,endnode=en,econns=ec)
 }
 for(n in seq_len(Nedges)) {
  if(!st$B[n,n]&&!empty(n)) {
   z<-info(n)
   if(z$sconns==1)mergenodes(which(st$A[z$startnode,]==1),z$startnode)
   if(!empty(n)&&z$econns==1)mergenodes(which(st$A[z$endnode,]==1),z$endnode)
  }
 }
 if(minlength>0) {
  for(n in seq_len(Nedges)) {
   z<-info(n)
   if(!empty(n)&&len(st$E[[n]])<minlength) {
    if(!st$Ec[n]||(st$Ec[n]==1&&st$B[n,n]==1)) removeedge(n)
    else if(z$spurdegree==2) {
     linking<-which(st$B[n,]==1)
     if(length(linking)==1L) removeedge(n)
     else {
      spurs<-n;lens<-len(st$E[[n]])
      for(i in seq_along(linking)) if(info(linking[i])$spurdegree) {
       spurs<-c(spurs,linking[i]);lens<-c(lens,len(st$E[[linking[i]]]))}
      del<-spurs[which.min(lens)]
      sp<-info(del)$spurnode
      m2<-which(st$A[sp,]==1)
      if(length(m2)!=2L)stop("attempt to merge other than 2 nodes",call.=FALSE)
      removeedge(del);mergenodes(m2[1],m2[2])
     }
    } else if(z$spurdegree==3) removeedge(n)
   }
  }
  for(n in seq_len(Nedges)) if(!empty(n)&&len(st$E[[n]])<minlength)
   if(!st$Ec[n]||(st$Ec[n]==1&&st$B[n,n]==1))removeedge(n)
 }
 Filter(Negate(is.null),st$E)
}
.trace_edgelink <- function(skeleton,minlength=10) {
 thin<-skeletonize(skeleton,method="hmt8")
 e<-cpp_edgelink_trace(thin)
 if(length(e)&&minlength>0)e<-cpp_cleanedgelist(e,minlength)
 lapply(e,function(z){storage.mode(z)<-"integer";z})
}

# Tangent samples in exactly the order of mapMATtoContour/getTangentPointsContour:
# [all first-side points; all second-side points], each ordered by segment i,
# then r1 in R[i]+(-.1,0,.1), then r2 in R[i+1]+(-.1,0,.1); coincident points
# contribute the skeleton point itself. Source = first point i of the segment.
.tangent_points_mlv <- function(branch) {
 X<-branch$row;Y<-branch$column;R<-branch$radius;n<-length(X)
 if(n<2L)return(list(points=matrix(numeric(),0,2),source=integer()))
 i<-rep(seq_len(n-1L),each=9L);j1<-rep(rep(1:3,each=3L),n-1L);j2<-rep(rep(1:3,3L),n-1L)
 dlt<-c(-.1,0,.1)
 x1<-X[i];y1<-Y[i];x2<-X[i+1L];y2<-Y[i+1L]
 r1<-c(R[i]+dlt[j1]);r2<-c(R[i+1L]+dlt[j2])
 r1<-ifelse(j1==2L,R[i],r1);r2<-ifelse(j2==2L,R[i+1L],r2)
 d<-sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
 zero<-d==0
 ca<-ifelse(zero,0,(r2-r1)/d)   # coincident points are replaced by the skeleton point below
 big<-abs(ca)>1;ca[big]<-ca[big]/(abs(ca[big])+0.00001)
 sa<-sqrt(1-ca*ca)
 mx<-.5*x1+.5*x2;my<-.5*y1+.5*y2;mr<-.5*r1+.5*r2;vx<-x1-mx;vy<-y1-my
 f1x<-ca*vx+sa*vy;f1y<--sa*vx+ca*vy;f2x<-ca*vx-sa*vy;f2y<-sa*vx+ca*vy
 s1<-mr/sqrt(f1x*f1x+f1y*f1y);s2<-mr/sqrt(f2x*f2x+f2y*f2y)
 p1<-cbind(mx+f1x*s1,my+f1y*s1);p2<-cbind(mx+f2x*s2,my+f2y*s2)
 p1[zero,]<-cbind(x1,y1)[zero,];p2[zero,]<-cbind(x1,y1)[zero,]
 list(points=rbind(p1,p2),source=c(i,i))
}
# MLV mirror scores: fitLineSegments on [X,Y] = (row, column) with the published
# terminal-span rule, then movmean(3) and the exponent, as computeMATpropertyPerBranch.
.mirror_mlv <- function(branch,power) {
 n<-nrow(branch);out<-numeric(n)
 if(n>3L)out<-cpp_fit_lines_mlv(cbind(as.double(branch$row),as.double(branch$column)))
 .moving_mean(out,3L)^power
}
