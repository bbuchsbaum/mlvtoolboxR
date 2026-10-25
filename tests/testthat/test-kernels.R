brute_distance <- function(mask) {
 sites<-which(mask);rc<-arrayInd(seq_along(mask),dim(mask));pt<-arrayInd(sites,dim(mask))
 d<-matrix(Inf,nrow(mask),ncol(mask));idx<-matrix(NA_integer_,nrow(mask),ncol(mask))
 if(length(sites))for(i in seq_along(mask)) {
  ds<-rowSums(sweep(pt,2,rc[i,],"-")^2);at<-which.min(ds);d[i]<-sqrt(ds[at]);idx[i]<-sites[at]
 }
 list(distance=d,index=idx)
}
test_that("EDT matches independent exhaustive search including ties", {
 set.seed(12)
 for(k in seq_len(20)) {
  m<-matrix(runif(8*11)<.1,8,11)
  expect_equal(distance_transform(m),brute_distance(m),tolerance=1e-12)
 }
 for(m in list(matrix(FALSE,3,4),matrix(TRUE,3,4)))
  expect_equal(distance_transform(m),brute_distance(m))
 expect_error(distance_transform(matrix(NA,2,2)),"binary")
})
scalar_aof <- function(D,index,samples=60,epsilon=1) {
 h<-nrow(D);w<-ncol(D);out<-matrix(0,h,w)
 circle<-cbind(cos((2*pi/samples)*(0:(samples-1))),sin((2*pi/samples)*(0:(samples-1))))
 for(r in seq_len(h))for(c in seq_len(w)) {
  if(!is.finite(D[r,c])||D[r,c]<=-1.5||r<=epsilon||c<=epsilon||r>=h-epsilon||c>=w-epsilon)next
  for(s in seq_len(samples)) {
   point<-c(r,c)+circle[s,]*epsilon+.5;cell<-floor(point)
   if(any(cell<1)||cell[1]>h||cell[2]>w)next
   label<-index[cell[1],cell[2]];if(is.na(label))next
   v<-arrayInd(label,dim(D))[1,]-point;norm<-sqrt(sum(v^2))
   if(norm>0)out[r,c]<-out[r,c]+sum(v/norm*circle[s,])
  }
 };out
}
test_that("AOF matches scalar geometry and leaves inputs intact", {
 m<-matrix(FALSE,9,13);m[5,7]<-TRUE;dt<-distance_transform(m);old<-dt
 expect_equal(average_outward_flux(dt$distance,dt$index),scalar_aof(dt$distance,dt$index),tolerance=1e-11)
 expect_identical(dt,old)
 expect_error(average_outward_flux(dt$distance,matrix(0L,9,13)),"indices")
 expect_error(average_outward_flux(dt$distance,dt$index,samples=0),"samples")
 p<-rbind(c(0,0),c(2,0),c(0,2),c(2,2));q<-rbind(c(1,1),c(.1,.2),c(4,2))
 n<-mlvtoolbox:::cpp_nearest(p,q)
 for(i in seq_len(nrow(q))) {
  ds<-rowSums(sweep(p,2,q[i,],"-")^2)
  expect_equal(n$index[i],which.min(ds));expect_equal(n$distance[i],sqrt(min(ds)))
 }
})
test_that("raster clipping and loop tracing are deterministic", {
 d<-mlv_drawing(list(matrix(c(-100,3,100,3),1)),7,5)
 m<-rasterize(d);expect_equal(sum(m),7);expect_true(all(m[3,]))
 loop<-matrix(FALSE,9,11);loop[3,3:8]<-TRUE;loop[7,3:8]<-TRUE;loop[3:7,3]<-TRUE;loop[3:7,8]<-TRUE
 traced<-trace_contours(loop,thin=FALSE)
 expect_length(traced$contours,1);expect_identical(rasterize(traced),loop)
 expect_error(trace_rgb(array(0,c(3,4,3)),function(z)matrix(FALSE,2,2)),"dimensions")
 m<-matrix(FALSE,15,19);m[4:12,5:15]<-TRUE;s<-skeletonize(m)
 expect_true(any(s));expect_identical(skeletonize(s),s)
 isolated<-matrix(FALSE,9,9);isolated[3,3]<-TRUE
 expect_equal(sum(skeletonize(isolated,2)),0)
})
test_that("branch scores have analytical and degenerate behavior", {
 b<-data.frame(x=1:21,y=rep(4,21),radius=rep(3,21))
 expect_equal(score_branch(b,"parallelism")[4:18],rep(1,15))
 expect_equal(score_branch(b,"separation"),rep((1-1/3)^10,21))
 expect_equal(score_branch(b,"taper")[4:18],rep(1,15))
 expect_equal(score_branch(b,"mirror"),rep((1-1/21)^10,21))
 expect_length(score_branch(b[FALSE,],"parallelism"),0)
 expect_equal(score_branch(b[1,],"parallelism"),0)
 expect_equal(score_branch(transform(b,radius=1),"separation"),rep(0,21))
 expect_true(all(is.na(score_branch(transform(b,radius=.5),"separation"))))
 expect_error(score_branch(transform(b,radius=-1)),"nonnegative")
})
test_that("MAT retains flux and valid zero mappings", {
 mask<-matrix(FALSE,35,45);mask[5:30,10]<-TRUE;mask[5:30,30]<-TRUE
 mat<-compute_mat(mask,min_branch_points=2);original<-mat$AOF
 expect_true(any(mat$skeleton));expect_true(all(is.finite(mat$distance_map)))
 rated<-compute_mat_properties(mat,c("parallelism","mirror","separation","taper"))
 expect_identical(rated$AOF,original)
 empty<-compute_mat(matrix(FALSE,5,7))
 expect_false(any(empty$skeleton));expect_true(all(is.infinite(empty$distance_map)))
 d<-mlv_drawing(list(matrix(c(10,3,10,12),1)),20,20)
 branch<-data.frame(row=3:12,column=rep(8L,10),x=rep(8,10),y=3:12,radius=rep(2,10),AOF=rep(1,10),separation=rep(0,10))
 custom<-structure(list(mask=rasterize(d),skeleton=matrix(FALSE,20,20),branches=list(branch),
   score_images=list(separation=matrix(0,20,20)),parameters=list()),class="mlv_mat")
 mapped<-map_mat_properties(d,custom)
 expect_true(any(mapped$mat_mapping_masks$separation))
 expect_equal(average_property(mapped,"separation"),0)
 expect_equal(sum(contour_statistics(mapped,"separation")$value),10)
})

test_that("AOF indexing agrees for spatially varying nearest-boundary fields", {
 set.seed(72)
 for(k in 1:3) {
  mask<-matrix(runif(11*13)<.15,11,13);dt<-distance_transform(mask)
  expect_equal(average_outward_flux(dt$distance,dt$index),
               scalar_aof(dt$distance,dt$index),tolerance=1e-11)
 }
})

test_that("degenerate masks and segmentation labels have explicit conventions", {
 full<-compute_mat(matrix(TRUE,5,7))
 expect_false(any(full$skeleton));expect_true(all(full$distance_map==-Inf))
 expect_true(all(full$AOF==0))
 labels<-matrix(c(1,1,2,2,2,2),2,3)
 edges<-labels_to_edges(labels)
 expect_identical(edges,matrix(c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE),2,3))
 expect_identical(labels_to_edges(matrix(1,1,1)),matrix(FALSE,1,1))
 expect_error(labels_to_edges(matrix(.1,1,1)),"integer")
 expect_equal(nrow(batch_statistics(list())),0)
})

test_that("deferred branch tracing preserves scoring results", {
 mask<-matrix(FALSE,25,35);mask[4:22,8]<-TRUE;mask[4:22,27]<-TRUE
 eager<-compute_mat(mask,min_branch_points=2)
 lazy<-compute_mat(mask,min_branch_points=2,trace_branches=FALSE)
 expect_length(lazy$branches,0)
 expect_equal(compute_mat_properties(lazy)$score_images,
              compute_mat_properties(eager)$score_images)
})
