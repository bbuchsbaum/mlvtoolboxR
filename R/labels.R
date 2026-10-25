#' Convert segmentation labels to an ink boundary mask
#' @param labels Finite integer-valued label matrix.
#' @return Logical matrix marking both sides of horizontal/vertical label
#'   transitions. Outer image borders are not added automatically.
#' @export
labels_to_edges <- function(labels) {
 if(!is.matrix(labels)||!is.numeric(labels)||any(!is.finite(labels)) ||
    any(labels!=floor(labels))||any(dim(labels)<1L))stop("labels must be a finite integer-valued matrix",call.=FALSE)
 h<-nrow(labels);w<-ncol(labels);out<-matrix(FALSE,h,w)
 if(h>1L) {
  diff<-labels[-1,,drop=FALSE]!=labels[-h,,drop=FALSE]
  out[-1,]<-out[-1,,drop=FALSE]|diff
  out[-h,]<-out[-h,,drop=FALSE]|diff
 }
 if(w>1L) {
  diff<-labels[,-1,drop=FALSE]!=labels[,-w,drop=FALSE]
  out[,-1]<-out[,-1,drop=FALSE]|diff
  out[,-w]<-out[,-w,drop=FALSE]|diff
 }
 out
}
