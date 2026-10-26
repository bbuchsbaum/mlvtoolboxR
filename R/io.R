#' Convert an upstream MATLAB drawing representation
#' @param x A named list with contours and imsize fields.
#' @param zero_length Reject or explicitly drop zero-length segments.
#' @return An mlv_drawing.
#' @export
as_mlv_drawing <- function(x,zero_length=c("error","drop")) {
 zero_length<-match.arg(zero_length)
 if(inherits(x,"mlv_drawing"))return(x)
 unwrap<-function(v){while(is.list(v)&&length(v)==1L&&(is.null(names(v))||all(names(v)=="")))v<-v[[1]];v}
 if(!is.list(x)||is.null(x$contours)||is.null(x$imsize))stop("legacy drawing needs contours and imsize",call.=FALSE)
 dimensions<-as.numeric(unwrap(x$imsize))
 if(length(dimensions)!=2L)stop("legacy imsize must contain width,height",call.=FALSE)
 cells<-x$contours
 # Cell arrays can have one additional scalar-list wrapper around each matrix.
 if(is.matrix(cells)&&is.numeric(cells))cells<-list(cells)
 cells<-lapply(as.vector(cells),unwrap)
 if(!is.null(x$numContours)) {
  n<-as.numeric(unwrap(x$numContours))
  if(length(n)!=1L||n!=length(cells))stop("legacy numContours disagrees with contour cells",call.=FALSE)
 }
 mlv_drawing(cells,dimensions[1],dimensions[2],zero_length=zero_length,
  source=list(original_image=if(!is.null(x$originalImage))unwrap(x$originalImage) else NULL,
              method=if(!is.null(x$lineMethod))unwrap(x$lineMethod) else "MATLAB import"))
}
#' Read MATLAB v5 vector drawing structures
#' @param file MAT file path.
#' @param variable Name of the drawing structure array.
#' @param zero_length Zero-length input policy.
#' @return One drawing for a scalar structure, otherwise a list of drawings.
#' @export
read_mlv_mat <- function(file,variable="vecLD",zero_length=c("error","drop")) {
 if(!requireNamespace("R.matlab",quietly=TRUE))stop("install optional package 'R.matlab'",call.=FALSE)
 connection<-base::file(file,"rb");on.exit(close(connection),add=TRUE)
 header<-readBin(connection,"raw",n=128L)
 if(length(header)>=8L&&identical(header[1:8],as.raw(c(137,72,68,70,13,10,26,10))))
  stop("MAT v7.3/HDF5 is not supported; export a v5 MAT structure",call.=FALSE)
 if(grepl("MATLAB 7.3 MAT-file",paste(rawToChar(header,multiple=TRUE),collapse=""),fixed=TRUE))
  stop("MAT v7.3/HDF5 is not supported; export a v5 MAT structure",call.=FALSE)
 values<-R.matlab::readMat(file,fixNames=FALSE)
 object<-values[[variable]]
 if(is.null(object))stop("MAT variable not found: ",variable,call.=FALSE)
 dims<-dim(object);fields<-dimnames(object)[[1]]
 if(is.list(object)&&length(dims)>=2L&&!is.null(fields)) {
  nfield<-dims[1];n<-length(object)/nfield
  out<-lapply(seq_len(n),function(i) {
   z<-object[(i-1L)*nfield+seq_len(nfield)];names(z)<-fields
   as_mlv_drawing(z,zero_length)
  })
  if(length(out)==1L)out[[1]] else out
 } else as_mlv_drawing(object,zero_length)
}
#' Export image or contour statistics
#' @param x A drawing.
#' @param file CSV destination.
#' @param ... Arguments to contour_statistics.
#' @return The written data frame invisibly. A companion file named
#'   `<file>.schema.rds` preserves the measurement schema, calibration and provenance.
#' @export
write_feature_csv <- function(x,file,...) {
 result<-contour_statistics(x,...)
 utils::write.csv(result,file,row.names=FALSE,na="")
 saveRDS(attr(result,"measurement_metadata",exact=TRUE),paste0(file,".schema.rds"))
 invisible(result)
}
