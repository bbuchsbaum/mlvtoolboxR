#' Plot several property panels
#' @param x A drawing.
#' @param properties Properties passed to plot.
#' @param columns Panel column count.
#' @param ... Additional drawing plot arguments.
#' @return The drawing invisibly; graphics parameters are restored.
#' @export
plot_properties <- function(x,properties=c("orientation","length","curvature","junctions"),
                             columns=2L,...) {
 .check_drawing(x)
 if(!length(properties)||length(columns)!=1L||!is.finite(columns)||columns<1||columns!=floor(columns))
  stop("invalid panel layout",call.=FALSE)
 old<-graphics::par(mfrow=c(ceiling(length(properties)/columns),columns))
 on.exit(graphics::par(old),add=TRUE)
 for(prop in properties)graphics::plot(x,property=prop,main=prop,...)
 invisible(x)
}
#' Summarize a collection of drawings
#' @param drawings List of drawings.
#' @param image_ids Unique drawing identifiers; defaults to names or sequence.
#' @param ... Arguments to contour_statistics.
#' @return Long statistics data frame with image_id.
#' @export
batch_statistics <- function(drawings,image_ids=names(drawings),...) {
 if(!is.list(drawings)||inherits(drawings,"mlv_drawing"))stop("drawings must be a list of drawings",call.=FALSE)
 if(is.null(image_ids))image_ids<-as.character(seq_along(drawings))
 if(length(image_ids)!=length(drawings)||anyNA(image_ids)||anyDuplicated(image_ids))
  stop("invalid image_ids",call.=FALSE)
 out<-lapply(seq_along(drawings),function(i) {
  z<-contour_statistics(drawings[[i]],...)
  result<-data.frame(image_id=rep(image_ids[i],nrow(z)),z)
  attr(result,"measurement_metadata")<-attr(z,"measurement_metadata",exact=TRUE)
  result
 })
 if(!length(out))return(data.frame(image_id=character(),contour_statistics(mlv_drawing(list(),1,1),properties=character())))
 metadata<-lapply(out,attr,which="measurement_metadata",exact=TRUE)
 for(i in seq_along(metadata)).compare_feature_schemas(metadata[[1]]$schema,metadata[[i]]$schema)
 result<-do.call(rbind,out)
 attr(result,"measurement_metadata")<-list(schema=metadata[[1]]$schema,inputs=stats::setNames(metadata,image_ids))
 result
}
