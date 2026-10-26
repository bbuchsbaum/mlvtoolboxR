.bins_count <- function(bins) {
 if(length(bins)!=1L || !is.finite(bins) || bins<1 || bins!=floor(bins) ||
    bins>.Machine$integer.max) stop("bins must be a positive integer",call.=FALSE)
 as.integer(bins)
}
.empty_statistics <- function() data.frame(contour_id=character(),property=character(),
 bin=integer(),label=character(),lower=double(),upper=double(),center=double(),
 value=double(),units=character())
.grouped_bins <- function(values,weights,groups,edges,n_groups,index=NULL) {
 k<-length(edges)-1L
 if(as.double(n_groups)*k>.Machine$integer.max)stop("histogram exceeds index capacity",call.=FALSE)
 out<-matrix(0,n_groups,k)
 if(!length(values)||!n_groups)return(out)
 if(is.null(index))index<-findInterval(values,edges,all.inside=TRUE)
 key<-groups+n_groups*(index-1L)
 # Preserve R's sum() and the within-bin input order, including its precision.
 totals<-vapply(split(weights,key),sum,double(1),na.rm=TRUE)
 out[as.integer(names(totals))]<-totals
 out
}
.statistics_blocks <- function(x,properties,bins,by,normalize,length_range,curvature_range) {
 .check_drawing(x);bins<-.bins_count(bins)
 allowed<-c("orientation","length","curvature","horizontal_vertical",
            "junction_type","junction_angle","mirror","parallelism","separation","taper")
 if(any(!properties %in% allowed))stop("unknown statistic property",call.=FALSE)
 for(r in list(length_range,curvature_range))
  if(length(r)!=2L||any(!is.finite(r))||r[1]<0||r[2]<=r[1])
   stop("log-bin ranges must be finite and increasing",call.=FALSE)
 ng<-if(by=="drawing")1L else length(x$ids)
 if(!ng||!length(properties))return(list())
 st<-.segment_table(x)
 lens<-vapply(x$contours,function(z)sum(sqrt(rowSums((z[,3:4,drop=FALSE]-z[,1:2,drop=FALSE])^2))),double(1))
 totals<-if(by=="drawing")sum(lens) else lens
 segment_group<-if(by=="drawing")rep.int(1L,nrow(st)) else st$contour_index
 contour_group<-if(by=="drawing")rep.int(1L,length(lens)) else seq_along(lens)
 if(any(grepl("^junction_",properties))) {
  j<-if(is.null(x$junctions))detect_junctions(x) else x$junctions
  membership<-lapply(j$members,function(m)unique(m$contour_index))
  ji<-if(by=="drawing")seq_len(nrow(j)) else rep.int(seq_len(nrow(j)),lengths(membership))
  jg<-if(by=="drawing")rep.int(1L,nrow(j)) else unlist(membership,use.names=FALSE)
 }
 blocks<-lapply(properties,function(prop) {
  labels<-NULL;unit<-"pixels";index<-NULL
  if(prop=="orientation") {
   width<-180/bins;edges<-seq(0,180,length.out=bins+1L)
   values<-st$orientation;weights<-st$length;groups<-segment_group
   # Scale the axial angle directly. Shifting/subtracting half a bin can put
   # 90 and 270 degrees on different sides of the same edge for odd bin counts.
   index<-(floor((values%%180)*bins/180+.5)%%bins)+1L
   centers<-seq(0,180-width,length.out=bins)
  } else if(prop=="horizontal_vertical") {
   edges<-seq(-1,1,length.out=bins+1L)
   values<-pmax(-1,pmin(1,abs(cos(st$orientation*pi/180))-abs(sin(st$orientation*pi/180))))
   weights<-st$length;groups<-segment_group;centers<-(head(edges,-1)+tail(edges,-1))/2
  } else if(prop %in% c("length","curvature")) {
   range<-if(prop=="length")length_range else curvature_range
   edges<-seq(log10(range[1]+1),log10(range[2]+1),length.out=bins+1L)
   values<-log10(if(prop=="length")lens+1 else st$curvature+1)
   weights<-if(prop=="length")lens else st$length
   groups<-if(prop=="length")contour_group else segment_group
   centers<-10^((head(edges,-1)+tail(edges,-1))/2)-1
  } else if(grepl("^junction_",prop)) {
   unit<-"count";groups<-jg;weights<-rep.int(1,length(ji))
   if(prop=="junction_type") {
    labels<-c("T","Y","X","Arrow","Star");edges<-seq(.5,5.5)
    values<-match(j$type[ji],labels);centers<-seq_len(5L)
   } else {
    edges<-seq(0,120,length.out=bins+1L);values<-j$min_angle[ji]
    centers<-(head(edges,-1)+tail(edges,-1))/2
   }
  } else {
   if(is.null(x$mat_properties[[prop]]))stop("MAT property not mapped: ",prop,call.=FALSE)
   pieces<-x$mat_properties[[prop]][seq_along(x$ids)]
   values<-unlist(pieces,use.names=FALSE)
   groups<-if(by=="drawing")rep.int(1L,length(values)) else rep.int(seq_along(pieces),lengths(pieces))
   finite<-is.finite(values);values<-values[finite];groups<-groups[finite]
   weights<-rep.int(1,length(values));edges<-seq(0,1,length.out=bins+1L)
   centers<-(head(edges,-1)+tail(edges,-1))/2;unit<-"count"
  }
  counts<-.grouped_bins(values,weights,groups,edges,ng,index)
  if(normalize) {
   counts<-counts/totals*10000
   counts[totals==0,]<-0
  }
  shown_edges<-if(prop %in% c("length","curvature"))10^edges-1 else edges
  if(prop=="orientation")shown_edges<-edges-90/bins
  list(counts=counts,property=prop,bin=seq_len(ncol(counts)),
       label=if(is.null(labels))as.character(signif(centers,8)) else labels,
       lower=head(shown_edges,-1),upper=tail(shown_edges,-1),center=centers,
       units=if(normalize)paste0(unit," per 10000 contour pixels") else unit)
 })
 blocks
}
#' Histograms of drawing and contour properties
#'
#' Orientation, curvature and horizontal/vertical histograms are weighted by
#' segment length. Length histograms collect total contour length. Junctions
#' count events; per-contour junctions count each participating contour once.
#' MAT histograms count mapped pixels, including valid zero scores. Values
#' outside the logarithmic range enter the first/last bin. Orientation is axial;
#' exact bin boundaries enter the upper bin, including odd bin counts.
#' @param x A drawing.
#' @param properties Property names.
#' @param bins Positive number of bins (junction types always have five).
#' @param by Summarize the drawing or each contour.
#' @param normalize Scale by contour length to units per 10,000 pixels.
#' @param length_range,curvature_range Increasing nonnegative log-bin ranges.
#' @return A long data frame with bin boundaries, centers, values and units.
#' @export
contour_statistics <- function(x,properties=c("orientation","length","curvature",
 "horizontal_vertical","junction_type","junction_angle"),bins=8,
 by=c("drawing","contour"),normalize=FALSE,length_range=c(2,max(3,x$width+x$height)),
 curvature_range=c(0,90)) {
 by<-match.arg(by)
 blocks<-.statistics_blocks(x,properties,bins,by,normalize,length_range,curvature_range)
 if(!length(blocks))return(.empty_statistics())
 counts<-do.call(cbind,lapply(blocks,`[[`,"counts"))
 widths<-vapply(blocks,function(b)ncol(b$counts),integer(1))
 flatten<-function(field)unlist(lapply(blocks,`[[`,field),use.names=FALSE)
 ng<-nrow(counts);per_group<-ncol(counts)
 data.frame(contour_id=rep(if(by=="drawing")NA_character_ else x$ids,each=per_group),
  property=rep(rep(properties,widths),ng),bin=rep(flatten("bin"),ng),
  label=rep(flatten("label"),ng),lower=rep(flatten("lower"),ng),
  upper=rep(flatten("upper"),ng),center=rep(flatten("center"),ng),
  value=as.vector(t(counts)),units=rep(rep(flatten("units"),widths),ng),row.names=NULL)
}
#' Wide per-contour feature table
#' @param x A drawing.
#' @param properties Properties to include.
#' @param bins Number of bins.
#' @param normalize Normalize by contour length.
#' @return One row per contour, with named numeric predictors.
#' @export
feature_table <- function(x,properties=c("orientation","length","curvature","junction_type"),
                          bins=8,normalize=TRUE) {
 blocks<-.statistics_blocks(x,unique(properties),bins,"contour",normalize,
                            c(2,max(3,x$width+x$height)),c(0,90))
 out<-data.frame(contour_id=x$ids)
 if(!length(blocks))return(out)
 counts<-do.call(cbind,lapply(blocks,function(b) {
  colnames(b$counts)<-paste0(b$property,"_",b$bin);b$counts
 }))
 data.frame(out,counts,check.names=FALSE,row.names=NULL)
}
#' Predict one score per contour
#' @param x A drawing.
#' @param model An R model with a predict method, or a function of a feature table.
#' @param properties,bins Feature-table settings.
#' @param ... Additional model prediction arguments.
#' @return A finite numeric vector in contour order.
#' @export
predict_contours <- function(x,model,properties=c("orientation","length","curvature","junction_type"),
                             bins=8,...) {
 f <- feature_table(x,properties,bins)
 predictors <- f[,setdiff(names(f),"contour_id"),drop=FALSE]
 result <- if(is.function(model))model(predictors,...) else stats::predict(model,newdata=predictors,...)
 if(!is.numeric(result) || length(result)!=length(x$ids) || any(!is.finite(result)))
  stop("model must return one finite numeric score per contour",call.=FALSE)
 as.numeric(result)
}
