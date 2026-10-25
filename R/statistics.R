.bins_count <- function(bins) {
 if(length(bins)!=1L || !is.finite(bins) || bins<1 || bins!=floor(bins)) stop("bins must be a positive integer",call.=FALSE)
 as.integer(bins)
}
.weighted_bins <- function(values,weights,edges) {
 k <- length(edges)-1L
 if(!length(values)) return(numeric(k))
 idx <- findInterval(values,edges,all.inside=TRUE)
 out <- numeric(k)
 for(i in seq_len(k)) out[i] <- sum(weights[idx==i],na.rm=TRUE)
 out
}
#' Histograms of drawing and contour properties
#'
#' Orientation, curvature and horizontal/vertical histograms are weighted by
#' segment length. Length histograms collect total contour length. Junctions
#' count events; per-contour junctions count each participating contour once.
#' MAT histograms count mapped pixels, including valid zero scores. Values
#' outside the logarithmic range enter the first/last bin.
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
 .check_drawing(x); bins <- .bins_count(bins); by <- match.arg(by)
 allowed <- c("orientation","length","curvature","horizontal_vertical",
              "junction_type","junction_angle","mirror","parallelism","separation","taper")
 if(any(!properties %in% allowed)) stop("unknown statistic property",call.=FALSE)
 for(r in list(length_range,curvature_range))
  if(length(r)!=2L || any(!is.finite(r)) || r[1]<0 || r[2]<=r[1])
   stop("log-bin ranges must be finite and increasing",call.=FALSE)
 st <- .segment_table(x)
 lens <- vapply(x$contours,function(z)sum(sqrt(rowSums((z[,3:4,drop=FALSE]-z[,1:2,drop=FALSE])^2))),double(1))
 j <- if(any(grepl("^junction_",properties))) {if(is.null(x$junctions))detect_junctions(x) else x$junctions} else NULL
 groups <- if(by=="drawing") list(seq_along(x$ids)) else lapply(seq_along(x$ids),identity)
 result <- list()
 for(ci in groups) {
  ss <- st[st$contour_index %in% ci,,drop=FALSE]
  total <- sum(lens[ci]); id <- if(by=="drawing") NA_character_ else x$ids[ci]
  for(prop in properties) {
   labels <- NULL; unit <- "pixels"
   if(prop=="orientation") {
    width <- 180/bins; edges <- seq(0,180,length.out=bins+1L)
    values <- (ss$orientation+width/2)%%180; weights <- ss$length
    centers <- seq(0,180-width,length.out=bins)
   } else if(prop=="horizontal_vertical") {
    edges <- seq(-1,1,length.out=bins+1L)
    values <- pmax(-1,pmin(1,abs(cos(ss$orientation*pi/180))-abs(sin(ss$orientation*pi/180))))
    weights <- ss$length; centers <- (head(edges,-1)+tail(edges,-1))/2
   } else if(prop %in% c("length","curvature")) {
    range <- if(prop=="length") length_range else curvature_range
    edges <- seq(log10(range[1]+1),log10(range[2]+1),length.out=bins+1L)
    values <- log10(if(prop=="length")lens[ci]+1 else ss$curvature+1)
    weights <- if(prop=="length") lens[ci] else ss$length
    centers <- 10^((head(edges,-1)+tail(edges,-1))/2)-1
   } else if(grepl("^junction_",prop)) {
    jj <- j
    if(by=="contour" && nrow(jj)) jj <- jj[vapply(jj$members,function(m)ci %in% m$contour_index,logical(1)),,drop=FALSE]
    unit <- "count"
    if(prop=="junction_type") {
     labels <- c("T","Y","X","Arrow","Star"); edges <- seq(0.5,5.5)
     values <- match(jj$type,labels); weights <- rep(1,length(values)); centers <- seq_len(5L)
    } else {
     edges <- seq(0,120,length.out=bins+1L)
     values <- jj$min_angle; weights <- rep(1,length(values))
     centers <- (head(edges,-1)+tail(edges,-1))/2
    }
   } else {
    if(is.null(x$mat_properties[[prop]]))stop("MAT property not mapped: ",prop,call.=FALSE)
    values <- unlist(x$mat_properties[[prop]][ci],use.names=FALSE)
    values <- values[is.finite(values)]
    weights <- rep(1,length(values)); edges <- seq(0,1,length.out=bins+1L)
    centers <- (head(edges,-1)+tail(edges,-1))/2; unit <- "count"
   }
   counts <- .weighted_bins(values,weights,edges)
   if(normalize) counts <- if(total>0) counts/total*10000 else rep(0,length(counts))
   shown_edges <- if(prop %in% c("length","curvature")) 10^edges-1 else edges
   # Orientation boundaries are shifted relative to the displayed bin centers.
   if(prop=="orientation") shown_edges <- edges-90/bins
   z <- data.frame(contour_id=rep(id,length(counts)),property=prop,bin=seq_along(counts),
       label=if(is.null(labels))as.character(signif(centers,8)) else labels,
       lower=head(shown_edges,-1),upper=tail(shown_edges,-1),center=centers,
       value=counts,units=if(normalize)paste0(unit," per 10000 contour pixels") else unit)
   result[[length(result)+1L]] <- z
  }
 }
 if(!length(result)) return(data.frame(contour_id=character(),property=character(),bin=integer(),
   label=character(),lower=double(),upper=double(),center=double(),value=double(),units=character()))
 do.call(rbind,result)
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
 .check_drawing(x)
 z <- contour_statistics(x,properties,bins,by="contour",normalize=normalize)
 out <- data.frame(contour_id=x$ids)
 if(!nrow(z)) return(out)
 names <- paste0(z$property,"_",z$bin)
 for(n in unique(names)) {
  zz <- z[names==n,,drop=FALSE]
  out[[n]] <- zz$value[match(x$ids,zz$contour_id)]
 }
 out
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
