.feature_schema_from_blocks <- function(blocks, normalize, spec=NULL) {
 columns <- lapply(blocks,function(b) {
  prop <- b$property
  data.frame(name=paste0(prop,"_",b$bin),property=prop,bin=b$bin,
   label=b$label,lower=b$lower,upper=b$upper,center=b$center,
   measurement_units=if(!is.null(b$measurement_units))b$measurement_units else if(prop=="orientation"||prop=="junction_angle")"degrees" else
    if(prop=="length")"pixels" else if(prop=="curvature")"degrees/pixel" else "dimensionless",
   weight=if(prop=="length")"contour_length" else
    if(prop %in% c("orientation","curvature","horizontal_vertical"))"segment_length" else
    if(grepl("^junction_",prop))"event_count" else "finite_mapped_pixel_count",
   boundary=if(prop=="orientation")"axial_180_upper_edge_wrap" else
    if(prop=="junction_type")"category" else "upper_edge_clamp_outside",
   units=b$units,stringsAsFactors=FALSE)
 })
 out<-structure(list(schema_version=1L,estimator="mlv2023-port-0.1.1",
  normalization=if(normalize)"per_10000_contour_pixels" else "none",
  columns=if(length(columns))do.call(rbind,columns) else data.frame()),
  class="mlv_feature_schema")
 if(!is.null(spec)) {
  out$schema_version<-2L;out$estimator<-spec$settings$estimator
  if(normalize)out$normalization<-paste0("per_10000_contour_",spec$settings$units)
  out$measurement_spec<-.spec_semantics(spec)
 }
 out
}
.check_feature_schema <- function(schema) {
 version2<-is.list(schema)&&identical(schema$schema_version,2L)
 if(!inherits(schema,"mlv_feature_schema") ||
    !identical(names(schema),c("schema_version","estimator","normalization","columns",if(version2)"measurement_spec")) ||
    (!identical(schema$schema_version,1L)&&!version2) ||
    !is.character(schema$estimator)||length(schema$estimator)!=1L||is.na(schema$estimator)||
    !is.character(schema$normalization)||length(schema$normalization)!=1L||
    is.na(schema$normalization)||!schema$normalization %in% c("none","per_10000_contour_pixels","per_10000_contour_visual_degrees") ||
    !is.data.frame(schema$columns))stop("invalid or unsupported feature schema",call.=FALSE)
 if(nrow(schema$columns) &&
    (!identical(names(schema$columns),c("name","property","bin","label","lower","upper","center",
      "measurement_units","weight","boundary","units")) ||
     anyNA(schema$columns)||anyDuplicated(schema$columns$name)||
     any(!nzchar(schema$columns$name))))stop("invalid feature schema columns",call.=FALSE)
 if(version2&&(!is.list(schema$measurement_spec)||!identical(schema$measurement_spec$version,1L)))
  stop("invalid measurement specification in schema",call.=FALSE)
 if(version2) {
  semantic<-schema$measurement_spec
  args<-semantic$settings
  if(identical(args$units,"visual_degrees"))args$pixels_per_degree<-1
  rebuilt<-tryCatch(.spec_semantics(do.call(feature_spec,args)),error=function(e)NULL)
  if(!identical(semantic,rebuilt))stop("invalid or unsupported measurement specification in schema",call.=FALSE)
 }
 schema
}
#' Inspect the measurement schema of a feature table
#'
#' A schema records ordered predictor names, actual bin boundaries, weights,
#' units, normalization and the measurement definition version. Canvas size is
#' represented only through the resulting bin boundaries, not as a separate
#' compatibility condition. Row subsets retain the schema; renamed, reordered
#' or removed predictor columns are rejected.
#' @param x A feature table returned by `feature_table()`.
#' @return An `mlv_feature_schema`, suitable for binding or checked prediction.
#' @export
feature_schema <- function(x) {
 if(!is.data.frame(x))stop("x must be a feature table",call.=FALSE)
 schema <- attr(x,"feature_schema",exact=TRUE)
 if(is.null(schema))stop("missing training feature schema; use feature_table()",call.=FALSE)
 schema <- .check_feature_schema(schema)
 if(!identical(names(x),c("contour_id",schema$columns$name)))
  stop("feature schema mismatch: predictor names/order",call.=FALSE)
 schema
}
.compare_feature_schemas <- function(expected,actual) {
 expected <- .check_feature_schema(expected);actual <- .check_feature_schema(actual)
 fields <- names(expected)
 bad <- fields[!vapply(fields,function(k)identical(expected[[k]],actual[[k]]),logical(1))]
 if("columns" %in% bad) {
  columns <- union(names(expected$columns),names(actual$columns))
  bad <- c(setdiff(bad,"columns"),paste0("columns.",columns[!vapply(columns,
   function(k)identical(expected$columns[[k]],actual$columns[[k]]),logical(1))]))
 }
 if(length(bad))stop("feature schema mismatch: ",paste(bad,collapse=", "),call.=FALSE)
 invisible(TRUE)
}
#' Bind a prediction model to its training feature schema
#'
#' Bind using the exact feature table used for training. The wrapper does not
#' train the model or infer its training history. It retains the schema through
#' R serialization without relying on a model to preserve data-frame attributes.
#' @param model An R prediction model or function accepting a predictor table.
#' @param features The training feature table, or its `feature_schema()`.
#' @return An `mlv_feature_model` containing the model and training schema.
#' @export
bind_feature_model <- function(model,features) {
 schema <- if(is.data.frame(features))feature_schema(features) else .check_feature_schema(features)
 if(inherits(model,"mlv_feature_model"))stop("model is already bound",call.=FALSE)
 structure(list(model=model,schema=schema),class="mlv_feature_model")
}
