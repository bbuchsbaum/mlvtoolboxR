schema_drawing <- function(size=200)mlv_drawing(list(matrix(c(1,1,101,1),1)),size,size)

test_that("canvas-dependent predictor meanings are rejected before calling a model", {
 small<-schema_drawing();large<-schema_drawing(2000)
 f<-feature_table(small,"length");g<-feature_table(large,"length")
 expect_equal(which(as.numeric(f[1,-1])>0),6L)
 expect_equal(which(as.numeric(g[1,-1])>0),4L)
 calls<-0L;model<-function(z){calls<<-calls+1L;rep(7,nrow(z))}
 bound<-bind_feature_model(model,f)
 expect_error(predict_contours(large,bound,properties="length"),"columns.lower")
 expect_error(predict_contours(small,model,properties="length"),"missing training")
 expect_equal(calls,0L)
 expect_identical(predict_contours(small,bound,properties="length"),7)
 expect_equal(calls,1L)
 expect_error(predict_contours(small,bound,properties="length",bins=4),"columns")
 expect_error(predict_contours(small,bound,properties="length",normalize=FALSE),"normalization")
 expect_error(predict_contours(small,bound,properties="length",length_range=c(2,500)),"columns.lower")
 expect_equal(calls,1L)
 # A common explicit range makes the meaning independent of unused canvas area.
 expect_identical(predict_contours(large,bound,properties="length",length_range=c(2,400)),7)
 expect_warning(expect_identical(predict_contours(large,model,properties="length",check_spec=FALSE),7),"unchecked prediction")
})

test_that("schema binds ordered predictor semantics and survives serialization", {
 d<-schema_drawing();f<-feature_table(d,c("length","curvature"))
 s<-feature_schema(f)
 expect_identical(s$columns$name,names(f)[-1])
 expect_identical(s$columns$weight,rep(c("contour_length","segment_length"),each=8))
 expect_identical(s$normalization,"per_10000_contour_pixels")
 expect_identical(feature_schema(f[1,,drop=FALSE]),s)
 bad<-f;names(bad)[2]<-"renamed"
 expect_error(feature_schema(bad),"names/order")
 expect_error(feature_schema(f[,c(1,3,2,4:ncol(f))]),"missing training|names/order")
 expect_error(feature_schema(f[,-2]),"missing training|names/order")
 reordered<-f[,c(1,3,2,4:ncol(f))];attr(reordered,"feature_schema")<-s
 expect_error(feature_schema(reordered),"names/order")
 expect_error(bind_feature_model(function(z)1,bad),"names/order")
 raw<-f;attr(raw,"feature_schema")<-NULL
 expect_error(feature_schema(raw),"missing training")
 model<-bind_feature_model(function(z)sum(z$length_6),f)
 restored<-unserialize(serialize(model,NULL))
 expect_identical(restored$schema,s)
 expect_identical(predict_contours(d,restored,properties=c("length","curvature")),10000)
 expect_error(predict_contours(d,restored,properties=c("curvature","length")),"columns.name")
 expect_error(bind_feature_model(model,f),"already bound")
 expect_error(predict_contours(d,model,expected_spec=feature_schema(feature_table(d,"length"))),"columns")
})

test_that("actual boundaries, axial conventions and empty tables have schemas", {
 d<-schema_drawing()
 f<-feature_table(d,"orientation",bins=5)
 s<-feature_schema(f)
 expect_equal(s$columns$lower,seq(-18,126,36))
 expect_equal(s$columns$upper,seq(18,162,36))
 expect_identical(s$columns$boundary,rep("axial_180_upper_edge_wrap",5))
 h<-contour_statistics(d,"length",bins=5,length_range=c(2,500))
 f<-feature_table(d,"length",bins=5,length_range=c(2,500))
 expect_identical(feature_schema(f)$columns$lower,h$lower)
 empty<-feature_table(mlv_drawing(list(),200,200),"length",bins=5,length_range=c(2,500))
 expect_identical(feature_schema(empty),feature_schema(f))
 expect_equal(nrow(empty),0)
 no_properties<-feature_table(d,character())
 expect_s3_class(feature_schema(no_properties),"mlv_feature_schema")
 expect_equal(nrow(feature_schema(no_properties)$columns),0)
})

test_that("schema validation rejects unsupported versions and missing fields", {
 d<-schema_drawing();f<-feature_table(d,"length");s<-feature_schema(f)
 for(field in names(s)) {
  bad<-s;bad[[field]]<-NULL
  expect_error(bind_feature_model(function(z)1,bad),"invalid|unsupported")
 }
 bad<-s;bad$schema_version<-2L
 expect_error(bind_feature_model(function(z)1,bad),"unsupported")
 bad<-s;bad$columns$name[2]<-bad$columns$name[1]
 expect_error(bind_feature_model(function(z)1,bad),"invalid")
 expect_error(feature_schema(1),"feature table")
 expect_error(feature_table(d,"length",normalize=NA),"normalize")
 expect_error(predict_contours(d,function(z)1,check_spec=NA),"check_spec")
})
