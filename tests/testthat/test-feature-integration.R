test_that("curvature histograms integrate the affine profile, not cell means", {
 d<-contract_drawing(contract_xy("elbow"))
 s<-feature_spec("curvature",bins=2,normalize=FALSE,curvature_range=c(0,80),
                 estimator="arclength",curvature_scale=2)
 f<-feature_table(d,spec=s)
 # Boundary sqrt(81)-1 = 8; the tent exceeds 8 on 4*(1-8/45) units.
 high<-4*(1-8/45)
 expect_equal(unname(as.numeric(f[1,-1])),c(30-high,high),tolerance=1e-12)
 for(xy in list(contract_xy("elbow_split"),contract_xy("elbow_split")[5:1,])) {
  z<-feature_table(contract_drawing(xy),spec=s)
  expect_equal(z,f,tolerance=1e-10)
 }
 st<-contour_statistics(d,spec=s)
 expect_equal(st$value,c(30-high,high),tolerance=1e-12)
 expect_equal(mlvtoolbox:::.linear_histogram(c(-2,2,10),c(2,2,20),c(4,3,5),
  c(1,1,1),c(0,2,4),1),matrix(c(4,8),1))
 t<-feature_table(d,spec=feature_spec("curvature",bins=2,normalize=FALSE,
                    estimator="turning",curvature_range=c(0,80)))
 expect_equal(as.numeric(t[1,-1]),c(30,0))
})

test_that("physical features, schema and model bindings survive calibration changes", {
 d<-contract_drawing(contract_xy("elbow")*20)
 e<-contract_drawing(contract_xy("elbow")*40);e$width<-2000L;e$height<-2000L
 s<-feature_spec(c("length","curvature","orientation"),units="visual_degrees",
  pixels_per_degree=20,estimator="arclength",curvature_scale=2)
 f<-feature_table(d,spec=s);g<-feature_table(e,spec=s,pixels_per_degree=40)
 expect_equal(as.matrix(f[,-1]),as.matrix(g[,-1]),tolerance=1e-10)
 expect_identical(feature_schema(f),feature_schema(g))
 expect_equal(attr(g,"measurement_metadata")$calibration$pixels_per_degree,40)
 model<-bind_feature_model(function(z)rowSums(z),f)
 expect_equal(predict_contours(e,model,spec=s,pixels_per_degree=40),sum(as.numeric(f[1,-1])))
 bad<-model;bad$schema$measurement_spec$definitions['contour']<-"unknown"
 expect_error(predict_contours(d,bad,spec=s),"unsupported")
 mismatch<-s;mismatch$settings$curvature_scale<-3
 expect_error(predict_contours(d,model,spec=mismatch),"mismatch")
 expect_error(feature_table(d,spec=s,bins=3),"spec or")
 expect_error(contour_statistics(d,spec=s,normalize=FALSE),"spec or")
 expect_error(predict_contours(d,model,spec=s,bins=3),"spec or")
 expect_error(feature_table(d,spec=s,pixels_per_degree=0),"pixels_per_degree")
 enriched<-compute_contour_properties(e,spec=s,pixels_per_degree=40)
 expect_equal(average_property(enriched,"curvature"),3)
 expect_equal(sum(as.data.frame(enriched)$length),30)
 expect_equal(enriched$feature_spec$settings$pixels_per_degree,40)
 rotated<-rotate_drawing(enriched,30)
 expect_null(rotated$segment_properties)
 expect_identical(rotated$feature_spec,enriched$feature_spec)
 expect_equal(average_property(rotated,"curvature"),3)
 file<-tempfile();on.exit(unlink(file))
 saveRDS(list(enriched,model),file)
 restored<-readRDS(file)
 expect_identical(restored[[1]],enriched)
 expect_identical(restored[[2]]$schema,model$schema)
 expect_equal(predict_contours(restored[[1]],restored[[2]]),
              predict_contours(enriched,model))
})

test_that("batch and CSV exports retain input metadata and reject incompatible schemas", {
 a<-contract_drawing(contract_xy("elbow"));b<-a;b$width<-2000L
 expect_error(batch_statistics(list(a,b),properties="length"),"mismatch")
 s<-feature_spec(c("length","curvature"),estimator="arclength",curvature_scale=2)
 z<-batch_statistics(list(a=a,b=b),spec=s)
 expect_identical(names(attr(z,"measurement_metadata")$inputs),c("a","b"))
 expect_equal(z$value[z$image_id=="a"],z$value[z$image_id=="b"])
 file<-tempfile(fileext=".csv");on.exit(unlink(c(file,paste0(file,".schema.rds"))))
 expected<-write_feature_csv(a,file,spec=s)
 expect_equal(read.csv(file)$value,expected$value,tolerance=1e-12)
 expect_identical(readRDS(paste0(file,".schema.rds")),attr(expected,"measurement_metadata"))
 empty<-feature_table(mlv_drawing(list(),10,10),spec=s)
 expect_equal(nrow(empty),0)
 expect_identical(feature_schema(empty),feature_schema(feature_table(a,spec=s)))
})

test_that("MAT settings and validity propagate through the complete pipeline", {
 d<-mlv_drawing(list(matrix(c(15,5,15,35),1),matrix(c(25,5,25,35),1)),40,40)
 s<-feature_spec(c("parallelism","separation"),branch_estimator="arclength",
   branch_window=3,branch_power=2,reference_radius=2,mapping_radius=2,
   mat=list(min_component=1,min_branch_points=2))
 m<-compute_mat(d,spec=s);m<-compute_mat_properties(m)
 expect_equal(m$parameters$window,3)
 expect_equal(m$parameters$power,2)
 b<-m$branches[[1]]
 expect_equal(b$separation,score_branch(b,"separation",window=3,power=2,
   estimator="arclength",reference_radius=2))
 x<-map_mat_properties(d,m)
 expect_equal(x,analyse_mat(d,spec=s)$drawing)
 f<-feature_table(x)
 coverage<-attr(f,"measurement_metadata")$coverage
 expect_true(all(coverage$sampled==31))
 expect_true(all(coverage$coverage>=0 & coverage$coverage<=1))
 expect_true(all(is.finite(render_drawing(x,"separation"))))
 expect_true(all(is.finite(feature_density(x,"separation"))))
 changed<-s;changed$settings$branch_window<-4
 expect_error(feature_table(x,spec=changed),"mismatch")
 expect_error(map_mat_properties(d,m,spec=changed),"mismatch")
 changed<-s;changed$settings$mat$epsilon<-2
 expect_error(compute_mat_properties(m,spec=changed),"raster specification")
 shifted<-shift_contours(x,matrix(c(1,1,0,0),2))
 expect_null(shifted$mat_properties)
 expect_error(feature_table(shifted),"lack a feature specification")
 x$contours[[1]][1,1]<-14
 expect_error(feature_table(x),"stale")
 expect_error(render_drawing(x,"separation"),"stale")
 expect_error(feature_density(x,"separation"),"stale")
 expect_error(map_mat_properties(x,m),"source geometry")
 zero<-map_mat_properties(d,m)
 zero$mat_properties$parallelism<-lapply(zero$mat_properties$parallelism,function(v)rep(0,length(v)))
 zero$mat_properties$parallelism[[2]][]<-NA_real_
 z<-feature_table(zero)
 cv<-attr(z,"measurement_metadata")$coverage
 expect_equal(cv$mapped[cv$property=="parallelism"],c(31L,0L))
 expect_gt(z$parallelism_1[1],0)
 expect_equal(sum(as.numeric(z[2,grepl('^parallelism_',names(z))])),0)
})

test_that("tracing provenance and derived rendering use the selected recipe", {
 mask<-matrix(FALSE,40,40);mask[10,5:20]<-TRUE;mask[10:30,20]<-TRUE
 s<-feature_spec("curvature",estimator="arclength",curvature_scale=2,
                 extraction=list(thin=FALSE,simplify=.1))
 x<-trace_contours(mask,spec=s)
 expect_equal(x$source$parent$simplify,.1)
 expect_false(x$source$parent$thin)
 expect_identical(x$feature_spec,s)
 expect_equal(sum(feature_density(x,"curvature")),90,tolerance=1e-10)
 expect_equal(dim(render_drawing(x,"curvature")),c(40,40,3))
 image<-array(1,c(40,40,3));provider<-function(image)mask
 expect_error(trace_rgb(image,provider,spec=s),"provider")
 s$settings$extraction$provider<-"test";s$settings$extraction$provider_version<-"1"
 rgb<-trace_rgb(image,provider,spec=s)
 expect_equal(rgb$source$extraction_provider,"test")
 expect_equal(rgb$source$provider_version,"1")
 expect_error(trace_rgb(image,provider,spec=s,min_points=3),"spec or")
})
