graph_fixture <- function() {
 env<-new.env(parent=environment())
 sys.source(system.file("fixtures","graph-v1.R",package="mlvtoolbox"),envir=env)
 env$graph_v1_fixtures()
}
fixture_graph <- function(x,f=graph_fixture(),...)do.call(contour_graph,
  c(list(x=x,spec=f$spec),utils::modifyList(f$settings,list(...))))

test_that("arrangement topology matches hand-enumerated frozen expectations", {
 f<-graph_fixture()
 for(i in seq_len(nrow(f$expected))) {
  expected<-f$expected[i,];x<-f$topology[[expected$case]]
  g<-fixture_graph(x,f);actual<-graph_features(g)
  expect_equal(unname(as.numeric(actual[1,names(expected)[-1]])),
               unname(as.numeric(expected[1,-1])),info=expected$case)
  expect_setequal(g$edges$contour_id,x$ids)
  expect_equal(sum(g$edges$length),sum(as.data.frame(x,spec=f$spec)$length),tolerance=1e-8)
  expect_identical(g$source_geometry,mlvtoolbox:::.geometry_snapshot(x))
 }
 bowtie<-fixture_graph(f$topology$bowtie,f)
 expect_true(bowtie$contours$closed);expect_false(bowtie$contours$simple_closed)
 expect_error(fixture_graph(f$overlap,f),"overlapping spans")
})

test_that("relationships distinguish matched-marginal organizations", {
 f<-graph_fixture()
 for(pair in f$pairs) {
  a<-contour_statistics(pair[[1]],spec=f$spec);b<-contour_statistics(pair[[2]],spec=f$spec)
  expect_equal(a$value,b$value,tolerance=1e-8)
 }
 summarise<-function(pair)lapply(f$pairs[[pair]],function(x)graph_features(fixture_graph(x,f)))
 topology<-summarise("connectivity")
 expect_equal(vapply(topology,`[[`,integer(1),"components"),c(1L,4L))
 expect_equal(vapply(topology,`[[`,integer(1),"cycles"),c(1L,0L))
 continuation<-summarise("continuation")
 expect_equal(vapply(continuation,`[[`,integer(1),"continuation_pairs"),c(1L,0L))
 parallel<-summarise("parallel")
 expect_equal(vapply(parallel,`[[`,integer(1),"parallel_contour_pairs"),c(1L,0L))
 expect_equal(parallel[[1]]$parallel_shared_extent,10)
 expect_equal(parallel[[1]]$mean_parallel_separation,5)
 expect_true(is.na(parallel[[2]]$mean_parallel_separation))
 enclosure<-summarise("enclosure")
 expect_equal(vapply(enclosure,`[[`,integer(1),"enclosures"),c(1L,0L))
 expect_equal(vapply(enclosure,`[[`,integer(1),"cycles"),c(2L,2L))
 g<-fixture_graph(f$pairs$enclosure[[1]],f)
 expect_equal(g$enclosure$outer_area,400)
 expect_equal(g$enclosure$inner_area,25)
 expect_equal(g$enclosure$clearance,5)
})

test_that("graph descriptors obey geometric and physical-unit invariances", {
 f<-graph_fixture()
 for(pair in f$pairs) {
  x<-pair[[1]];ref<-graph_features(fixture_graph(x,f))
  shifted<-shift_contours(x,matrix(rep(c(3,7),each=length(x$ids)),ncol=2))
  expect_equal(graph_features(fixture_graph(shifted,f)),ref,tolerance=1e-8)
  expect_equal(graph_features(fixture_graph(rotate_drawing(x,37),f)),ref,tolerance=1e-8)
  reversed<-mlv_drawing(lapply(x$contours,function(z)z[nrow(z):1,c(3,4,1,2),drop=FALSE]),x$width,x$height,x$ids)
  expect_equal(graph_features(fixture_graph(reversed,f)),ref,tolerance=1e-8)
  physical<-f;physical$spec<-feature_spec(c("length","curvature","orientation"),
    normalize=FALSE,units="visual_degrees",pixels_per_degree=1,estimator="arclength")
  scaled<-mlv_drawing(lapply(x$contours,function(z)z*2),200,200,x$ids)
  a<-fixture_graph(x,physical);b<-fixture_graph(scaled,physical,pixels_per_degree=2)
  expect_equal(as.matrix(graph_features(a)),as.matrix(graph_features(b)),tolerance=1e-8)
  expect_identical(attr(graph_features(a),"graph_contract"),attr(graph_features(b),"graph_contract"))
  expect_equal(attr(graph_features(b),"input_pixel_size"),.5)
 }
 x<-f$pairs$parallel[[1]];base<-graph_features(fixture_graph(x,f))
 x$contours[[1]]<-rbind(c(0,0,4,0),c(4,0,10,0))
 split<-graph_features(fixture_graph(x,f))
 keep<-setdiff(names(base),c("vertices","edges"))
 expect_equal(split[,keep],base[,keep],tolerance=1e-8)
 expect_equal(split$vertices,base$vertices+1)
 expect_equal(split$edges,base$edges+1)
})

test_that("domain boundaries, identity and stale-graph checks are explicit", {
 f<-graph_fixture();x<-f$pairs$parallel[[1]];g<-fixture_graph(x,f)
 expect_equal(nrow(as.data.frame(g)),2)
 expect_setequal(as.data.frame(g)$contour_id,x$ids)
 expect_identical(graph_features(g,x=x,spec=f$spec),graph_features(g))
 file<-tempfile();on.exit(unlink(file));saveRDS(g,file)
 expect_identical(graph_features(readRDS(file)),graph_features(g))
 expect_error(graph_features(g,x=rotate_drawing(x,10)),"stale")
 spec<-f$spec;spec$settings$curvature_scale<-2
 expect_error(graph_features(g,spec=spec),"mismatch")
 physical<-feature_spec(units="visual_degrees",pixels_per_degree=1)
 physical_graph<-contour_graph(x,spec=physical)
 physical$settings$pixels_per_degree<-2
 expect_error(graph_features(physical_graph,spec=physical),"mismatch")
 physical_graph$spec<-physical
 expect_error(graph_features(physical_graph),"specification changed")
 changed<-g;changed$parameters$tolerance<-2
 expect_error(graph_features(changed),"specification changed")
 changed<-g;changed$method<-"unknown"
 expect_error(graph_features(changed),"unsupported")
 expect_error(fixture_graph(x,f,angle_tolerance=90),"less than 90")
 expect_error(fixture_graph(x,f,tolerance=0),"tolerance")
 expect_equal(nrow(fixture_graph(x,f,parallel_distance=4.999)$parallel),0)
 expect_equal(nrow(fixture_graph(x,f,parallel_distance=5)$parallel),1)
 expect_equal(nrow(fixture_graph(f$pairs$continuation[[1]],f,gap_max=1.999)$continuation),0)
 expect_equal(nrow(fixture_graph(f$pairs$continuation[[1]],f,gap_max=2)$continuation),1)
 # A strict enclosure is not asserted when boundaries touch or are within tolerance.
 touch<-f$pairs$enclosure[[1]]
 touch$contours[[2]]<-touch$contours[[2]]-rep(c(5,0,5,0),each=4)
 expect_error(fixture_graph(touch,f),"overlapping")
 close<-f$pairs$enclosure[[1]]
 close$contours[[2]]<-sweep(close$contours[[2]],2,c(4.9995,0,4.9995,0),"-")
 expect_equal(nrow(fixture_graph(close,f)$enclosure),0)
 tiny<-mlv_drawing(list(matrix(c(0,0,.0001,0),1)),10,10)
 small<-fixture_graph(tiny,f)
 expect_equal(small$collapsed_spans,1)
 expect_equal(nrow(small$edges),0)
 expect_equal(nrow(small$nodes),1)
 many<-mlv_drawing(rep(list(matrix(c(0,0,1,0),1)),1001),10,10)
 expect_error(fixture_graph(many,f),"1000 source segments")
})

test_that("grid topology and transitive welding have independent expectations", {
 f<-graph_fixture()
 grid<-function(n)mlv_drawing(c(lapply(0:(n-1),function(k)matrix(c(0,k,n-1,k),1)),
   lapply(0:(n-1),function(k)matrix(c(k,0,k,n-1),1))),200,200)
 z<-graph_features(fixture_graph(grid(10),f))
 expect_equal(c(z$vertices,z$edges,z$components,z$cycles),c(100,180,1,81))
 expect_equal(z$branch_nodes,96)
 expect_equal(z$endpoints,0)
 expect_error(fixture_graph(grid(100),f),"20000 candidate")
 chain<-mlv_drawing(list(matrix(c(0,0,-1,-1),1),matrix(c(.0009,0,0,1),1),
  matrix(c(.0018,0,1,-1),1)),10,10)
 g<-fixture_graph(chain,f)
 expect_equal(sum(g$nodes$degree==3),1)
 expect_equal(g$nodes$x[g$nodes$degree==3],.0009,tolerance=1e-10)
 default<-contour_graph(f$topology$cross)
 expect_equal(graph_features(default)$cycles,0)
 enriched<-compute_contour_properties(f$topology$cross,spec=f$spec)
 expect_identical(contour_graph(enriched)$spec,f$spec)
})
