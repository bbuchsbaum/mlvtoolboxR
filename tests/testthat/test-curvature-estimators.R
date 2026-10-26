test_that("new curvature estimators pass the frozen independent contract", {
 expect_curvature_contract(curvature_profile,contour_turning)
})
test_that("clipped tents and opposing turns conserve their stated masses", {
 clipped<-contract_drawing(rbind(c(0,0),c(1,0),c(1,-10)))
 p<-curvature_profile(clipped,scale=2)
 expect_equal(sum(p$signed_curvature*p$length),90,tolerance=1e-10)
 expect_equal(p$signed_start[1],90*2/7,tolerance=1e-10)
 expect_equal(max(p$curvature_start),90*4/7,tolerance=1e-10)
 # Two opposite 90-degree turns, separated by h=2, away from ends.
 d<-contract_drawing(rbind(c(0,0),c(10,0),c(10,-2),c(20,-2)))
 p<-curvature_profile(d,scale=2)
 expect_equal(sum(p$signed_curvature*p$length),0,tolerance=1e-10)
 expect_equal(sum(p$curvature*p$length),135,tolerance=1e-10)
 expect_equal(sum(contour_turning(d)$absolute_turn),180,tolerance=1e-10)
 # Stable clipped normalization when the requested support dwarfs the curve.
 p<-curvature_profile(clipped,scale=1e15)
 expect_equal(sum(p$signed_curvature*p$length),90,tolerance=1e-10)
})
test_that("closed seams, geometric units and straight/empty cases are explicit", {
 xy<-rbind(c(0,0),c(10,0),c(10,-10),c(0,-10))
 d<-contract_drawing(xy,TRUE)
 p<-curvature_profile(d,scale=10)
 expect_equal(p$curvature,rep(9,nrow(p)),tolerance=1e-10)
 seam<-contract_drawing(rbind(c(5,0),xy[2:4,],xy[1,]),TRUE)
 for(method in c("arclength","turning"))
  expect_equal(contract_integrals(curvature_profile(seam,method,scale=10)),
               contract_integrals(curvature_profile(d,method,scale=10)),tolerance=1e-10)
 physical<-curvature_profile(d,scale=.5,units="visual_degrees",pixels_per_degree=20)
 expect_equal(sum(physical$length),2,tolerance=1e-12)
 expect_equal(physical$curvature,rep(180,nrow(physical)),tolerance=1e-10)
 expect_equal(sum(physical$signed_curvature*physical$length),360,tolerance=1e-10)
 straight<-contract_drawing(rbind(c(0,0),c(3,0),c(10,0)))
 expect_equal(curvature_profile(straight)$curvature,0)
 expect_equal(nrow(contour_turning(straight)),0)
 expect_equal(nrow(curvature_profile(mlv_drawing(list(),10,10))),0)
 expect_error(curvature_profile(d,scale=0),"scale")
 expect_error(curvature_profile(d,units="visual_degrees"),"pixels_per_degree")
 tiny<-contract_drawing(rbind(c(0,0),c(1e-15,0)))
 expect_equal(sum(curvature_profile(tiny)$length),1e-15,tolerance=1e-30)
})
test_that("the frozen raster tracing and simplification grid meets its sensitivity bound", {
 for(res in c(1,2))for(offset in c(0,.25))for(tol in c(0,.25,.5,1)) {
  t<-seq(0,2*pi,length.out=513)
  xy<-cbind((50+offset+40*cos(t))*res,(50+offset-40*sin(t))*res);xy[513,]<-xy[1,]
  d<-mlv_drawing(list(cbind(xy[-513,],xy[-1,])),101*res,101*res)
  traced<-trace_contours(rasterize(d),simplify=tol)
  expect_length(traced$contours,1)
  z<-traced$contours[[1]]
  expect_equal(unname(z[1,1:2]),unname(tail(z,1)[,3:4]))
  p<-curvature_profile(traced,scale=8,units="visual_degrees",pixels_per_degree=res)
  expect_equal(abs(sum(p$signed_curvature*p$length)),360,tolerance=1e-7)
  mean<-sum(p$curvature*p$length)/sum(p$length)
  expect_lt(abs(mean/(180/(pi*40))-1),.25)
 }
 # Small geometric displacement alone cannot bound total absolute turning.
 zigzag<-contract_drawing(contract_xy("zigzag"))
 expect_gt(sum(contour_turning(zigzag)$absolute_turn),360)
})
