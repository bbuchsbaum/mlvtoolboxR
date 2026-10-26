test_that("acceptance fixtures preserve the known legacy counterexample", {
 xy<-contract_xy("elbow")
 d<-contract_drawing(xy)
 expect_equal(average_property(d,"curvature"),6)
 # This split repeats the terminal direction, exposing the duplicated turn.
 split<-contract_drawing(rbind(xy[1:2,],c(10,-10),xy[3,]))
 expect_equal(average_property(split,"curvature"),3)
 s<-as.data.frame(d)
 expect_equal(sum(s$curvature*s$length),180)
 expect_equal(sum(as.data.frame(split)$length),sum(s$length))
 expect_equal(sum(as.data.frame(contract_circle())$length),2*pi*10,tolerance=1e-4)
 expect_error(expect_curvature_contract(NULL,NULL))
})
