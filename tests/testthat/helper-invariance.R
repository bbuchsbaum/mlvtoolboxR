# Frozen before implementing the new estimators. Values below follow elementary
# geometry; no production curvature helper generates the references.
contract_drawing <- function(xy,closed=FALSE) {
 if(closed)xy<-rbind(xy,xy[1,])
 mlv_drawing(list(cbind(xy[-nrow(xy),,drop=FALSE],xy[-1,,drop=FALSE])),200,200)
}
contract_xy <- function(name) switch(name,
 elbow=rbind(c(0,0),c(10,0),c(10,-20)),
 elbow_split=rbind(c(0,0),c(3,0),c(10,0),c(10,-2),c(10,-20)),
 triangle=rbind(c(0,0),c(13,0),c(8,-7)),
 retrace=rbind(c(0,0),c(10,0),c(0,0)),
 zigzag=cbind(0:40,rep(c(-.1,.1),length.out=41)),
 stop("unknown contract fixture"))
contract_circle <- function(n=512,radius=10) {
 theta<-seq(0,2*pi,length.out=n+1L)[-(n+1L)]
 contract_drawing(cbind(20+radius*cos(theta),20-radius*sin(theta)),closed=TRUE)
}
contract_integrals <- function(p) c(length=sum(p$length),
 signed=sum(p$signed_curvature*p$length),absolute=sum(p$curvature*p$length),
 square=sum(p$length*(p$curvature_start^2+p$curvature_start*p$curvature_end+p$curvature_end^2)/3))
expect_curvature_contract <- function(profile,turns) {
 # API-independent gate: a candidate must supply both operations explicitly.
 stopifnot(is.function(profile),is.function(turns))
 d<-contract_drawing(contract_xy("elbow"));split<-contract_drawing(contract_xy("elbow_split"))
 for(method in c("turning","arclength")) {
  a<-profile(d,estimator=method,scale=2);b<-profile(split,estimator=method,scale=2)
  expect_equal(contract_integrals(a),contract_integrals(b),tolerance=1e-10)
  expect_equal(sum(a$length),30,tolerance=1e-12)
  expect_equal(sum(a$signed_curvature*a$length),90,tolerance=1e-10)
  expect_equal(sum(a$curvature*a$length),90,tolerance=1e-10)
  rev<-contract_drawing(contract_xy("elbow_split")[5:1,])
  z<-contract_integrals(profile(rev,estimator=method,scale=2));z['signed']<--z['signed']
  expect_equal(z,contract_integrals(a),tolerance=1e-10)
  rotated<-rotate_drawing(split,37,center=c(0,0))
  expect_equal(contract_integrals(profile(rotated,estimator=method,scale=2)),contract_integrals(a),tolerance=1e-10)
  scaled<-contract_drawing(contract_xy("elbow_split")*3)
  v<-contract_integrals(profile(scaled,estimator=method,scale=6))
  expect_equal(v,c(length=90,signed=90,absolute=90,square=contract_integrals(a)['square']/3),ignore_attr=TRUE,tolerance=1e-10)
 }
 # A single 90-degree tent of half-width 2 has height 45, area 90,
 # and squared integral 45^2 * 4 / 3 = 2700.
 p<-profile(d,estimator="arclength",scale=2)
 expect_equal(contract_integrals(p)['square'],c(square=2700),tolerance=1e-10)
 expect_equal(max(c(p$curvature_start,p$curvature_end)),45,tolerance=1e-10)
 # Unequal-arm dual cell is 15 units wide, with density 90/15 = 6.
 p<-profile(d,estimator="turning")
 expect_equal(sum(p$length[p$curvature>0]),15,tolerance=1e-10)
 expect_equal(unique(p$curvature[p$curvature>0]),6,tolerance=1e-10)
 xy<-contract_xy("triangle")
 for(method in c("turning","arclength")) {
  ref<-contract_integrals(profile(contract_drawing(xy,TRUE),estimator=method,scale=2))
  expect_equal(unname(ref['signed']),360,tolerance=1e-9)
  for(k in 1:3) {
   order<-c(k:3,if(k>1)seq_len(k-1) else integer())
   expect_equal(contract_integrals(profile(contract_drawing(xy[order,],TRUE),estimator=method,scale=2)),ref,tolerance=1e-9)
  }
 }
 circle<-profile(contract_circle(),estimator="arclength",scale=2)
 expected<-180/(pi*10)
 expect_lt(sum(abs(circle$curvature-expected)*circle$length)/sum(circle$length)/expected,.04)
 expect_equal(sum(circle$signed_curvature*circle$length),360,tolerance=1e-8)
 expect_error(profile(contract_drawing(contract_xy("retrace")),estimator="arclength"),"180|revers")
 expect_equal(sum(turns(d)$absolute_turn),90,tolerance=1e-12)
 expect_equal(nrow(turns(mlv_drawing(list(),10,10))),0)
}
