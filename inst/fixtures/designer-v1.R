# Frozen synthetic acceptance inputs; no optimizer is used to construct these.
designer_v1_fixtures <- function() {
 elbow<-function(length,angle,i) {
  a<-c(10+15*i,50);b<-a+c(length/2,0)
  c<-b+length/2*c(cos(angle*pi/180),-sin(angle*pi/180))
  rbind(c(a,b),c(b,c))
 }
 lengths<-c(10,4,4,4,10,4,4,4)
 angles<-c(100,90,90,90,0,1,1,1)
 contours<-Map(elbow,lengths,angles,seq_along(lengths))
 drawing<-mlv_drawing(contours,200,100,ids=sprintf("c%02d",1:8))
 spec<-feature_spec(c("length","orientation","curvature"),bins=2,normalize=FALSE,
                    estimator="arclength",curvature_scale=.1,length_range=c(0,100))
 common<-list(spec=spec,target="integrated_curvature",length_bounds=c(8,12),
  count_bounds=c(1,3),nuisance=c(length=2,orientation=12),count_imbalance=2,
  mode="set",time_limit=5,max_nodes=100000)
 case<-function(drawings,overrides=list(),expected,required_win=FALSE)
  list(drawings=drawings,arguments=utils::modifyList(common,overrides),
       expected=expected,required_win=required_win)
 small<-mlv_drawing(Map(elbow,c(5,5),c(90,0),1:2),100,100,ids=c("corner","line"))
 tiny<-lapply(1:2,function(i)subset_contours(small,i));names(tiny)<-c("a","b")
 tie<-mlv_drawing(Map(elbow,rep(5,4),rep(0,4),1:4),100,100,ids=letters[1:4])
 split<-list(first=subset_contours(drawing,c(1,2,5,6)),second=subset_contours(drawing,c(3,4,7,8)))
 split$second$width<-2000L;split$second$height<-1000L
 list(
  knapsack=case(list(image=drawing),expected=list(optimum=270,ranked=100),required_win=TRUE),
  set_knapsack=case(split,expected=list(optimum=270,ranked=100),required_win=TRUE),
  ties=case(list(tie=tie),list(length_bounds=c(5,5),count_bounds=c(1,1),
    nuisance=c(length=0,orientation=0),count_imbalance=0),expected=list(optimum=0,ranked=0)),
  small_set=case(tiny,list(length_bounds=c(5,5),count_bounds=c(1,1),
    nuisance=c(length=0,orientation=10),count_imbalance=0),expected=list(optimum=90,ranked=90)),
  per_drawing_infeasible=case(tiny,list(length_bounds=c(5,5),count_bounds=c(1,1),
    nuisance=c(length=0,orientation=10),count_imbalance=0,mode="per_drawing"),expected=list(status="infeasible")),
  smaller_budget=case(list(image=drawing),list(length_bounds=c(6,8)),expected=list(optimum=178,ranked_feasible=FALSE)),
  exhausted=case(list(image=drawing),list(max_nodes=1),expected=list(status="budget_exhausted")))
}
