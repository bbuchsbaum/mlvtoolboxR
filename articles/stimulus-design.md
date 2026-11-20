# Matched conditions across a stimulus set

Matching across a stimulus set can be feasible when each drawing is too
small to split. This minimal example has a corner in one drawing and a
straight contour in another, both ten pixels long.

``` r

drawings <- list(
  corner = mlv_drawing(list(rbind(c(5,5,10,5), c(10,5,10,10))),20,20),
  straight = mlv_drawing(list(matrix(c(5,5,15,5),1)),20,20))
spec <- feature_spec(c("length","orientation","curvature"), bins=2,
                     normalize=FALSE, estimator="arclength", curvature_scale=.5)
result <- design_subsets(drawings,spec,length_bounds=c(10,10),
  count_bounds=c(1,1),nuisance=c(length=0,orientation=10))
result
#> Whole-contour subset design: optimal 
#> Target difference: 90  | upper bound: 90 
#> Search nodes: 7  | elapsed: 0 seconds
as.data.frame(result)
#>   image_id contour_id condition
#> 1   corner          1         A
#> 2 straight          1         B
result$diagnostics$groups
#> [[1]]
#> [[1]]$group
#> [1] "set"
#> 
#> [[1]]$feasible
#> [1] TRUE
#> 
#> [[1]]$length_A
#> [1] 10
#> 
#> [[1]]$length_B
#> [1] 10
#> 
#> [[1]]$count_A
#> [1] 1
#> 
#> [[1]]$count_B
#> [1] 1
#> 
#> [[1]]$objective
#> [1] 90
#> 
#> [[1]]$residual
#>      length orientation 
#>           0          10
result$diagnostics$per_drawing
#> [[1]]
#> [[1]]$group
#> [1] "corner"
#> 
#> [[1]]$feasible
#> [1] FALSE
#> 
#> [[1]]$length_A
#> [1] 10
#> 
#> [[1]]$length_B
#> [1] 0
#> 
#> [[1]]$count_A
#> [1] 1
#> 
#> [[1]]$count_B
#> [1] 0
#> 
#> [[1]]$objective
#> [1] 90
#> 
#> [[1]]$residual
#>      length orientation 
#>          10          10 
#> 
#> 
#> [[2]]
#> [[2]]$group
#> [1] "straight"
#> 
#> [[2]]$feasible
#> [1] FALSE
#> 
#> [[2]]$length_A
#> [1] 0
#> 
#> [[2]]$length_B
#> [1] 10
#> 
#> [[2]]$count_A
#> [1] 0
#> 
#> [[2]]$count_B
#> [1] 1
#> 
#> [[2]]$objective
#> [1] 0
#> 
#> [[2]]$residual
#>      length orientation 
#>          10          10
```

The optimizer maximizes the difference in integrated absolute curvature
at the selected scale. Orientation matching uses raw length-weighted bin
masses: the tolerated L1 distance here is ten pixels. The resulting
target difference is 90 geometric degrees. This is a geometric
construction, not evidence of a perceptual effect. A and B are disjoint
and contain only original whole contours.

``` r

within_drawing <- design_subsets(drawings,spec,length_bounds=c(10,10),
  count_bounds=c(1,1),nuisance=c(length=0,orientation=10),mode="per_drawing")
within_drawing$status
#> [1] "infeasible"
```

Each source has only one contour, so it cannot provide two nonempty
conditions. The result reports proven infeasibility. A time or node
limit has a different meaning:

``` r

limited <- design_subsets(drawings,spec,length_bounds=c(10,10),
  count_bounds=c(1,1),nuisance=c(length=0,orientation=10),time_limit=0)
limited$status
#> [1] "budget_exhausted"
limited$solver[c("upper_bound","nodes","search_seconds")]
#> $upper_bound
#> [1] 90
#> 
#> $nodes
#> [1] 0
#> 
#> $search_seconds
#> [1] 0
```

Zero time performs no search. `budget_exhausted` makes no feasibility
claim. If a limited run does find a solution, `feasible` retains its
verified residuals and a conservative objective bound. Only completed
searches certify optimality or infeasibility, within the declared
numerical tolerance.

``` r

path <- tempfile(fileext=".rds")
saveRDS(result,path)
identical(readRDS(path)$selection,result$selection)
#> [1] TRUE
unlink(path)
```

The result preserves input geometry, provenance, source/contour IDs,
specification and reconstructed condition drawings. Physical
measurements use visual degrees with explicit input calibration.
MAT/junction features can be requested as fresh diagnostics after
selection; they are not additive optimization constraints.
