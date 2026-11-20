# Summarize topology and geometric contour relations

Summarize topology and geometric contour relations

## Usage

``` r
graph_features(graph, x = NULL, spec = NULL)
```

## Arguments

- graph:

  A contour graph.

- x:

  Optional source drawing to check for stale geometry.

- spec:

  Optional requested specification for compatibility validation.

## Value

A one-row data frame with graph sizes, components, cycle rank, degree
counts, closure, continuation, parallel extent/separation and enclosure.
Mean parallel separation is NA when no parallel relations exist. Raw
vertex and edge counts depend on retained degree-two vertices; cycle
rank does not.
