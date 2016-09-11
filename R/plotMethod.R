plot.elasticIsing <- function(x,..., optArgs = list() ){
  qgraph(do.call(optimalGraph,c(list(x),optArgs))$graph, ...)
}