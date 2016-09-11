# Rescales a function to (-1, 1) parameterization and makes the graph positive definite
positiveGraph <- function(x){
  resRescaled <- IsingSampler::LinTransform(x$graph, x$thresholds, c(0,1), c(-1,1))
  resRescaled$graph <- resRescaled$graph - (diag(nrow(resRescaled$graph)) * (min(eigen(resRescaled$graph)$values)-.Machine$double.neg.eps))
  return(resRescaled)
}