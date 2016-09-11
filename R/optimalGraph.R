# Selects an optimal graph, which is recoded to the (-1, 1) parameterization and made positive definite.
optimalGraph <- function(object, cost){
  stopifnot(is(object, "elasticIsing"))
  if (missing(cost)){
    cost <- dimnames(object$costs)[[3]][1]
  }
  # Find optimals:
  opt <- which(object$costs[,,cost] == min(object$costs[,,cost]), arr.ind = TRUE)
  
  # Warning for multiple optimums:
  if (nrow(opt)>1){
    warning("Multiple optimums found, returning a list!")
  }
  
  cat("Optimal graph(s) found at:\n")
  print(data.frame(lambda = object$lambdaMatrix[opt[,1],opt[,2]], alpha = object$alpha[opt[,2]]))
  
  Res <- list()
  for (i in seq_len(nrow(opt))){
    Res[[i]] <- elasticIsingInner(object$data, object$lambdaMatrix[opt[i,1],opt[i,2]], object$alpha[opt[i,2]], and = object$and)
    Res[[i]]$thresholds <- Res[[i]]$thresholds[[1]]
    Res[[i]]$graph <- Res[[i]]$networks[[1]]
    Res[[i]] <- Res[[i]][names(Res[[i]])!="networks"]
  }
  
  # Rescale all graphs:
  Res <- lapply(Res, positiveGraph)
  
  if (length(Res)==1){
    Res <- Res[[1]]
  }
  
  return(Res)
}