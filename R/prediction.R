predictIsing <- function(net, thresh, data)
{
  # Remove diagonal:
  diag(net) <- 0
  
  # Number of observations:
  n <- nrow(data)
  k <- ncol(data)
  
  # Check:
  if (k != ncol(net))
  {
    stop("data does not have same number of variables as network")
  }
  
  # predicted values:
  Res <- matrix(0,n,k)
  
  # Start loop:
  for (i in seq_len(n))
  {
    for (j in seq_len(k))
    {
      e <- exp(net[,j] %*% data[i,] + thresh[j])
      Res[i,j] <-  e / (1 + e)
    }
  }
  
  # Return:
  return(Res)
}