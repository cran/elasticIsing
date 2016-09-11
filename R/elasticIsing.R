elasticIsing <- function(
  data, # Binary data
  nLambda = 100,
  lambda.min.ratio = 0.01, # min lambda is selected by taking min.ratio * lambda_max(alpha=1)
  # lambda, # Vector of penalty lambdas to test (can be missing)
  alpha = seq(0,1,length=10), # Vector of alphas to test.
  cost = c('mspe','rmspe','mape','tmspe','rtmspe'), # Cost functions to apply
  K = 10, # number of folds
  and = TRUE
){
  data <- as.matrix(data)
  
  # Parameters:
  Np <- nrow(data) # Number of persons
  Ni <- ncol(data) # Number of items
  Nc <- length(cost) # Number of cost functions
  
  # Sort alpha 1 to 0:
  alpha <- sort(alpha,decreasing = TRUE)
  
  # Lengths:
  Nl <- nLambda # Number of lambdas
  Na <- length(alpha) # Number of alphas
  
  # Lambda matrix (to fill):
  lambdaMatrix <- matrix(NA, Nl, Na)
  
  # Create folds:
  Folds <- cvTools::cvFolds(Np, K)
  
  # Test if data is binary, stop if not:
  if (!all(sort(unique(c(unlist(data)))) == c(0,1))){
    stop("Binary data required")
  }
  
  # Matrix to store prediction cost:
  predictionCost <- array(NA, c(Nl, Na,Nc))
  dimnames(predictionCost) <- list(NULL,alpha,cost)
  
  # Initialize progress bar:
  pb <- utils::txtProgressBar(min = 0, max = Na * K, initial = 0, style = 3)
  
  
  # For every alpha:
  for (a in seq_len(Na)){
    # Compute lambdas:
    maxLambda <- max(sapply(seq_len(Ni),function(i){
      Res <- glmnet(data[,-i,drop=FALSE],data[,i],family = "binomial", alpha = alpha[a],nlambda = 4)
      max(Res$lambda)
      }))
    
    if (a==1){
      # Compute min lambda:
      minLambda <- maxLambda * lambda.min.ratio
    }
    
    # Lambda sequence:
    lambda <- 10^seq(log10(maxLambda), log10(minLambda), len=Nl)
      
    # Make lambda decreasing:
    lambda <- sort(lambda, decreasing = TRUE)
    
    # Fill matrix:
    lambdaMatrix[,a] <- lambda
    
    # List with predicted values for every lambda:
    PredictedValues <- rep(list(matrix(NA,Np, Ni)), Nl)
    
    # For every fold:
    for (k in seq_len(K)){
      # Which in block:
      inBlock <- sort(Folds$subsets[Folds$which==k,1])
      
      # Estimte Ising:
      foldEstimates <- elasticIsingInner(data[-inBlock,], lambda, alpha[a],   and)
      
      # For every lambda, predict and store:
      for (l in seq_len(Nl)){
        PredictedValues[[l]][inBlock,] <- predictIsing(foldEstimates$networks[[l]], foldEstimates$thresholds[[l]], data[inBlock,])
      }
      
      utils::setTxtProgressBar(pb, (a-1)*K + k)
    }
    
    
    # Compute prediction cost:
    for (c in seq_len(Nc)){
      predictionCost[,a,c] <- sapply(PredictedValues, function(p) do.call(cost[c],list(data, p))) 
    }
  }
  close(pb)
  
  # Select values with minimal predictive cost:
  Mininals <- list()
  for (c in seq_len(Nc)){
    Mininals[[c]] <- which(predictionCost[,,c] == min(predictionCost[,,c]), arr.ind = TRUE)
    Mininals[[c]] <- as.data.frame(Mininals[[c]])
    Mininals[[c]][,1] <- lambda[Mininals[[c]][,1]]
    Mininals[[c]][,2] <- alpha[Mininals[[c]][,2]]
    colnames(Mininals[[c]]) <- c("lambda","alpha")
    rownames(Mininals[[c]]) <- NULL
  }
  
  #   
  #   cat("Minimal prediction costs found at:\n")
  #   print(Minimal)
  #   
  Res <- list(
    minimal = Mininals,
    costs = predictionCost,
    lambdaMatrix = lambdaMatrix,
    alpha = alpha,
    data = data,
    and = and)
  class(Res) <- "elasticIsing"
  return(Res)
}


## Inner function to compute elastic net Ising model for series lambda and single alpha:
elasticIsingInner <- function(
  data, # Binary data
  lambda, # Vector of penalty lambdas to test
  alpha, # Vector of alphas to test.
  and = TRUE
){
  data <- as.matrix(data)
  # Parameters:
  Np <- nrow(data)
  Ni <- ncol(data)
  Nl <- length(lambda)
  stopifnot(length(alpha)==1)
  
  # Run glm for every node:
  ResNodes <- list()
  
  for (n in seq_len(Ni)){
    ResNodes[[n]] <- glmnet(as.matrix(data[,-n]),data[,n],"binomial",alpha = alpha, lambda = lambda)
  }
  
  # Construct networks for every lambda:
  Coefs <- lapply(ResNodes, coef)
  Nets <- rep(list(matrix(0,Ni,Ni)),Nl)
  Thresh <- rep(list(vector("numeric", Ni)), Nl)
  
  for (i in seq_len(Ni)){
    for (l in seq_len(Nl)){
      Nets[[l]][-i,i] <- Coefs[[i]][-1,l]
      Thresh[[l]][i] <- Coefs[[i]][1,l]
    }
  }
  
  # Symmetrize:
  if (and){
    Nets <- lapply(Nets, function(x) (x!=0 & t(x)!= 0) * (x + t(x))/2  )
  } else {
    Nets <- lapply(Nets, function(x) (x!=0 | t(x)!= 0) * (x + t(x))/2  )
  }
  
  # Return:
  return(list(networks = Nets, thresholds = Thresh))
}


