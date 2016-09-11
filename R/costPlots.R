### Plotting function
costPlots <- function(object, filename = "elasticIsing.pdf", width= 8, height= 5,  theta = 25, phi = 30, ticktype = "simple",
                      accuracy = TRUE # If TRUE, plot accuracy instead
                      )
{
  # Check for akima:
  if (!"akima" %in% rownames(installed.packages())){
    stop("'akima' package not installed")
  }
  
  lambdaMatrix <- object$lambdaMatrix
  alpha <- object$alpha
    
  stopifnot(is(object, "elasticIsing"))
  pdf(filename,width=width,height=height)
  ###
  for (c in seq_len(dim(object$costs)[3]))
  {
    costName <- dimnames(object$costs)[[3]][c]
    
    Cost <- object$costs[,,c]
    if (accuracy){
      Cost <- -Cost
      lab <- "Accuracy"
    } else {
      lab <- "Cost"
    }
    rownames(Cost) <- colnames(Cost) <- NULL
    
    melted <- melt(Cost)
 
    names(melted) <- c("lambda","alpha","cost")
    # melted$lambda <- lambdavec[melted$lambda] %>% round(3) %>% as.character
    # melted$alpha <- alphavec[melted$alpha]%>% round(3) %>% as.character
    for (i in seq_len(nrow(melted))){
      melted$lambda[i] <- lambdaMatrix[melted$lambda[i],melted$alpha[i]]
    }
    melted$alpha <- alpha[melted$alpha]
    
    
 
    #   precision <- -Cost+max(Cost)
    #   precision <- precision / max(precision)
    #   precision <- round(precision*99) + 1
    cols <-  terrain.colors(100)
#       preCols <- matrix(NA,nrow(precision),ncol(precision))
#       preCols[] <- cols[c(precision)]

    # preCols <- preCols[nrow(preCols):1,]
    #   library(rgl)
    # surface3d(lambda,alpha,ElasticRes$cost, color = terrain.colors(100)[precision], xlim=c(0,1))
    # axes3d( xlab = "lambda")
    # aspect3d(1,1,0.5)
    # 
    # persp3d(lambda[order(lambda)],alpha,Cost[order(lambda),], color = preCols[order(lambda),], 
    #         xlab = "Lambda", ylab = "Alpha", zlab = "Cost", aspect = c(1,1,0.5),theta=300)
    # rgl.viewpoint(30,30)
    # min <- which(Cost == min(Cost), arr.ind = TRUE)
    # lines3d(x= lambda[min[1]], y = alpha[min[2]], z = range(Cost), lwd = 3, col = "red")
    
#     persp(lambda[order(lambda)],alpha,Cost[order(lambda),], col = preCols, 
#           expand = 0.5 , shade = 0, border = NA,
#           xlab = "Lambda", ylab = "Alpha", zlab = "Cost", theta = 130, phi = 30,
#           main = costName)
#     
    # On 100 by 100 grid, compute tiles:
    # interpolation
    interpolated <- akima::interp(x = log(melted$lambda), y = melted$alpha, z = melted$cost,
                                  xo = seq(log(min(melted$lambda)),log(max(melted$lambda)),length=100),
                                  yo = sort(unique(melted$alpha)))

    # Fill in NAs of strictly regularized models with lowest accuracy or highest cost.
    if (accuracy){
      interpolated$z[is.na(interpolated$z)] <- min(interpolated$z[!is.na(interpolated$z)])
    } else {
      interpolated$z[is.na(interpolated$z)] <- max(interpolated$z[!is.na(interpolated$z)])
    }

    
    # Colors:
    z <- interpolated$z[order(interpolated$x),]
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    
    preCols <- cols[cut(zfacet,100)]
    
    persp(interpolated$x,interpolated$y,interpolated$z, col = preCols, 
          expand = 0.5 , shade = 0, border = NA,
          xlab = "Lambda (log)", ylab = "Alpha", zlab = lab, theta = theta, phi = phi,
          main = costName, ticktype =ticktype) -> res
#     
#     persp(lambda[order(lambda)],alpha,Cost[order(lambda),], col = preCols, 
#           expand = 0.5 , shade = 0, border = NA,
#           xlab = "Lambda", ylab = "Alpha", zlab = lab, theta = theta, phi = phi,
#           main = costName, ticktype =ticktype) -> res

  }
  dev.off()
  message(paste0("Output stored in ",getwd(),"/",filename))
}
