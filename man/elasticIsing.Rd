\name{elasticIsing}
\alias{elasticIsing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Ising model using elastic-net and cross-validation.
}
\description{
Uses the glmnet package for elastic-net computation and the cvTools package for cross-validation error. Use \code{\link{optimalGraph}} to select the optimal graph.
}
\usage{
elasticIsing(data, nLambda = 100, lambda.min.ratio = 0.01, 
          alpha = seq(0, 1, length = 10), cost = c("mspe", 
          "rmspe", "mape", "tmspe", "rtmspe"), K = 10, 
          and = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A binary dataset
}
  \item{nLambda}{
Number of lambda tuning parameters
}
  \item{lambda.min.ratio}{
Lambda min ratio, see details.
}
  \item{alpha}{
Vector with values of alpha to test
}
  \item{cost}{
Cost functions from the cvTools package to use.
}
  \item{K}{
The number of splits in k-fold cross-validation.
}
  \item{and}{
Should an AND-rule be used? If \code{TRUE}, both A -> B and B <- A need to be non-zero to obtain the egde A -- B. If \code{FALSE}, an OR-rule is used.
}
}
\details{
For each alpha, the maximum lambda is obtained from glmnet. The minimum lambda for ALL levels of alpha is lambda.min.ratio * lambda.max obtained when alpha = 1.
}
\value{
An \code{elasticIsing} object, with the following elements:
\item{minimal}{Values with minimal predictive cost}
\item{costs}{Predictive cost}
\item{lambdaMatrix}{Matrix indicating lambda values used. Columns correspond to the alpha values.}
\item{alpha}{Alpha values used}
\item{data}{Dataset used}
\item{and}{AND-rule}
}
\references{
Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for
Generalized Linear Models via Coordinate Descent. Journal of Statistical Software,
33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
  
Andreas Alfons (2012). cvTools: Cross-validation tools for regression models. R package
version 0.3.2. https://CRAN.R-project.org/package=cvTools
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{optimalGraph}}, \code{\link{costPlots}}
}
\examples{
library("IsingSampler")

# Input:
P <- 5 # Number of nodes
nSample <- 250 # Number of samples

# Chain graph:
Graph <- matrix(0, P, P)
for (i in 1:P){
  Graph[i,i\%\%P+1] <- Graph[i\%\%P+1,i] <- 0.5
}

# Thresholds:
Thresh <- rep(0, P)

# Response options (0,1 or -1,1):
Resp <- c(0L,1L)

# Simulate with metropolis:
Data <- IsingSampler(nSample, Graph, Thresh)

\dontrun{
# Estimate:
Res <- elasticIsing(Data)

# Optimal graph:
optimalGraph(Res)

# Plot result:
plot(Res)

# Cost plots:
costPlots(Res)
}
}
